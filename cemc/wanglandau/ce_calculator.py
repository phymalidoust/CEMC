from ase.calculators.calculator import Calculator
from ase.ce.corrFunc import CorrFunction
from ase.ce import BulkCrystal
from ase.ce import BulkSpacegroup
from ase.build import bulk
import unittest
from itertools import product, combinations
import os
import numpy as np
import copy
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
from matplotlib import pyplot as plt
from ase.visualize import view
from cemc.mcmc import linear_vib_correction as lvc
from mpi4py import MPI
try:
    from cemc.ce_updater import ce_updater as ce_updater
    use_cpp = True
except Exception as exc:
    use_cpp = False
    print (str(exc))
    print ("Could not find C++ version, falling back to Python version")

def get_ce_calc( small_bc, bc_kwargs, eci=None, size=[1,1,1], free_unused_arrays_BC=False ):
    """
    Constructs a CE calculator by first computing the correlation function
    from a small cell

    Arguments
    -----------
    small_bc - Instance of BulkCrystal or BulkSpacegroup with a relatively small unitcell
    bc_kwargs - dictionary of the keyword arguments used to construct small_bc
    eci - Effective Cluster Interactions
    size - The atoms in small_bc will be extended by this amount
    """
    rank = MPI.COMM_WORLD.Get_rank()
    unknown_type = False
    large_bc = small_bc # Just for the other processes
    init_cf = None
    if ( rank == 0 ):
        calc1 = CE( small_bc, eci )
        init_cf = calc1.get_cf()
        cell_lenghts = small_bc.atoms.get_cell_lengths_and_angles()[:3]
        min_length = np.min(cell_lenghts)/2.0

        bc_kwargs["size"] = size
        bc_kwargs["max_cluster_dia"] = min_length

        db_name = "temporary_db.db"
        if ( os.path.exists(db_name) ):
            os.remove( db_name )
        bc_kwargs["db_name"] = db_name

        if ( isinstance(small_bc,BulkCrystal) ):
            large_bc = BulkCrystal(**bc_kwargs)
        elif ( isinstance(small_bc,BulkSpacegroup) ):
            large_bc = BulkSpacegroup(**bc_kwargs)
        else:
            unknown_type = True

    unknown_type = MPI.COMM_WORLD.bcast( unknown_type, root=0 )
    if ( unknown_type ):
        raise TypeError( "The small_bc argument has to by of type BulkCrystal or BulkSpacegroup" )
    large_bc = MPI.COMM_WORLD.bcast( large_bc, root=0 )
    init_cf = MPI.COMM_WORLD.bcast( init_cf, root=0 )
    calc2 = CE( large_bc, eci, initial_cf=init_cf, free_unused_arrays_BC=free_unused_arrays_BC )
    return calc2

class CE( Calculator ):
    """
    Class for updating the CE when symbols change
    """

    implemented_properties = ["energy"]
    def __init__( self, BC, eci=None, initial_cf=None, free_unused_arrays_BC=False ):
        Calculator.__init__( self )
        self.BC = BC
        self.corrFunc = CorrFunction(self.BC)
        if ( initial_cf is None ):
            self.cf = self.corrFunc.get_cf( self.BC.atoms )
        else:
            self.cf = initial_cf

        if ( eci is None ):
            eci = {name:1.0 for name in self.cf.keys()}
        self.eci = eci
        # Make supercell
        self.atoms = self.BC.atoms
        symbols = [atom.symbol for atom in self.BC.atoms] # Keep a copy of the original symbols

        # Make sure that the database information fits
        if ( len(self.BC.atoms) != self.BC.trans_matrix.shape[0] ):
            msg = "The number of atoms and the dimension of the translation matrix is inconsistent"
            msg = "Dimension of translation matrix: {}. Number of atoms: {}".format(self.BC.trans_matrix.shape,len(self.BC.atoms))
            raise ValueError( msg )

        #print (self.basis_elements)
        #if ( len(BC.basis_elements) > 1 ):
        #    raise ValueError( "At the moment only one site type is supported!" )
        self.old_cfs = []
        self.old_atoms = self.atoms.copy()
        self.changes = []
        self.ctype = {}
        #self.create_ctype_lookup()
        self.convert_cluster_indx_to_list()
        self.permutations = {}
        self.create_permutations()
        self.BC.trans_matrix = np.array(self.BC.trans_matrix).astype(np.int32)
        self.updater = None
        if ( use_cpp ):
            self.updater = ce_updater.CEUpdater()
            self.updater.init( self.BC, self.cf, self.eci, self.permutations )

            if ( not self.updater.ok() ):
                raise RuntimeError( "Could not initialize C++ CE updater" )

            if ( free_unused_arrays_BC ):
                self.BC.dist_matrix = None
                self.BC.trans_matrix = None
                msg = "ClusterExpansion setting is potentially invalid.\n"
                msg += "dist_matrix and trans_matrix is set to None to save memory"
                print (msg)

        if ( use_cpp ):
            self.clear_history = self.updater.clear_history
            self.undo_changes = self.updater.undo_changes
            self.update_cf = self.updater.update_cf
        else:
            self.clear_history = self.clear_history_pure_python
            self.undo_changes = self.undo_changes_pure_python
            self.update_cf = self.update_cf_pure_python

        # Set the symbols back to their original value
        self.set_symbols(symbols)
        self._linear_vib_correction = None

    def get_full_cluster_names( self, cnames ):
        """
        Returns the full cluster names with decoration info in the end
        """
        full_names = self.cf.keys()
        print (full_names)
        only_prefix = [name.rpartition("_")[0] for name in full_names]
        full_cluster_names = []

        # First insert the one body names, nothing to be done for them
        for name in cnames:
            if ( name.startswith("c1") ):
                full_cluster_names.append(name)
        for name in cnames:
            if ( name.startswith("c1") ):
                continue
            indx = only_prefix.index(name)
            full_cluster_names.append( full_names[indx] )
        return full_cluster_names

    @property
    def linear_vib_correction( self ):
        return self._linear_vib_correction

    @linear_vib_correction.setter
    def linear_vib_correction( self, linvib ):
        if ( not isinstance(linvib,lvc.LinearVibCorrection) ):
            raise TypeError( "Linear vib correction has to be of type LinearVibCorrection!" )
        if ( self.linear_vib_correction is not None ):
            orig_eci = self.linear_vib_correction.reset(self.eci)
            if ( orig_eci is not None ):
                self.eci = orig_eci
            self.update_ecis(self.eci)
        self._linear_vib_correction = linvib
        if ( self.updater is not None ):
            # This just initialize a LinearVibCorrection object, it does not change the ECIs
            self.updater.add_linear_vib_correction( ce_updater.map_str_dbl(linvib.eci_per_kbT) )

    def include_linvib_in_ecis( self, T ):
        """
        Includes the effect of linear vibration correction in the ECIs
        """
        if ( self.linear_vib_correction is None ):
            return
        orig_eci = self.linear_vib_correction.reset( self.eci )

        # Reset the ECIs to the original
        self.update_ecis(self.eci)
        self.ecis = self.linear_vib_correction.include( self.eci, T )
        self.update_ecis(self.eci)

    def vib_energy( self, T ):
        """
        Returns the vibration energy per atom
        """
        if ( self.updater is not None ):
            return self.updater.vib_energy(T)

        if ( self.linear_vib_correction is not None ):
            return self.linear_vib_correction.energy(T,self.cf)
        return 0.0

    def initialize_correlation_functions( self ):
        """
        Initialize the correlation functions by characterizing a 4x4x4 structure
        """
        temp_db_name = "temporary_database{}.db".format(MPI.COMM_WORLD.Get_rank())
        conc_args = {
            "conc_ratio_min_1":[[0 for i in range(len(self.BC.basis_elements[0]))]],
            "conc_ratio_max_1":[[0 for i in range(len(self.BC.basis_elements[0]))]]
        }
        conc_args["conc_ratio_min_1"][0][0] = 1
        conc_args["conc_ratio_max_1"][0][-1] = 1
        clat = None
        bc = BulkCrystal( crystalstructure=self.BC.crystalstructure, a=self.BC.a, c=self.BC.c,
        size=[4,4,4], basis_elements=self.BC.basis_elements, conc_args=conc_args, db_name=temp_db_name,
        max_cluster_size=4)
        bc._get_cluster_information()
        cf = CorrFunction(bc)

        # TODO: This only works for one site type
        for atom in bc.atoms:
            atom.symbol = bc.basis_elements[0][0]

        for atom in self.BC.atoms:
            atom.symbol = bc.basis_elements[0][0]
        corr_funcs = cf.get_cf(bc.atoms)
        os.remove(temp_db_name)
        return corr_funcs

    def convert_cluster_indx_to_list( self ):
        """
        Converts potentials arrays to lists
        """
        for symm in range(len(self.BC.cluster_indx)):
            for i in range(len(self.BC.cluster_indx[symm])):
                if ( self.BC.cluster_indx[symm][i] is None ):
                    continue
                for j in range(len(self.BC.cluster_indx[symm][i])):
                    if ( self.BC.cluster_indx[symm][i][j] is None ):
                        continue
                    for k in range(len(self.BC.cluster_indx[symm][i][j])):
                        if ( isinstance(self.BC.cluster_indx[symm][i][j][k],np.ndarray) ):
                            self.BC.cluster_indx[symm][i][j][k] = self.BC.cluster_indx[symm][i][j][k].tolist()
                        else:
                            self.BC.cluster_indx[symm][i][j][k] = list(self.BC.cluster_indx[symm][i][j][k])

                    if ( isinstance(self.BC.cluster_indx[symm][i][j],np.ndarray) ):
                        self.BC.cluster_indx[symm][i][j] = self.BC.cluster_indx[symm][i][j].tolist()
                    else:
                        self.BC.cluster_indx[symm][i][j] = list(self.BC.cluster_indx[symm][i][j])

                if ( isinstance(self.BC.cluster_indx[symm][i],np.ndarray) ):
                    self.BC.cluster_indx[symm][i] = self.BC.cluster_indx[symm][i].tolist()
                else:
                    self.BC.cluster_indx[symm][i] = list(self.BC.cluster_indx[symm][i])

    def create_permutations( self ):
        """
        Creates a list of permutations of basis functions that should be passed
        to the C++ module
        """
        bf_list = list(range(len(self.BC.basis_functions)))
        for num in range(2,len(self.BC.cluster_names)):
            perm = list(product(bf_list, repeat=num))
            #perm = list(combinations(bf_list, repeat=num))
            self.permutations[num] = perm


    def get_energy( self ):
        """
        Returns the energy of the system
        """
        energy = 0.0
        for key,value in self.eci.iteritems():
            energy += value*self.cf[key]
        return energy*len(self.atoms)

    def create_ctype_lookup( self ):
        """
        Creates a lookup table for cluster types based on the prefix
        """
        for n in range(2,len(self.BC.cluster_names)):
            for ctype in range(len(self.BC.cluster_names[n])):
                name = self.BC.cluster_names[n][ctype]
                prefix = name#name.rpartition('_')[0]
                self.ctype[prefix] = (n,ctype)

    def update_cf_pure_python( self, single_change ):
        """
        Changing one element and update the correlation functions
        """
        indx = single_change[0]
        old_symb = single_change[1]
        new_symb = single_change[2]
        self.old_cfs.append( copy.deepcopy(self.cf) )
        if ( old_symb == new_symb ):
            return self.cf
        natoms = len(self.atoms)
        bf_list = list(range(len(self.BC.basis_functions)))

        self.atoms[indx].symbol = new_symb

        bf = self.BC.basis_functions
        for name in self.eci.keys():
            if ( name == "c0" ):
                continue
            elif ( name.startswith("c1") ):
                dec = int(name[-1]) - 1
                self.cf[name] += (bf[dec][new_symb]-bf[dec][old_symb])/natoms
                continue
            prefix = name.rpartition('_')[0]
            dec = int(name.rpartition('_')[-1]) - 1

            res = self.ctype[prefix]
            num = res[0]
            ctype = res[1]
            #for n in range(2, len(self.BC.cluster_names)):
            #    try:
            #        ctype = self.BC.cluster_names[n].index(prefix)
            #        num = n
            #        break
            #    except ValueError:
            #        continue
            perm = list(product(bf_list, repeat=num))
            count = len(self.BC.cluster_indx[num][ctype])*natoms
            sp = self.spin_product_one_atom( indx, self.BC.cluster_indx[num][ctype], perm[dec] )
            sp /= count
            bf_indx = perm[dec][0]
            self.cf[name] += num*( bf[bf_indx][new_symb] - bf[bf_indx][old_symb] )*sp
        return self.cf

    def spin_product_one_atom( self, ref_indx, indx_list, dec ):
        """
        Spin product for a single atom
        """
        num_indx = len(indx_list)
        bf = self.BC.basis_functions
        sp = 0.0
        for i in range(num_indx):
            sp_temp = 1.0
            for j, indx in enumerate(indx_list[i][:]):
                trans_indx = self.corrFunc.trans_matrix[ref_indx, indx]
                sp_temp *= bf[dec[j+1]][self.atoms[trans_indx].symbol]
            sp += sp_temp
        return sp

    def undo_changes_pure_python( self ):
        """
        This function undo all changes stored in all symbols starting from the
        last one
        """
        for i in range(len(self.changes),0,-1):
            entry = self.changes[i-1]
            self.atoms[entry[0]].symbol = entry[1]
            self.cf = self.old_cfs[i-1]
        self.clear_history()

    def clear_history_pure_python( self ):
        """
        Clears the history of the calculator
        """
        self.changes = []
        self.old_cfs = []

    def calculate( self, atoms, properties, system_changes ):
        """
        Calculates the energy. The system_changes is assumed to be a list
        of tuples of the form (indx,old_symb,new_symb)
        """
        if ( use_cpp ):
            energy = self.updater.calculate(system_changes)
            self.cf = self.updater.get_cf()
        else:
            self.changes += system_changes
            for entry in system_changes:
                self.update_cf( entry )
            energy = self.get_energy()
        self.results["energy"] = energy
        return self.results["energy"]

    def get_cf( self ):
        """
        Returns the correlation functions
        """
        if ( self.updater is None ):
            return self.cf
        else:
            return self.updater.get_cf()

    def update_ecis( self, new_ecis ):
        """
        Updates the ecis
        """
        self.eci = new_ecis
        if ( not self.updater is None ):
            self.updater.set_ecis(self.eci)

    def get_singlets( self, array ):
        if ( self.updater is None ):
            indx = 0
            for key,value in self.cf.iteritems():
                if ( key.startswith("c1") ):
                    singlets[indx] = value
                    indx += 1
            return array
        else:
            self.updater.get_singlets( array )
            return array

    def set_composition( self, comp ):
        """
        Change composition of an object.
        """
        # Verify that the sum of the compositions is one
        tot_conc = 0.0
        for key,conc in comp.iteritems():
            tot_conc += conc

        if ( np.abs(tot_conc-1.0) > 1E-6 ):
            raise ValueError( "The specified concentration does not sum to 1!" )
        # Change all atoms to the first one
        init_elm = comp.keys()[0]
        for i in range( len(self.atoms) ):
            #self.update_cf( (i,self.atoms[i].symbol,init_elm) ) # Set all atoms to init element
            self.calculate( self.atoms, ["energy"], [(i,self.atoms[i].symbol,init_elm)] )
        start = 0
        for elm,conc in comp.iteritems():
            if ( elm == init_elm ):
                continue
            n_at = int( round(conc*len(self.atoms)) )
            for i in range(start,start+n_at):
                self.update_cf( (i,init_elm,elm) )
            start += n_at
        self.clear_history()
        print ("Composition set")

    def set_symbols( self, symbs ):
        """
        Change the symbols of the entire atoms object
        """
        if ( len(symbs) != len(self.atoms ) ):
            raise ValueError( "Length of the symbols array has to match the length of the atoms object.!" )
        for i,symb in enumerate(symbs):
            self.update_cf( (i,self.atoms[i].symbol,symb) )
        self.clear_history()

    def write( self, fname ):
        """
        Stores all nessecary information required to restart the calculation
        from the state it ended
        """
        backup_data = {}
        backup_data["cf"] = self.get_cf()

################################################################################
##                           UNIT TESTS                                       ##
################################################################################
