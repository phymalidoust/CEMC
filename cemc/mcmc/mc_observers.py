import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import copy
import numpy as np
from ase.io.trajectory import TrajectoryWriter
from cemc.ce_updater import ce_updater

class MCObserver( object ):
    def __init__( self ):
        self.name = "GenericObserver"

    def __call__( self, system_changes ):
        """
        Gets information about the system changes and can perform some action
        """
        pass

    def reset(self):
        """
        Resets all values of the MC observer
        """
        pass

class CorrelationFunctionTracker( MCObserver ):
    """
    Class that tracks the history of the Correlation function
    Only relevant if the calculator is a CE calculator
    """
    def __init__( self, ce_calc ):
        self.cf = []
        self.ce_calc = ce_calc
        self.name = "CorrelationFunctionTracker"

    def __call__( self, system_changes ):
        """
        Updates the correlation functions
        """
        self.cf.append( copy.deepcopy(self.ce_calc.cf) )

    def plot_history( self, max_size=10 ):
        """
        Creates a plot of the history (only if history is tracked)
        """
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for key in self.cf[0].keys():
            size = int(key[1])
            if ( size > max_size ):
                continue
            cf_history = [cf[key] for cf in self.cf]
            ax.plot( cf_history, label=key, ls="steps" )
        ax.legend( loc="best", frameon=False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        return fig

class PairCorrelationObserver( MCObserver ):
    """
    Class that computes the average value of all the ECIs
    """
    def __init__( self, ce_calc ):
        self.cf = {}
        self.cf_squared = {}
        self.ce_calc = ce_calc
        if ( self.ce_calc.updater is None ):
            raise RuntimeError( "This observer can only be used with the C++ version of the CF updater" )
        self.n_entries = 0
        self.name = "PairCorrelationObserver"

        for key,value in self.ce_calc.eci.iteritems():
            if ( key.startswith("c2_") ):
                self.cf[key] = 0.0
                self.cf_squared[key] = 0.0

    def __call__( self, system_changes ):
        """
        Updates the correlation functions
        """
        new_cf = self.ce_calc.updater.get_cf()
        self.n_entries += 1
        for key in self.cf.keys():
            self.cf[key] += new_cf[key]
            self.cf_squared[key] += new_cf[key]**2

    def get_average( self ):
        """
        Returns the average
        """
        avg_cf = copy.deepcopy(self.cf)
        for key in avg_cf.keys():
            avg_cf[key] /= self.n_entries
        return avg_cf

    def get_std( self ):
        """
        Returns the standard deviation
        """
        std_cf = {key:0.0 for key in self.cf.keys()}
        for key in self.cf.keys():
            std_cf[key] = np.sqrt( self.cf_squared[key]/self.n_entries - (self.cf[key]/self.n_entries)**2 )#/np.sqrt(self.n_entries)
        return std_cf

class LowestEnergyStructure(MCObserver):
    def __init__( self, ce_calc, mc_obj ):
        self.ce_calc = ce_calc
        self.mc_obj = mc_obj
        self.lowest_energy = np.inf
        self.lowest_energy_atoms = None
        self.lowest_energy_cf = None
        self.atoms = None
        self.name = "LowestEnergyStructure"

    def __call__( self, system_changes ):
        if ( self.lowest_energy_atoms is None or self.lowest_energy_cf is None ):
            self.lowest_energy_cf = self.ce_calc.get_cf()
            self.lowest_energy_atoms = self.ce_calc.atoms.copy()
            self.lowest_energy = self.mc_obj.current_energy
            self.atoms = self.mc_obj.atoms.copy()
            return

        if ( self.mc_obj.current_energy < self.lowest_energy ):
            self.lowest_energy = self.mc_obj.current_energy
            self.lowest_energy_atoms = self.ce_calc.atoms.copy()
            self.lowest_energy_cf = self.ce_calc.get_cf()

class SGCObserver(MCObserver):
    def __init__( self, ce_calc, mc_obj, n_singlets ):
        super(SGCObserver,self).__init__()
        self.name = "SGCObersver"
        self.ce_calc = ce_calc
        self.mc = mc_obj

        self.quantities = {
            "singlets":np.zeros( n_singlets, dtype=np.float64 ),
            "singlets_sq":np.zeros( n_singlets, dtype=np.float64 ),
            "energy":0.0,
            "energy_sq":0.0,
            "singl_eng":np.zeros( n_singlets, dtype=np.float64 ),
            "counter":0
        }

        """
        # Track average value of the singlet terms
        self.singlets = np.zeros( n_singlets, dtype=np.float64 )

        # Track average value of the energy
        self.energy = 0.0

        # Track average value of energy squared
        self.energy_sq = 0.0

        # Track average value of particle-energy correlation
        self.singl_eng = np.zeros_like( self.singlets )
        self.counter = 0
        """

    def reset(self):
        """
        Resets all variables to zero
        """
        self.quantities["singlets"][:] = 0.0
        self.quantities["singlets_sq"][:] = 0.0
        self.quantities["energy"] = 0.0
        self.quantities["energy_sq"] = 0.0
        self.quantities["singl_eng"][:] = 0.0
        self.quantities["counter"] = 0
        """
        self.singlets[:] = 0.0
        self.energy = 0.0
        self.energy_sq = 0.0
        self.singl_eng[:] = 0.0
        self.counter = 0
        """

    def __call__( self, system_changes ):
        self.quantities["counter"] += 1
        new_singlets = np.zeros_like( self.singlets )
        self.ce_calc.get_singlets(  new_singlets )

        self.quantities["singlets"] += new_singlets
        self.quantities["singlets_sq"] += new_singlets**2
        self.quantities["energy"] += self.mc.current_energy_without_vib()
        self.quantities["energy_sq"] += self.mc.current_energy_without_vib()**2
        self.quantities["singl_eng"] += new_singlets*self.mc.current_energy_without_vib()
        """
        self.singlets += new_singlets
        self.energy += self.mc.current_energy
        self.energy_sq += self.mc.current_energy**2
        self.singl_eng += new_singlets*self.mc.current_energy
        """

    @property
    def energy(self):
        return self.quantities["energy"]

    @property
    def energy_sq(self):
        return self.quantities["energy_sq"]

    @property
    def singlets(self):
        return self.quantities["singlets"]

    @property
    def singl_eng(self):
        return self.quantities["singl_eng"]

    @property
    def counter(self):
        return self.quantities["counter"]

class Snapshot( MCObserver ):
    def __init__(self, trajfile="default.traj", atoms=None ):
        super(Snapshot,self).__init__()
        self.name = "Snapshot"
        if ( not trajfile.endswith(".traj") ):
            raise ValueError( "This object stores all images in a trajectory file. File extension should be .traj" )
        if ( atoms is None ):
            raise ValueError( "No atoms object given!" )
        self.atoms = atoms
        self.traj = TrajectoryWriter( trajfile, mode="a" )


    def __call__( self, system_changes ):
        self.traj.write(self.atoms)

class NetworkObserver( MCObserver ):
    def __init__( self, calc=None, cluster_name=None, element=None ):
        if ( calc is None ):
            raise ValueError( "No calculator given. Has to be a CE calculator (with C++ support)" )
        if ( cluster_name is None ):
            raise ValueError( "No cluster name given!" )
        if ( element is None ):
            raise ValueError( "No element given!" )
        self.fast_cluster_tracker = ce_updater.ClusterTracker( calc.updater, cluster_name, element )
        super(NetworkObserver,self).__init__()
        self.name = "NetworkObserver"
        self.res = {
            "avg_size":0.0,
            "avg_size_sq":0.0,
            "number_of_clusters":0
        }
        self.max_size = 0
        self.indx_max_cluster = []
        self.atoms_max_cluster = None

    def __call__( self, system_changes ):
        new_res = self.fast_cluster_tracker.get_cluster_statistics_python()
        for key in self.res.keys():
            self.res[key] += new_res[key]
        if ( new_res["max_size"] > self.max_size ):
            self.max_size = new_res["max_size"]
            self.atoms_max_cluster = self.calc.atoms.copy()
            clust_indx = self.fast_cluster_tracker.atomic_clusters2group_indx_python()
            self.indx_max_cluster = clust_indx

    def get_atoms_with_largest_cluster( self, highlight_element="Na" ):
        """
        Returns the atoms object which had the largest cluster and change the element
        of the atoms in the cluster to *highlight_element*
        """
        explored_grp_indices = []
        largest_cluster = []
        # Locate the largest cluster
        for i in range(0,len(self.atoms_max_cluster)):
            if ( self.atoms_max_cluster[i] in explored_grp_indices ):
                continue
            current_cluster = []
            for j in range(0,len(self.atoms_max_cluster)):
                if ( self.atoms_max_cluster[j] == self.atoms_max_cluster[i] ):
                    current_cluster.append(j)
            if ( len(current_cluster) > len(largest_cluster) ):
                largest_cluster = current_cluster

        for indx in largest_cluster:
            self.atoms_max_cluster[indx].symbol = highlight_element
        return self.atoms_max_cluster
