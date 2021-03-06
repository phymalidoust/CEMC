from cemc.mcmc import Montecarlo, TooFewElementsError
from itertools import combinations_with_replacement, combinations
import numpy as np
import time

class ChemicalPotentialROI(object):
    """
    Class that identifies interesting chemical potentials to study.
    The algoritthm performs the following steps
    """
    def __init__(self, atoms, conc_step=0.1, temperature=100, symbols=[]):
        self.atoms = atoms
        self.conc_step = conc_step
        self.temperature = temperature
        self.symbols =symbols
        self.status_every_sec = 30

    def _log(self, msg):
        """
        Logging
        """
        print(msg)

    def _estimate_internal_energy(self, conc, sweeps=2):
        """
        Estimates the internal energy of one structure
        """
        try:
            self.atoms._calc.set_composition(conc)
            mc = Montecarlo(self.atoms, self.temperature)
            mc.runMC(mode="fixed", steps=sweeps*len(self.atoms), equil=False)
            energy = mc.get_thermodynamic()["energy"]
        except TooFewElementsError as exc:
            energy = 1.0
        return energy

    def _find_energies(self, sweeps=2):
        """
        Estimate the energy at all compositions
        """
        singlets = self.atoms._calc.get_singlets()
        if len(singlets) != len(self.symbols)-1:
            msg = "The number of symbols does not match the number of basis singlets\n"
            msg += "Number of symbols: {}\n".format(len(self.symbols))
            msg += "Number of singlet terms: {}\n".format(len(singlets))
            msg += "It should be one more symbol compared to the number of singlets"
            raise ValueError(msg)

        n_concs = len(singlets)
        template_conc = np.linspace(0.0, 1.0, int(1.0/self.conc_step))
        result = []
        now = time.time()
        counter = 0
        n_comb = len(template_conc)**len(singlets)
        for comp in combinations_with_replacement(template_conc, n_concs):
            counter += 1
            if np.sum(comp) > 1.0:
                continue

            if time.time()-now > self.status_every_sec:
                self._log("Running composition {} of {}".format(counter,n_comb))
                now = time.time()
            conc_dict = {symb:value for symb,value in zip(self.symbols[1:],comp)}
            conc_dict[self.symbols[0]] = 1.0-np.sum(comp)
            energy = self._estimate_internal_energy(conc_dict, sweeps=sweeps)/len(self.atoms)
            res = {}
            res["conc"] = conc_dict
            res["energy"] = energy
            cf = self.atoms._calc.get_cf()
            singl = {key:value for key,value in cf.items() if key.startswith("c1")}
            res["singlets"] = singl
            result.append(res)
        return result

    def chemical_potential_roi(self, sweeps=2):
        """
        Finds the chemical potential that makes all structures have the same
        energy as the structure with lowest enthalpy of formation
        """
        ref_energies = {}
        singlets = {}
        for ref_symb in self.symbols:
            comp = {key:0.0 for key in self.symbols}
            comp[ref_symb] = 1.0
            self.atoms._calc.set_composition(comp)
            ref_energies[ref_symb] = self.atoms._calc.calculate(self.atoms, ["energy"], [])/len(self.atoms)
            singl = self.atoms._calc.get_cf()
            singl = {key:value for key,value in singl.items() if key.startswith("c1")}
            singlets[ref_symb] = singl

        energies = self._find_energies(sweeps=sweeps)
        e_form = np.zeros(len(energies))
        for i,entry in enumerate(energies):
            e_form[i] = entry["energy"]
            for symb in self.symbols:
                e_form[i] -= entry["conc"][symb]*ref_energies[symb]

        min_e_form = np.argmin(e_form)
        lowest_energy_form = energies[min_e_form]

        chemical_potentials = []
        N = len(self.symbols)-1
        A = np.zeros((N,N))
        rhs = np.zeros(N)
        key_indx = {key:i for i,key in enumerate(singlets[self.symbols[0]].keys())}
        mu_roi = []
        for comb in combinations(self.symbols,len(self.symbols)-1):
            row = 0
            for symb in comb:
                for key2,indx2 in key_indx.items():
                    A[row,indx2] = singlets[symb][key2]-lowest_energy_form["singlets"][key2]
                rhs[row] = ref_energies[symb]-lowest_energy_form["energy"]
                row += 1

            used_pseudo_inverse = False
            try:
                mu = np.linalg.solve(A,rhs)
            except np.linalg.LinAlgError:
                inv = np.linalg.pinv(A)
                mu = inv.dot(rhs)
                used_pseudo_inverse = True
            mu_dict = {}
            for key,indx in key_indx.items():
                mu_dict[key] = mu[indx]
            mu_dict["symbs"] = comb
            mu_dict["pseudo_inv"] = used_pseudo_inverse
            mu_roi.append(mu_dict)
        return mu_roi

    @staticmethod
    def list_match(list1, list2):
        """
        Checks if the elements in two lists match

        :param list1: First entry
        :param list2: Second entry
        """
        if len(list1) != len(list2):
            return False

        for entry in list1:
            if entry not in list2:
                return False
        return True

    @staticmethod
    def suggest_mu(mu_roi=None, N=10, extend_fraction=0.1, elements=None):
        """
        This function suggests mu that can be used for exploring the
        parameter space based on the region of interest found by
        the function chemical_potential_roi
        """
        suggested_sampling_lines = []
        names = [key for key in mu_roi[0].keys() if key.startswith("c1")]

        if elements is None:
            phase_combinations = combinations(mu_roi,2)
        else:
            if len(elements) != 2:
                raise ValueError("You have to specify start and end points!")
            phase_combinations = [[]]
            for entry in mu_roi:
                if ChemicalPotentialROI.list_match(entry["symbs"], elements[0]) or \
                    ChemicalPotentialROI.list_match(entry["symbs"], elements[1]):
                    phase_combinations[0].append(entry)

        for phase_comb in phase_combinations:
            mu_array1 = np.array([phase_comb[0][key] for key in names])
            mu_array2 = np.array([phase_comb[1][key] for key in names])
            suggested_sampling_lines.append(ChemicalPotentialROI.linear_sampling(mu_array1, mu_array2, N, extend_fraction=extend_fraction))
        return suggested_sampling_lines, names

    @staticmethod
    def linear_sampling(start, stop, N, extend_fraction=0.1):
        """
        Construct one array with linear sampling from start to stop
        """
        diff = stop-start
        unit_vector = diff/np.sqrt(np.sum(diff**2))
        mu = np.zeros((N,len(diff)))
        start = start-extend_fraction*diff
        stop = stop+extend_fraction*diff

        # Update the difference
        diff = stop-start
        step = diff/(N-1)
        for i in range(N):
            mu[i,:] = start + step*i
        return mu
