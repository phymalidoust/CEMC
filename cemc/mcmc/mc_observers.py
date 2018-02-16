import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import copy
import numpy as np

class MCObserver( object ):
    def __init__( self ):
        pass

    def __call__( self, system_changes ):
        """
        Gets information about the system changes and can perform some action
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

    def __call__( self, system_changes ):
        if ( self.lowest_energy_atoms is None or self.lowest_energy_cf is None ):
            self.lowest_energy_cf = self.ce_calc.get_cf()
            self.lowest_energy_atoms = self.ce_calc.atoms.copy()
            self.lowest_energy = self.mc_obj.current_energy
            return

        if ( self.mc_obj.current_energy < self.lowest_energy ):
            self.lowest_energy = self.mc_obj.current_energy
            self.lowest_energy_atoms = self.ce_calc.atoms.copy()
            self.lowest_energy_cf = self.ce_calc.get_cf()

class SGCObserver(MCObserver):
    def __init__( self, ce_calc, mc_obj, n_singlets ):
        self.ce_calc = ce_calc
        self.mc = mc_obj

        self.quantities = {
            "singlets":np.zeros( n_singlets, dtype=np.float64 ),
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
        self.quantities["energy"] += self.mc.current_energy
        self.quantities["energy_sq"] += self.mc.current_energy**2
        self.quantities["singl_eng"] += new_singlets*self.mc.current_energy
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
