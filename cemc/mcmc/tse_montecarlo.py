from cemc.mcmc import SGCMonteCarlo

class TSEMonteCarlo(SGCMonteCarlo):
    def __init__( self, atoms, temp, mpicomm=None, symbols=None, mpicomm=None, logfile="", size_window_width=10, \
    chem_pot=None, max_size_reactant=2, min_size_product=20 ):
        super(TSEMonteCarlo,self).__init__(atoms,temp,mpicomm=mpicomm,symbols=symbols)

        self.size_window_width = size_window_width
        self.max_size_reactant = max_size_reactant
        self.min_size_product = min_size_product
        self.chemical_potential = chem_pot
        
