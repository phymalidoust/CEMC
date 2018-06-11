from cemc.mcmc.nucleation_sampler import NucleationSampler, Mode
from cemc.mcmc import SGCMonteCarlo
from cemc.mcmc.mc_observers import NetworkObserver

class SGCNucleation( SGCMonteCarlo ):
    """
    Class to perform Monte Carlo simulations where the main objective is to
    compute the free energy for different network sizes

    See py:class:`cemc.mcmc.SGCMonteCarlo` for parameter explination

    :Keyword arguments:
        * *nucleation_sampler* Instance of :py:class:`cemc.mcmc.NucleationSampler`
        * *network_name* Name of the network to sutdy see `cemc.mcmc.NetworkObserver`
        * *network_element* Element in the network see `cemc.mcmc.NetworkObserver`
        * *chem_pot* The chemical potential at which to perform the simulations see :py:meth:`cemc.mcmc.SGCMonteCarlo.runMC`
        * *allow_solutes* Allow solute atoms outside the network (i.e. dispersed in the matrix). Default is *True*
    """
    def __init__( self, atoms, temp, **kwargs ):
        self.nuc_sampler = kwargs.pop("nucleation_sampler")
        kwargs["mpicomm"] = None
        self.network_name = kwargs.pop("network_name")
        self.network_element = kwargs.pop("network_element")
        chem_pot = kwargs.pop("chem_pot")
        self.allow_solutes = True
        if ( "allow_solutes" in kwargs.keys() ):
            self.allow_solutes = kwargs.pop("allow_solutes")
        super( SGCNucleation, self ).__init__( atoms, temp, **kwargs)
        self.chemical_potential = chem_pot

        self.network = NetworkObserver( calc=self.atoms._calc, cluster_name=self.network_name, element=self.network_element )
        self.attach( self.network )


        if ( self.allow_solutes ):
            self.log( "Solute atoms in cluster and outside is allowed" )
        else:
            self.log( "Solute atoms are only allowed in the cluster" )

    def accept( self, system_changes ):
        """
        Accept the trial move

        :param system_changes: See :py:meth:`cemc.mcmc.Montecarlo.accept`
        """
        move_accepted = SGCMonteCarlo.accept( self, system_changes )
        in_window,stat = self.nuc_sampler.is_in_window(self.network,retstat=True)
        if ( not self.allow_solutes ):
            new_size = stat["max_size"]
            cur_size = self.nuc_sampler.current_cluster_size
            if ( new_size != cur_size+1 and new_size != cur_size-1 ):
                in_window = False
            else:
                self.nuc_sampler.current_cluster_size = new_size
        return move_accepted and in_window

    def get_trial_move(self):
        """
        Perform a trial move
        """
        if ( not self.nuc_sampler.is_in_window(self.network) ):
            raise RuntimeError( "System is outside the window before the trial move is performed!" )
        return SGCMonteCarlo.get_trial_move(self)

    def run( self, nsteps=10000 ):
        """
        Run samples in each window until a desired precission is found

        :param nsteps: Number of MC steps in each window
        """
        if ( self.nuc_sampler.nucleation_mpicomm is not None ):
            self.nuc_sampler.nucleation_mpicomm.barrier()
        for i in range(self.nuc_sampler.n_windows):
            self.log( "Window {} of {}".format(i,self.nuc_sampler.n_windows) )
            self.nuc_sampler.current_window = i
            self.reset()
            self.nuc_sampler.bring_system_into_window(self.network)

            self.nuc_sampler.mode = Mode.equillibriate
            self.estimate_correlation_time()
            self.equillibriate()
            self.nuc_sampler.mode = Mode.sample_in_window

            current_step = 0
            while( current_step < nsteps ):
                current_step += 1
                self._mc_step()
                self.nuc_sampler.update_histogram(self)
                self.network.reset()

        if ( self.nuc_sampler.nucleation_mpicomm is not None ):
            self.nuc_sampler.nucleation_mpicomm.barrier()

    def remove_snapshot_observers(self):
        """
        Remove all Snapshot observers from the observers
        """
        self.observers = [obs for obs in self.observers if obs.name != "Snapshot"]

    def remove_network_observers(self):
        """
        Remove NetworkObservers
        """
        self.observers = [obs for obs in self.observers if obs[1].name != "NetworkObserver"]

    def is_reactant(self):
        """
        Returns true if the current state is in the reactant region
        """
        if ( self.max_size_reactant is None ):
            raise ValueError( "Maximum cluster size to be characterized as reactant is not set!" )

        stat = self.network.get_statistics()
        return stat["max_size"] < self.max_size_reactant

    def is_product(self):
        """
        Return True if the current state is a product state
        """
        if ( self.min_size_product is None ):
            raise ValueError( "Minimum cluster size to be characterized as product is not set!" )
        stat = self.network.get_statistics()
        return stat["max_size"] >= self.min_size_product

    def merge_product_and_reactant_path( self, reactant_traj, product_traj, reactant_symb, product_symb ):
        """
        Merge the product and reactant path into one file
        """
        folder = reactant_traj.rpartition("/")[0]
        symb_merged = folder+"/reaction2product.txt"
        reactant_symbols = []
        product_symbols = []

    def save_list_of_lists(self,fname,data):
        """
        Save a list of lists into a text file
        """
        with open(fname,'w') as outfile:
            for sublist in data:
                for entry in sublist:
                    outfile.write("{} ".format(entry))
                outfile.write("\n")

    def reset(self):
        """
        Overrides the parents method
        """
        super(SGCNucleation,self).reset()
        self.current_energy = 1E10
        self.network.reset()

    def read_list_of_lists(self,fname,dtype="str"):
        """
        Read list of lists
        """
        supported_dtypes = ["str"]
        if ( dtype not in supported_dtypes ):
            raise ValueError( "dtype hsa to be one of {}".format(supported_dtypes))

    def symbols2uint( self, symbols, description ):
        """
        Convert an array of symbols into a numpy array of indices to the desctiption array
        """
        nparray = np.zeros( len(symbols), dtype=np.uint8 )
        for i,symb in enumerate(symbols):
            nparray[i] = desctiption.index(symb)
        return nparray

    def uint2symbols( self, nparray, description ):
        """
        Convert uint8 array to symbols array
        """
        symbs = []
        for i in range(len(nparray)):
            symbs.append( description[nparray[i]] )
        return symbs

    def merge_reference_path( self, res_reactant, res_product ):
        """
        This store the reference path into a JSON
        """
        res_reactant["energy"] = res_reactant["energy"][::-1]
        res_reactant["symbols"] = res_reactant["symbols"][::-1]
        combined_path = {}
        combined_path["energy"] = res_reactant["energy"]+res_product["energy"]
        combined_path["symbols"] = res_reactant["symbols"]+res_product["symbols"]
        return combined_path

    def save_path( self, fname, res ):
        """
        Stores the path result to a JSON file
        """
        res["min_size_product"] = self.min_size_product
        res["max_size_reactant"] = self.max_size_reactant
        with open(fname,'w') as outfile:
            json.dump(res,outfile)

    def sweep(self):
        """
        Performs one MC sweep
        """
        for i in range(len(self.atoms)):
            self._mc_step()

    def set_mode( self, mode ):
        """
        Set the mode
        """
        known_modes = ["bring_system_into_window","sample_in_window","equillibriate","transition_path_sampling"]
        if ( mode not in known_modes ):
            raise ValueError( "Mode has to be one of {}".format(known_modes))

        if ( mode == "bring_system_into_window" ):
            self.nuc_sampler.mode = Mode.bring_system_into_window
        elif ( mode == "sample_in_window" ):
            self.nuc_sampler.mode = Mode.sample_in_window
        elif ( mode == "equillibriate" ):
            self.nuc_sampler.mode = Mode.sample_in_window
        elif ( mode == "transition_path_sampling" ):
            self.nuc_sampler.mode = Mode.transition_path_sampling

    def set_state( self, symbols ):
        """
        Sets the state of the system
        """
        self.atoms._calc.set_symbols(symbols)

    def show_statistics(self,path):
        """
        Show a plot indicating if the path is long enough
        """
        product_indicator = []
        reactant_indicator = []
        for state in path:
            self.network.reset()
            self.set_state(state)
            self.network(None)
            if ( self.is_product() ):
                product_indicator.append(1)
            else:
                product_indicator.append(0)

            if ( self.is_reactant() ):
                reactant_indicator.append(1)
            else:
                reactant_indicator.append(0)
        hB = np.cumsum(product_indicator)
        hA = np.cumsum(reactant_indicator)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( hB, label="Product indicator" )
        ax.plot( hA, label="Reactant indicator" )
        ax.set_xlabel( "MC sweeps" )
        ax.legend()
        return fig

    def find_transition_path( self, initial_cluster_size=None, max_size_reactant=None, min_size_product=None, path_length=1000, max_attempts=100, folder="." ):
        """
        Find one transition path
        """
        if ( initial_cluster_size is None ):
            raise ValueError( "Initial cluster size not given!" )
        if ( max_size_reactant is None ):
            raise ValueError( "The maximum cluster size allowed for the state to be characterized as reactant is not given!" )
        if ( min_size_product is None ):
            raise ValueError( "The minimum size of cluster allowed for the state to be characterized as product is not given!" )

        self.mode = Mode.transition_path_sampling
        self.max_size_reactant = max_size_reactant
        self.min_size_product = min_size_product

        found_reactant_origin = False
        found_product_origin = False
        self.remove_network_observers()
        self.attach( self.network, interval=len(self.atoms) )

        num_reactants = 0
        num_products = 0
        default_trajfile = folder+"/default_trajfile.traj"
        reactant_file = folder+"/trajectory_reactant.traj"
        product_file = folder+"/trajectory_product.traj"
        reference_path_file = folder+"/reference_path.json"

        self.network.reset()
        self.network.grow_cluster( initial_cluster_size )

        init_symbols = [atom.symbol for atom in self.atoms]
        target = "both"
        reactant_res = {}
        product_res = {}
        for attempt in range(max_attempts):
            self.reset()
            self.atoms._calc.set_symbols(init_symbols)
            try:
                res = self.find_one_transition_path( path_length=path_length/2, trajfile=default_trajfile, target=target )
            except DidNotReachProductOrReactantError as exc:
                self.log( str(exc) )
                self.log ( "Trying one more time" )
                continue

            if ( res["type"] == "reactant" ):
                num_reactants += 1
                target = "product" # Reactant is found, search only for products
                if ( not found_reactant_origin ):
                    os.rename( default_trajfile,reactant_file)
                    found_reactant_origin = True
                    reactant_res = copy.deepcopy(res)

            elif ( res["type"] == "product" ):
                num_products += 1
                target = "reactant" # Product is found search only for reactant
                if ( not found_product_origin ):
                    os.rename( default_trajfile,product_file)
                    found_product_origin = True
                    product_res = copy.deepcopy(res)

            if ( os.path.exists(default_trajfile) ):
                os.remove(default_trajfile)

            if ( found_product_origin and found_reactant_origin ):
                combined_path = self.merge_reference_path(reactant_res,product_res)
                self.save_path( reference_path_file, combined_path )
                self.log( "Found a path to the product region and a path to the reactant region" )
                self.log( "They are stored in {} and {}".format(product_file,reactant_file))
                self.log( "The reference path is stored in {}".format(reference_path_file) )
                self.show_statistics(combined_path["symbols"])
                return
            self.log( "Attempt: {} of {} ended in {} region".format(attempt,max_attempts,res["type"]) )
        msg = "Did not manage to find both a configuration in the product region and the reactant region\n"
        raise RuntimeError( msg )


    def find_one_transition_path( self, path_length=1000, trajfile="default.traj", target="both" ):
        """
        Finds a transition path by running random samples
        """
        supported_targets = ["reactant","product","both"]
        if ( target not in supported_targets ):
            raise ValueError( "Target has to be one of {}".format(supported_targets) )

        # Check if a snapshot tracker is attached
        traj = TrajectoryWriter( trajfile, mode="w" )
        current_step = 0
        result = {}
        symbs = []
        unique_symbols = []
        for atom in self.atoms:
            if ( atom.symbol not in unique_symbols ):
                unique_symbols.append(atom.symbol)

        output_every_sec = 30
        now = time.time()
        energies = []
        result = {}
        for sweep in range(path_length):
            self.network.reset()
            if ( time.time() - now > output_every_sec ):
                self.log( "Sweep {} of {}".format(sweep,path_length))
                now = time.time()
            self.sweep()
            self.network(None) # Explicitly enforce a construction of the network
            energies.append(self.current_energy)
            symbs.append( [atom.symbol for atom in self.atoms] )
            atoms = self.network.get_atoms_with_largest_cluster(prohibited_symbols=unique_symbols)
            if ( atoms is None ):
                traj.write(self.atoms)
            else:
                traj.write(atoms)

            if ( target == "reactant" ):
                if ( self.is_product() ):
                    # Terminate before the desired path length is reached
                    result["type"] = "product"
                    result["symbols"] = symbs
                    result["energy"] = energies
                    return result
            elif ( target == "product" ):
                if ( self.is_reactant() ):
                    result["type"] = "reactant"
                    result["symbols"] = symbs
                    result["energy"] = energies
                    # Terminate before the desired path length is reached
                    return result

        traj.close()
        if ( self.is_reactant() ):
            result["type"] = "reactant"
        elif( self.is_product() ):
            result["type"] = "product"
        else:
            stat = self.network.get_statistics()
            max_size = stat["max_size"]
            msg = "State did not end up in product or reactant region. Increase the number of sweeps.\n"
            msg += "Max. cluster size {}. Max cluster size reactants {}. Min cluster size products {}".format(max_size,self.max_size_reactant,self.min_size_product)
            raise DidNotReachProductOrReactantError( msg )

        result["symbols"] = symbs
        result["energy"] = energies
        return result
