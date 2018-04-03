import unittest
from ase.build import bulk
try:
    from ase.ce import BulkCrystal
    from cemc.wanglandau import WangLandauInit, WangLandau, WangLandauDBManager
    has_CE = True
except Exception as exc:
    print (exc)
    has_CE = False

db_name = "temp_db.db"
wl_db_name = "wanglandau_test_init.db"
bc_kwargs = {
    "crystalstructure":"fcc",
    "size":[3,3,3],
    "basis_elements":[["Al","Mg"]],
    "db_name":db_name,
    "conc_args":{"conc_ratio_min_1":[[1,0]],"conc_ratio_max_1":[[0,1]]},
    "max_cluster_dia":4,
    "a":4.05
}

eci = {
    "c1_0":1.0,
    "c2_1000_1_00":1.0
}
class TestInitWLSim( unittest.TestCase ):
    def test_no_throw( self ):
        no_throw = True
        if ( not has_CE ):
            self.skipTest( "ASE version does not have CE" )

        msg = ""
        try:
            initializer = WangLandauInit( wl_db_name )
            T = [1000,10]
            comp = {"Al":0.5,"Mg":0.5}
            initializer.insert_atoms( bc_kwargs, size=[5,5,5], T=T, n_steps_per_temp=10, eci=eci, composition=comp )
            initializer.prepare_wang_landau_run( [("id","=","1")] )
            atoms = initializer.get_atoms( 1, eci )
            db_manager = WangLandauDBManager( wl_db_name )
            runID = db_manager.get_next_non_converged_uid( 1 )
            if ( runID == -1 ):
                raise ValueError( "No new Wang Landau simulation in the database!" )
            simulator = WangLandau( atoms, wl_db_name, runID, fmin=1.8 )
            simulator.run_fast_sampler( mode="adaptive_windows" )
        except Exception as exc:
            msg = str(exc)
            no_throw = False
        self.assertTrue( no_throw, msg=msg )

if __name__ == "__main__":
    unittest.main()