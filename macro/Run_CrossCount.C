void Run_CrossCount(
  std::string muon_data_name = "../../Run13WnessTreeWithFvtxAndRpc_Data.root",
  std::string spin_data_base_name = "../../spin_data_local.txt",
  std::string out_dir = "/home/jollyhrothgar/root_trees/crossing_dist/macro/figures",
  std::string lib = "libCrossCount.so"
    ){
  gSystem->Load(lib.c_str());
  CrossCount cc;
  cc.AddMuonData(muon_data_name,"newsngmuons_basic_cut");
  cc.AddSpinDatabaseData(spin_data_base_name);
  //cc.AddSpinDatabaseData("alksjda");
  cc.Run();
  cc.SaveFigures(out_dir);
}
