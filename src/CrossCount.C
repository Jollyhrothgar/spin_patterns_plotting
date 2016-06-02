#include<cstdio>
#include<fstream>
#include<map>
#include<sstream>
#include "CrossCount.h"
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"

CrossCount::CrossCount(){
	std::cout << "CrossCount instantiated at" << this << std::endl;
  state_["MUON_DATA_READY"] = false;
  state_["SPIN_DATABASE_READY"] = false;
}

int CrossCount::AddMuonData(const std::string& muon_data_name, const std::string& tree_name = "newsngmuons_basic_cut"){
  muon_data_name_ = muon_data_name;
  TFile* f = new TFile(muon_data_name_.c_str(),"READ");
  TTree* t = (TTree*)f->Get(tree_name.c_str());
  if (!f) {
    std::cerr << "Problem loading ROOT file: " << muon_data_name_ << std::endl; 
    return 1;
  }
  if (!t){ 
    std::cerr << "Problem loading TTree from ROOT File."
      << std::endl << "TTree: " << tree_name 
      << std::endl << "TFile: " << muon_data_name_ 
      << std::endl;
    return 2;
  }

  Int_t run_number;
  Int_t clockcross;
  Float_t charge;
  Float_t pz;
  t->SetBranchAddress("charge",&charge);
  t->SetBranchAddress("clockcross",&clockcross);
  t->SetBranchAddress("pz",&pz);
  t->SetBranchAddress("Run_Number",&run_number);

  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    data d;
    d.q = charge;
    d.pz = pz;
    d.clockcross = clockcross;
    //std::cout << d.q << ", " << d.pz << ", " << d.clockcross << ", " << run_number << std::endl;
    muons_[run_number].push_back(d);
  }

  std::cout << muons_.size() << " runs loaded." << std::endl;
  delete t;
  state_["MUON_DATA_READY"] = true;
  return 0;
}


// Expects CSV
int CrossCount::AddSpinDatabaseData(const std::string& spin_data_base_name) {
  std::ifstream infile(spin_data_base_name.c_str());
  if(! infile){
    std::cerr << "Problem opening" << spin_data_base_name << std::endl;
    return 1;
  }
  std::string line = "";
  bool first_line = true;
  while(getline(infile,line)){
    // skip header
    if(first_line) {
      first_line = false;
      continue;
    }
    spin_db data;
    int run_number;
    std::stringstream ss(line);
    std::string token;
    int token_counter = 0;
    while(getline(ss,token,',')){
      std::stringstream convert(token);
      if( token_counter == 0 ){
        convert >> run_number;
      } else if(token_counter == 1){
        convert >> data.pol_blue;
      } else if(token_counter == 2){
        convert >> data.pol_yell;
      } else {
        int spin_pat;
        convert >> spin_pat;
        data.bunch_pol.push_back(spin_pat);
      }
      token_counter++;
    }

    if(data.bunch_pol.size() != 120) {
      std::cout << "not enough spin data" << std::endl;
      continue;
    }
    // exactly one entry per run
    spin_data_[run_number] = data;
  }
  infile.close();

  // Check if extraction worked...
  if(spin_data_.size() == 0) {
    std::cerr << "no data was loaded from " << spin_data_base_name << std::endl;
    return 2;
  }
  state_["SPIN_DATABASE_READY"] = true;
  std::cout << "Successfully loaded " << spin_data_base_name << std::endl;
  return 0;
}

int CrossCount::SaveFigures(const std::string& out_dir){
  char name[246];
  if(state_["MUON_DATA_READY"]) {
    sprintf(name,"%s/figures.root",out_dir.c_str());
    TFile * out_file = new TFile(name,"RECREATE");
    out_file->cd();
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        cross_count[i][j]->Write();
      }
    }
    cross_count_sum->Write();
    out_file->Close();
    delete out_file;
  }
  return 0;
}

int CrossCount::Run(){
  std::cout << "Generating and saving figures" << std::endl;
  char name[256];
  char title[256];
  if( state_["MUON_DATA_READY"] ) {
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        sprintf(name,"crossing_dist_arm%d_charge%d",i,j);
        sprintf(title,
            "Crossing Distribution A%d Q%d;Run Index;Bunch Index;Crossing Count",i,j);
        cross_count[i][j] = new TH2F(
            name,title,muons_.size(),-0.5,(float)muons_.size()+0.5,120,-0.5,119.5);  
      }
    }

    cross_count_sum = new TH2F(
        "crossing_dist_sum",
        "Crossing Distribution;Run Index;Bunch Index;Crossing Count",
        muons_.size(),-0.5,(float)muons_.size()+0.5,120,-0.5,119.5
        );

    int run_counter = 0;
    for(auto run = muons_.begin(); run != muons_.end(); ++run){
      std::vector<data> runs = run->second;
      for(auto datum = runs.begin(); datum != runs.end(); ++datum){
        data d = *datum;
        int arm_i;
        int charge_i;
        if(d.q > 0){
          charge_i = 1;
        } else {
          charge_i = 0;
        }

        if(d.pz > 0){
          arm_i = 1;
        } else {
          arm_i = 0;
        }
        d.clockcross = (d.clockcross+5)%120;
        cross_count[arm_i][charge_i]->Fill(run_counter,d.clockcross);
        cross_count_sum->Fill(run_counter,d.clockcross);
      }
      run_counter++;
    }
  }
  return 0;
}

CrossCount::~CrossCount(){
  std::cout << "CrossCount Instance Finished, destroying" << this << std::endl;
}
