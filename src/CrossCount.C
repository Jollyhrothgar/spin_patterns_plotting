#include "CrossCount.h"

#include<cstdio>
#include<fstream>
#include<map>
#include<sstream>
#include<memory>
#include<vector>

#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

CrossCount::CrossCount(){
	std::cout << "CrossCount instantiated at" << this << std::endl;
  state_["MUON_DATA_READY"] = false;
  state_["SPIN_DATABASE_READY"] = false;
  spin_pat_map_[0] = "++";
  spin_pat_map_[1] = "-+";
  spin_pat_map_[2] = "+-";
  spin_pat_map_[3] = "--";
  spin_pat_map_[4] = "**";
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
  Float_t Wness;
  t->SetBranchAddress("charge",&charge);
  t->SetBranchAddress("clockcross",&clockcross);
  t->SetBranchAddress("pz",&pz);
  t->SetBranchAddress("Run_Number",&run_number);
  t->SetBranchAddress("Wness",&Wness);

  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    MuonData d;
    d.q = charge;
    d.pz = pz;
    d.clockcross = clockcross;
    d.wness = Wness;
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
    SpinData data;
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
        // if no spin patterns, pattern defaults to -999
        if(spin_pat < 0){
          continue;
        }
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
  sprintf(name,"%s/figures.root",out_dir.c_str());
  TFile * out_file = new TFile(name,"RECREATE");
  out_file->cd();
  if(state_["MUON_DATA_READY"]) {
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        cross_count_[i][j]->Write();
      }
    }
    cross_count_sum_->Write();
  }
  if(state_["SPIN_DATABASE_READY"]){
    for(int i = 0; i < 5; i++){
      spin_patterns_[i]->Write();
    }
  }
  if(state_["MUON_DATA_READY"] && state_["SPIN_DATABASE_READY"]) {
    for(int i = 0; i < 5; i++){
      muon_track_spin_patterns_summed_[i]->Write();
      for(int arm_i = 0; arm_i < 2; arm_i++){
        for(int charge_i = 0; charge_i < 2; charge_i++){
          muon_track_spin_patterns_[arm_i][charge_i][i]->Write();
        }
      }
    }
  }
  out_file->Write();
  out_file->Close();
  delete out_file;
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
        cross_count_[i][j] = new TH2F(
            name,title,muons_.size(),-0.5,(float)muons_.size()+0.5,120,-0.5,119.5);  
        registry_.push_back(cross_count_[i][j]);
      }
    }
    cross_count_sum_ = new TH2F(
        "crossing_dist_sum",
        "Crossing Distribution;Run Index;Bunch Index;Crossing Count",
        muons_.size(),-0.5,(float)muons_.size()+0.5,120,-0.5,119.5
        );
    registry_.push_back(cross_count_sum_);

    int run_counter = 0;
    for(auto run = muons_.begin(); run != muons_.end(); ++run){
      std::vector<MuonData> runs = run->second;
      for(auto datum = runs.begin(); datum != runs.end(); ++datum){
        MuonData d = *datum;
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
        cross_count_[arm_i][charge_i]->Fill(run_counter,d.clockcross);
        cross_count_sum_->Fill(run_counter,d.clockcross);
      }
      run_counter++;
    }
  }
  if(state_["SPIN_DATABASE_READY"]){
    for(int i = 0; i < 5; i++){
      sprintf(name,"spin_patterns_%d",i);
      sprintf(title,
          "Spin Pattern %s Counts;Spin Pattern %s;Counts",
          spin_pat_map_[i].c_str(), spin_pat_map_[i].c_str());
      spin_patterns_[i] = new TH1F(name,title,5,-0.5,4.5);
      registry_.push_back(spin_patterns_[i]);
    }
    for (auto run = spin_data_.begin(); run != spin_data_.end(); ++run){
      SpinData spin_data = run->second;
      for(auto pat = spin_data.bunch_pol.begin(); pat != spin_data.bunch_pol.end();++pat)
      {
        //std::cout << *pat << std::endl; 
        spin_patterns_[*pat]->Fill(*pat);
      }
    }
  }
  if(state_["SPIN_DATABASE_READY"] && state_["MUON_DATA_READY"]){
    for(int i = 0; i < 5; i++){
      sprintf(name,"muon_track_spin_pattern_summed_p%d",i);
      sprintf(
          title,
          "Muon Track Spin Pattern Summed, %s;Spin Pattern %s,Count",
          spin_pat_map_[i].c_str(),
          spin_pat_map_[i].c_str()
          );
      muon_track_spin_patterns_summed_[i] = new TH1F(name,title,5,-0.5,4.5);
      registry_.push_back(muon_track_spin_patterns_summed_[i]);
      for(int arm_i = 0; arm_i < 2; arm_i++){
        for(int charge_i = 0; charge_i < 2; charge_i++){
          sprintf(name,"muon_track_spin_pattern_a%d_q%d_p%d",arm_i,charge_i,i);
          sprintf(
              title,
              "Muon Track Spin Pattern Arm %d Charge %d %s;Spin Pattern %s,Count",
              arm_i,
              charge_i,
              spin_pat_map_[i].c_str(),
              spin_pat_map_[i].c_str()
              );
          muon_track_spin_patterns_[arm_i][charge_i][i] = new TH1F(name,title,5,-0.5,4.5);
          registry_.push_back(muon_track_spin_patterns_[arm_i][charge_i][i]);
        }
      }
    }
    for(auto run = muons_.begin(); run != muons_.end(); ++run){
      std::vector<MuonData> runs = run->second;
      int run_number = run->first;
      for(auto datum = runs.begin(); datum != runs.end(); ++datum){
        MuonData d = *datum;
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
        if(spin_data_.find(run_number) == spin_data_.end()) continue;
        int pattern = spin_data_[run_number].bunch_pol[d.clockcross];
        muon_track_spin_patterns_[arm_i][charge_i][pattern]->Fill(pattern);
        muon_track_spin_patterns_summed_[pattern]->Fill(pattern);
      }

    }
  }
  return 0;
}

CrossCount::~CrossCount(){
  std::cout << "CrossCount Instance Finished, destroying " << this << std::endl;
  for(auto i = registry_.begin(); i != registry_.end(); ++i){
    if(*i){
      delete *i;
    }
  }
}
