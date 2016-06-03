#include<iostream>
#include<string>
#include<map>
#include<vector>

#include "TH2F.h"
#include "TObject.h"
#include "TCanvas.h"
class CrossCount{
	public:
		CrossCount();

    int AddMuonData(const std::string& muon_data_name,const std::string& tree_name);
    int AddSpinDatabaseData(const std::string& spin_data_base_name);

    int Run();

    int SaveFigures(const std::string& out_dir);
		~CrossCount();
	private:
    std::string muon_data_name_;
    std::map<std::string,bool> state_;

    struct MuonData{
      float q;
      float pz;
      int clockcross;
      float wness;
    };

    struct SpinData{
      float pol_blue;
      float pol_yell;
      std::vector<int> bunch_pol;
    };
    // cross_count_[arm][charge]
    // arm: 0 = south, 1 = north
    // charge: 0 = negative, 1 = positive
    TH2F* cross_count_[2][2];
    TH2F* cross_count_sum_;

    // spin_patterns_[spin_pat]
    TH1F* spin_patterns_[5];

    // muon_track_spin_patterns_[arm][charge][spin_pat]
    TH1F* muon_track_spin_patterns_[2][2][5];

    // muon_track_spin_patterns_summed_[spin_pat]
    TH1F* muon_track_spin_patterns_summed_[5];

    TCanvas* spin_pattern_count;
    std::map<int,std::vector<MuonData> > muons_;
    std::map<int,SpinData> spin_data_;
    std::vector<TObject*> registry_;

    // spin_pat_map_
    //
    // Maps the spin db interger to a spin pattern
    //
    // ex: +- : blue polarization +, yellow polarization -
    // 0 -> "++"
    // 1 -> "-+"
    // 2 -> "+-"
    // 3 -> "--"
    // 4 -> "EMPTY"
    std::map<int,std::string> spin_pat_map_;
};
