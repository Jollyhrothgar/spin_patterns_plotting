#include <iostream>
#include <string>
#include<map>

#include "TH2F.h"
class CrossCount
{
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

    struct data{
      float q;
      float pz;
      int clockcross;
    };

    struct spin_db{
      float pol_blue;
      float pol_yell;
      std::vector<int> bunch_pol;
    };
    // cross_count[dim][dim]
    // first dimension: arm: 0 = south, 1 = north
    // second dimension: charge: 0 = negative, 1 = positive
    TH2F* cross_count[2][2];
    TH2F* cross_count_sum;
    std::map<int,std::vector<data> > muons_;
    std::map<int,spin_db> spin_data_;
};
