#include <iostream>
#include <map>
#include <string>
#include <vector>
//#include "getProof.C"
//#include "TDSet.h"
//#include "TProof.h"
#include "TChain.h"
#include "TList.h"
#include "TH1F.h"
#include "TString.h"
#include "ttgamma3.h"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class std::map<std::string,std::string>+;
#pragma link C++ class std::pair<std::string,std::string>+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::map<std::string,std::vector<std::string>>+;
#pragma link C++ class std::pair<std::string,std::vector<std::string>>+;
#pragma link C++ class std::map<std::string,TH1*>+;
#endif

//typedef std::map<std::string, std::string> TStrStrMap;
//typedef std::pair<std::string, std::string> StrPair;
typedef std::vector<std::string> StrVector;
typedef std::map<std::string, StrVector> StrVecMap;
typedef std::pair<std::string, StrVector> StrVecPair;


void runttgamma(TString sample="all", TString ExtraOpts= "", int workers=8)  
{

  //TStrStrMap msamples;
  StrVecMap vsamples;
  StrVector vec;

  std::string prefixData="/eos/uscms/store/user/iraklis/ggNtuples/";
  std::string prefixMC = "/eos/uscms/store/user/iraklis/ggNtuples/";

  //sync exercise
  vec.push_back("/eos/uscms/store/user/yumiceva/ttgamma/sync/sync.root");
  vsamples.insert(StrVecPair("sync",vec));
  vec.clear();

  //Wjets
  vec.push_back("/eos/uscms/store/user/makouski/job_summer12_Wjets.root");
  vsamples.insert(StrVecPair("Wjets",vec));
  vec.clear();
  //WW_2l2nu
  vec.push_back(prefixMC+"job_summer12_WW_2l2nu.root");
  vsamples.insert(StrVecPair("WW_2l2nu",vec));
  vec.clear();
  //WWg
  vec.push_back(prefixMC+"job_summer12_WWg.root");
  vsamples.insert(StrVecPair("WWg",vec));
  vec.clear();
  //WZ_2l2q
  vec.push_back(prefixMC+"job_summer12_WZ_2l2q.root");
  vsamples.insert(StrVecPair("WZ_2l2q",vec));
  vec.clear();
  //WZ_3lnu
  vec.push_back(prefixMC+"job_summer12_WZ_3lnu.root");
  vsamples.insert(StrVecPair("WZ_3lnu",vec));
  vec.clear();
  //Wg
  vec.push_back(prefixMC+"job_summer12_Wg.root");
  vsamples.insert(StrVecPair("Wg",vec));
  vec.clear();
  //Wgg_FSR
  vec.push_back(prefixMC+"job_summer12_Wgg_FSR.root");
  vsamples.insert(StrVecPair("Wgg_FSR",vec));
  vec.clear();
  //ZZ_2e2mu
  vec.push_back(prefixMC+"job_summer12_ZZ_2e2mu.root");
  vsamples.insert(StrVecPair("ZZ_2e2mu",vec));
  vec.clear();
  //ZZ_2e2tau
  vec.push_back(prefixMC+"job_summer12_ZZ_2e2tau.root");
  vsamples.insert(StrVecPair("ZZ_2e2tau",vec));
  vec.clear();
  //ZZ_2mu2tau
  vec.push_back(prefixMC+"job_summer12_ZZ_2mu2tau.root");
  vsamples.insert(StrVecPair("ZZ_2mu2tau",vec));
  vec.clear();
  //ZZ_4e
  vec.push_back(prefixMC+"job_summer12_ZZ_4e.root");
  vsamples.insert(StrVecPair("ZZ_4e",vec));
  vec.clear();
  //ZZ_4mu
  vec.push_back(prefixMC+"job_summer12_ZZ_4mu.root");
  vsamples.insert(StrVecPair("ZZ_4mu",vec));
  vec.clear();
  //ZZ_4tau
  vec.push_back(prefixMC+"job_summer12_ZZ_4tau.root");
  vsamples.insert(StrVecPair("WWg",vec));
  vec.clear();
  //Zg
  vec.push_back(prefixMC+"job_summer12_Zg.root");
  vsamples.insert(StrVecPair("Zg",vec));
  vec.clear();
  //diphoton_box_10to25
  vec.push_back(prefixMC+"job_summer12_diphoton_box_10to25.root");
  vsamples.insert(StrVecPair("diphoton_box_10to25",vec));
  vec.clear();
  //diphoton_box_25to250
  vec.push_back(prefixMC+"job_summer12_diphoton_box_25to250.root");
  vsamples.insert(StrVecPair("diphoton_box_25to250",vec));
  vec.clear();
  //diphoton_box_250toInf
  vec.push_back(prefixMC+"job_summer12_diphoton_box_250toInf.root");
  vsamples.insert(StrVecPair("diphoton_box_250toInf",vec));
  vec.clear();
  //t_s
  vec.push_back(prefixMC+"job_summer12_t_s.root");
  vsamples.insert(StrVecPair("t_s",vec));
  vec.clear();
  //t_t
  vec.push_back(prefixMC+"job_summer12_t_t.root");
  vsamples.insert(StrVecPair("t_t",vec));
  vec.clear();
  //t_tW
  vec.push_back(prefixMC+"job_summer12_t_tW.root");
  vsamples.insert(StrVecPair("t_tW",vec));
  vec.clear();
  //tbar_s
  vec.push_back(prefixMC+"job_summer12_tbar_s.root");
  vsamples.insert(StrVecPair("tbar_s",vec));
  vec.clear();
  //tbar_t
  vec.push_back(prefixMC+"job_summer12_tbar_t.root");
  vsamples.insert(StrVecPair("tbar_t",vec));
  vec.clear();
  //tbar_tW
  vec.push_back(prefixMC+"job_summer12_tbar_tW.root");
  vsamples.insert(StrVecPair("tbar_tW",vec));
  vec.clear();
  //ttW
  vec.push_back(prefixMC+"job_summer12_ttW.root");
  vsamples.insert(StrVecPair("ttW",vec));
  vec.clear();
  //ttZ
  vec.push_back(prefixMC+"job_summer12_ttZ.root");
  vsamples.insert(StrVecPair("ttZ",vec));
  vec.clear();
  //ttg
  vec.push_back(prefixMC+"job_summer12_ttg.root");
  vsamples.insert(StrVecPair("ttg",vec));
  vec.clear();
  //ttgWz
  vec.push_back(prefixMC+"job_summer12_ttg_WHIZARD_v11.root");
  vsamples.insert(StrVecPair("ttgWz",vec));
  vec.clear();
  //ttjets
  vec.push_back(prefixMC+"job_summer12_ttinclusive.root");
  vsamples.insert(StrVecPair("ttjets",vec));
  vec.clear();
  //ttjers_2l
  //vec.push_back(prefixMC+"job_summer12_ttjets_2l.root");
  //vsamples.insert(StrVecPair("ttjets_2l",vec));
  //vec.clear();
  //DATA
  // electrons
  vec.push_back(prefixData+"job_electron_2012a_Jan22rereco.root");
  vec.push_back(prefixData+"job_electron_2012b_Jan22rereco.root");
  vec.push_back(prefixData+"job_electron_2012c_Jan2012rereco.root");
  vec.push_back("/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/iraklis/Makouski/CMSSW_5_3_12/src/ggAnalysis/ggNtuplizer/test/job_electron_2012d.root");
  // muons
  //vec.push_back(prefixData+"job_muon_2012a_Jan22rereco.root");
  //vec.push_back(prefixData+"job_muon_2012c_Jan22rereco.root");

  vsamples.insert(StrVecPair("data",vec));
  vec.clear();

  //TProof *proof = getProof("lite://",nwrks);
  TString proofOpt(Form("workers=%i",workers));
  //TProof *proof = TProof::Open(proofOpt);
  //if (!proof) {
  //  Printf("runProof: could not start/attach a PROOF session");
  //  return;
  //}

  //proof->SetLogLevel(2); // verbose

  TChain *chain = new TChain("ggNtuplizer/EventTree");
  StrVecMap::iterator ite;//
  for ( ite = vsamples.begin(); ite != vsamples.end(); ++ite) {//
    
    std::string name(sample.Data());
    std::string tmpname = (*ite).first;
    if ( name != "all" && name != tmpname ) continue;
    if ( name == "all" ) name = tmpname;
    StrVector location = vsamples[name];
    cout<< "Input sample: "<< name << endl;
    //StrVecMap::iterator ite;
    //for (ite = location.begin(); ite != location.end(); ++ite) {
    for ( size_t i=0; i < location.size(); ++i ) {
    
      //chain->Add(location.c_str());
      chain->Add( location[i].c_str() );
    }
    cout << "List of files:" << endl;
    chain->ls();
    //chain->SetProof(); // to run in PROOF mode
    TString opts(Form("sample=%s",name.c_str()));
    if (ExtraOpts != "" ) opts = ExtraOpts + " " + opts;

    ttgamma3 *myselector = new ttgamma3();

    //chain->Process("ttgamma3.C",opts.Data() );
    chain->Process(myselector, opts.Data() );
    cout << "Process: ttgamma.C done" << endl;
  }
  // logs
  
  //TList *logList = proof->GetManager()->GetSessionLogs()->GetListOfLogs();
  //for (int i=1; i< logList->GetSize(); ++i)
  //  {
  //    TProffLogElem *logElem = ( TProofLogElem* ) logList->At( i );
  //    macro = logElem->GetMacro();
  //    macro->SaveSource("data_muons_"+TString(Form("%i",i))+".stdout");
  //  }
  

  //chain->SetProof(0);
  //chain->Delete();
  //delete chain;
  //proof->ClearInput();
  //proof->ClearData();
  //proof->ClearInputData();
  //delete proof;
  cout << "done"<<endl;
}

# ifndef __CINT__
int main(int argc, char** argv)
{

  TString sample = "all";
  TString options = "";
  if ( argc > 0 ) sample = argv[1];
  if ( argc > 1 ) 
    {
      for (int i=2; i<= argc; ++i)
        options += argv[i];
    }
  
  runttgamma(sample, options);
  return 0;
}
# endif
