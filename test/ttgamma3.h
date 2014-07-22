//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 26 14:11:59 2013 by ROOT version 5.32/00
// from TChain ggNtuplizer/EventTree/
//////////////////////////////////////////////////////////

#ifndef ttgamma3_h
#define ttgamma3_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TString.h>
#include <TProofOutputFile.h>

// STD libraries
#include <map>
#include <string>
#include <vector>

// Tree reader
#include"EventTree.h"

using namespace std;

// Header file for the classes stored in the TTree if any.

class ttgamma3 : public TSelector {
public :

  void              ParseInput();
  //void              CreateHistograms();
  TString           fMyOpt;
  int               fChannel;
  bool              fVerbose;
  bool              fIsMC;
  bool              fPUreweighting;
  TString           fSample;
  TString           fOutdir;
  bool              fdoSync;
  bool              fdoJER;
  bool              fdoJERdown;
  bool              fdoJERup;
  bool              fdoBTAG;
  bool              fdoBTAGdown;
  bool              fdoBTAGup;
  bool              fdoTOPPT;
  bool              fdoTOPPTdown;
  bool              fdoTOPPTup;
  bool              fdoHLT;
  bool              fdoSkim;
  bool              fdoMuSF;
  bool              fKeepOnlyPhotons;
  int               fN_tt_filter;
  TProofOutputFile *fProofFile; // For optimized merging of the ntuple
  TH1F             *h1test;
  TH1F             *hcutflow;
  TH1F             *hPU_weights;
  vector< string >  fCutLabels;
  map< string, double > cutmap;
  map< string, TH1*> hmuons;
  map< string, TH1*> helectrons;
  map< string, TH1*> hphotons;
  map< string, TH1*> hjets;
  map< string, TH1*> hbtag;
  map< string, TH1*> hPVs;
  map< string, TH1*> hMET;
  map< string, TH1*> hM;
  map< string, TH1*> hMC;
  void               WriteHistograms(const char* name, map<string, TH1*> hcontainer);
  TFile          *fFile;

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  EventTree      *fReader;

  int phoRegion( double absEta);
  double phoEffArea03ChHad( double phoEta);
  double phoEffArea03NeuHad( double phoEta);
  double phoEffArea03Pho( double phoEta);
  float SFtop( float pt );

  ttgamma3(TTree * /*tree*/ =0) : fProofFile(0),h1test(0),hPU_weights(0),fFile(0),fChain(0) 
  { 
    fChannel =       1; //default 2=e+jets, 1=mu+jets
    fVerbose =       false;
    fIsMC    =       false;
    fPUreweighting = false;
    fSample  =       "";
    fOutdir  =       "";
    fdoSync  =       false;
    fdoJER   =       true;
    fdoJERdown =     false;
    fdoJERup =       false;
    fdoHLT   =       true;
    fdoSkim  =       false;
    fdoMuSF  =       true;
    fdoBTAG  =       true;
    fdoBTAGdown =    false;
    fdoBTAGup   =    false;
    fdoTOPPT =       true;
    fdoTOPPTdown =   false;
    fdoTOPPTup   =   false;
    fKeepOnlyPhotons=false;
    fN_tt_filter =   0;
  }
   virtual ~ttgamma3() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(ttgamma3,0);
};

#endif

#ifdef ttgamma3_cxx
void ttgamma3::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fReader = new EventTree( tree );
   fReader->InitSkim();
   //fChain = tree;
   //fChain->SetMakeClass(1);
   //fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin);
}

Bool_t ttgamma3::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ttgamma3_cxx
