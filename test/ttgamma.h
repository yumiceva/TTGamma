//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  6 20:19:50 2013 by ROOT version 5.32/00
// from TTree EventTree/Event data
// found on file: /eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012a_Aug6rereco_skim.root
//////////////////////////////////////////////////////////

#ifndef ttgamma_h
#define ttgamma_h

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

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ttgamma : public TSelector {

private:
  void              ParseInput();
  //void              CreateHistograms();
  TString           fMyOpt;
  int               fChannel;
  bool              fVerbose;
  bool              fIsMC;
  TString           fSample;
  TString           fOutdir;
  TProofOutputFile *fProofFile; // For optimized merging of the ntuple
  TH1F             *h1test;
  TH1F             *h2test;
  TH1F             *h3test;
  TH1F             *h4test;
  TH1F             *h5test;
  TH1F             *h6test;
  TH1F             *h7test;
  TH1F             *h8test;
  TH1F             *h9test;
  TH1F             *h10test;
  TH1F             *h11test;
  TH1F             *h12test;
  TH1F             *h13test;  
  TH1F             *h14test; 

public :

   TFile          *fFile;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[451];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[44][3];   //[nVtx]
   Int_t           IsVtxGood;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Int_t           nEle;
   Float_t         eleEcalEn[9];   //[nEle]
   Float_t         elePt[9];   //[nEle]
   Float_t         elePhi[9];   //[nEle]
   Float_t         eleSCEta[9];   //[nEle]
   Float_t         eleD0[9];   //[nEle]
   Float_t         eleDz[9];   //[nEle]
   Float_t         eleHoverE[9];   //[nEle]
   Float_t         eleEoverP[9];   //[nEle]
   Float_t         eledEtaAtVtx[9];   //[nEle]
   Float_t         eledPhiAtVtx[9];   //[nEle]
   Float_t         eleSigmaIEtaIEta[9];   //[nEle]
   Float_t         eleIsoTrkDR03[9];   //[nEle]
   Float_t         eleIsoEcalDR03[9];   //[nEle]
   Float_t         eleIsoHcalDR03[9];   //[nEle]
   Int_t           eleMissHits[9];   //[nEle]
   Float_t         eleIDMVATrig[9];   //[nEle]
   Float_t         elePFChIso03[9];   //[nEle]
   Float_t         elePFPhoIso03[9];   //[nEle]
   Float_t         elePFNeuIso03[9];   //[nEle]
   Int_t           nPho;
   Float_t         phoEt[9];   //[nPho]
   Float_t         phoEta[9];   //[nPho]
   Float_t         phoPhi[9];   //[nPho]
   Float_t         phoHoverE[9];   //[nPho]
   Float_t         phoSigmaIEtaIEta[9];   //[nPho]
   Float_t         phoPFChIso[9];   //[nPho]
   Float_t         phoPFPhoIso[9];   //[nPho]
   Float_t         phoPFNeuIso[9];   //[nPho]
   Int_t           phoIsConv[9];   //[nPho]
   Int_t           nMu;
   Float_t         muEta[10];   //[nMu]
   Float_t         muPhi[10];   //[nMu]
   Float_t         muPt[10];   //[nMu]
   Float_t         muPFIsoR04_CH[10];   //[nMu]
   Float_t         muPFIsoR04_NH[10];   //[nMu]
   Float_t         muPFIsoR04_Pho[10];   //[nMu]
   Float_t         rho25;
   Float_t         rho25_neu;
   Float_t         rho25_muPFiso;
   Float_t         rho25_elePFiso;
   Int_t           nJet;
   Float_t         jetPt[74];   //[nJet]
   Float_t         jetEta[74];   //[nJet]
   Float_t         jetPhi[74];   //[nJet]
   Float_t         jetCHF[74];   //[nJet]
   Float_t         jetNHF[74];   //[nJet]
   Float_t         jetCEF[74];   //[nJet]
   Float_t         jetNEF[74];   //[nJet]
   Int_t           jetNConstituents[74];   //[nJet]
   Float_t         jetCombinedSecondaryVtxBJetTags[74];   //[nJet]
   Float_t         jetNCharged[74];   //[nJet]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPFIsoR04_CH;   //!
   TBranch        *b_muPFIsoR04_NH;   //!
   TBranch        *b_muPFIsoR04_Pho;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho25_neu;   //!
   TBranch        *b_rho25_muPFiso;   //!
   TBranch        *b_rho25_elePFiso;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetCombinedSecondaryVtxBJetTags;   //!
   TBranch        *b_jetNCharged;   //!

  ttgamma(TTree * /*tree*/ =0) : fProofFile(0),h1test(0),fFile(0),fChain(0) 
  { 
    fChannel = 2; //default e+jets
    fVerbose = false;
    fIsMC    = false;
    fSample  = "";
    fOutdir  = "";
  }

   virtual ~ttgamma() { }
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

   ClassDef(ttgamma,0);
};

#endif

#ifdef ttgamma_cxx
void ttgamma::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   Info("Init","tree: %p", tree);
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleEcalEn", eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleIDMVATrig", eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("elePFChIso03", elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoPFChIso", phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muPFIsoR04_CH", muPFIsoR04_CH, &b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH", muPFIsoR04_NH, &b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho", muPFIsoR04_Pho, &b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho25_neu", &rho25_neu, &b_rho25_neu);
   fChain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso, &b_rho25_muPFiso);
   fChain->SetBranchAddress("rho25_elePFiso", &rho25_elePFiso, &b_rho25_elePFiso);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetNCharged", jetNCharged, &b_jetNCharged);
}

Bool_t ttgamma::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ttgamma_cxx
