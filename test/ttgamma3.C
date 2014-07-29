#define ttgamma3_cxx
// The class definition in ttgamma3.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("ttgamma3.C")
// Root > T->Process("ttgamma3.C","some options")
// Root > T->Process("ttgamma3.C+")
//

#include "ttgamma3.h"
#include "EvtCleaning.h"
#include "ElectronSelector.h"
#include "MuonSelector.h"
#include "MuonScaleFactor.h"
#include "BTagScaleFactor.h"

#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>


int ttgamma3::phoRegion(double absEta){
  int region = 0;
  if( absEta >= 1.0  ) region++;
  if( absEta >= 1.479) region++;
  if( absEta >= 2.0  ) region++;
  if( absEta >= 2.2  ) region++;
  if( absEta >= 2.3  ) region++;
  if( absEta >= 2.4  ) region++;
  return region;
}
double ttgamma3::phoEffArea03ChHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  return area[phoRegion(eta)];
}
double ttgamma3::phoEffArea03NeuHad(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  return area[phoRegion(eta)];
}
double ttgamma3::phoEffArea03Pho(double phoEta){
  double eta = TMath::Abs(phoEta);
  static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  return area[phoRegion(eta)];
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
float ttgamma3::SFtop(float pt){
  if(fdoTOPPT && fdoTOPPTdown==false && fdoTOPPTup==false) return exp(0.159 - 0.00141*pt);
  if(fdoTOPPT && fdoTOPPTup) return exp(0.148 - 0.00129*pt);
  return 1.0;
}

void ttgamma3::ParseInput()
{

  if (fMyOpt.Contains("muon") || fChannel == 1)
    {
      fChannel = 1;
      Info("Begin","Apply selection for MUON+jets events");
    }
  if (fMyOpt.Contains("electron") || fChannel == 2)
    {
      fChannel = 2;
      Info("Begin","Apply selection for ELECTRON+jets events");
    }
  if (fMyOpt.Contains("verbose"))
    {
      fVerbose = true;
    }
  if (fMyOpt.Contains("outdir")) {
    TString tmp = fMyOpt;
    tmp = tmp.Remove(0,fMyOpt.Index("outdir")+7);
    if (tmp.Index(" ") > -1 ) tmp = tmp.Remove(tmp.Index(" "),tmp.Sizeof());
    fOutdir = tmp+"/";
    Info("Begin","Output files will be written to directory: %s", fOutdir.Data());
  }
  if (fMyOpt.Contains("sample"))
    {
      TString tmp = fMyOpt;
      tmp = tmp.Remove(0,fMyOpt.Index("sample")+7);
      if (tmp.Index(" ") > -1 ) tmp = tmp.Remove(tmp.Index(" "),tmp.Sizeof());
      fSample = tmp;

      if (fMyOpt.Contains("onlyphotons") && 
          (fSample=="ttjets_0l" || fSample=="ttjets_1l" || fSample=="ttjets_2l"))
        {
          fSample = tmp+"_g";
          fKeepOnlyPhotons = true;
        }

      Info("Begin","Histogram names will have suffix: %s", fSample.Data());

      if ( fSample.Contains("data") )
        {
          fIsMC = false;
          Info("Begin","This sample is treated as DATA");
        }
      else
        {
          fIsMC = true;
          Info("Begin","This sample is treated as MC");
        }
    }

  if (fIsMC) 
    {
      fPUreweighting = true;
      fdoTOPPT = true;
      fdoBTAG = true;
    }
  else 
    {
      fPUreweighting = false;
      fdoMuSF = false;
      fdoTOPPT = false;
      fdoBTAG = false;
    }

  if (fMyOpt.Contains("sync"))
    {
      fdoSync = true;
      fdoJER = false;
      fdoHLT = false;
      fPUreweighting = false;
      fdoMuSF = false;
    }
  
  if (fMyOpt.Contains("skim") && !(fMyOpt.Contains("inputskim")) )
    {
      fdoSkim = true;
      Info("Begin","Running in SKIM mode");
      fPUreweighting = false;
      fdoJER = false;
    }

  if ( fMyOpt.Contains("noMuSF") ) fdoMuSF = false;
  
  if ( fPUreweighting ) Info("Begin","Apply PU reweighting.");
  if ( fdoMuSF ) Info("Begin","Apply muon scale factors.");
  if ( fdoTOPPT ) Info("Begin","Apply top pt reweighting.");
  if ( fdoBTAG ) Info("Begin","Apply b-tagging scale factors.");
}

void ttgamma3::WriteHistograms(const char* name, map<string, TH1*> hcontainer)
{
  fFile->cd();
  fFile->mkdir(name);
  fFile->cd(name);
  for ( map<string,TH1* >::const_iterator imap=hcontainer.begin(); imap!=hcontainer.end(); ++imap)
    {
      TH1 *temp = imap->second;
      if ( temp->GetEntries() > 0 )
        temp->Write();
    }
}

void ttgamma3::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   fMyOpt = option;
   Info("Begin", "Processing with options: %s", option.Data());

   ParseInput();
   //fVerbose = true;
   if (fPUreweighting) Info("Begin", "PU reweighting is ON");
}

void ttgamma3::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  fMyOpt = option;
  Info("SlaveBegin","Processing with options: %s", option.Data() );

  //ParseInput(); // Do we need this?

  //initialize the Tree branch addresses
  Init(tree);

  // get PU weights
  if (fPUreweighting) {
    TFile *PU_data_file     = TFile::Open( "HelperFiles/MyDataPileupHistogram.root" );
    if ( PU_data_file->IsZombie() ) {
      Info("SlaveBegin","Error opening file: %s", PU_data_file->GetName() );
    }
    hPU_weights = new TH1F( *(static_cast<TH1F*>(PU_data_file->Get( "pileup" )->Clone() )) );
    TFile *PU_MC_file = TFile::Open( "HelperFiles/MyMCPileupHistogram.root");
    if ( PU_MC_file->IsZombie() ) {
      Info("SlaveBegin","Error opening file: %s", PU_MC_file->GetName() );
    }
    TString hmc_name = "hPUTrue_";
    hmc_name += fSample;
    if (! PU_MC_file->GetListOfKeys()->Contains( hmc_name ) ) {
      Info("SlaveBegin","Error: Histogram name "+hmc_name+" does not exist in file "+PU_MC_file->GetName());
      assert("SlaveBegin STOP");
      
    }      
    TH1F* hPU_mc = new TH1F( *(static_cast<TH1F*>(PU_MC_file->Get( hmc_name )->Clone() )) );
    hPU_weights->Scale(1./(hPU_weights->Integral()));
    hPU_mc->Scale(1./(hPU_mc->Integral()));
    hPU_weights->Divide( hPU_mc );

  }

  // Create Output file
  TString tmpfilename = "results";
  TString outprefixfile = "results_";
  if ( fdoSkim ) {
    outprefixfile = "skim_";
  }
  if ( fSample != "" ) tmpfilename = fOutdir+outprefixfile+fSample+".root";
  else tmpfilename = fOutdir + outprefixfile+"data.root";
  fFile = TFile::Open( tmpfilename, "RECREATE");
  //Info("SlaveBegin", "Output filename: %s", tmpfilename );

  TString hname = "_"+fSample;
  // vertices
  hPVs["N_pre"] = new TH1F("NPV_pre"+hname,"Number of PVs",35,-0.5,34.5);
  hPVs["N"] = new TH1F("NPV"+hname,"Number of PVs",35,-0.5,34.5);
   // muons
  hmuons["pt_pre"] = new TH1F("muon_pt_pre"+hname,"p_{T}^{#mu} [GeV/c]", 50, 0, 500);
  hmuons["eta_pre"] = new TH1F("muon_eta_pre"+hname,"#eta^{#mu}", 20, -2.1, 2.1);
  hmuons["phi_pre"] = new TH1F("muon_phi_pre"+hname,"#phi^{#mu}", 30, -3.15, 3.15);
  hmuons["reliso_pre"] = new TH1F("muon_reliso_pre"+hname,"I_{rel}", 40, 0, 0.2);
  hmuons["relisocorr_pre"] = new TH1F("muon_relisocorr_pre"+hname,"I_{rel}^{corr}", 40, 0, 0.2);
  hmuons["pt"] = new TH1F("muon_pt"+hname,"p_{T}^{#mu} [GeV/c]", 50, 0, 500);
  hmuons["eta"] = new TH1F("muon_eta"+hname,"#eta^{#mu}", 20, -2.1, 2.1);
  hmuons["phi"] = new TH1F("muon_phi"+hname,"#phi^{#mu}", 30, -3.15, 3.15);
  hmuons["relisocorr"] = new TH1F("muon_relisocorr"+hname,"I_{rel}^{corr}", 40, 0, 0.2);
  // electrons
  helectrons["pt_pre"] = new TH1F("electron_pt_pre"+hname,"p_{T}^{e} [GeV/c]", 50, 20, 500);
  helectrons["eta_pre"] = new TH1F("electron_eta_pre"+hname,"#eta^{e}", 20, -2.1, 2.1);
  helectrons["phi_pre"] = new TH1F("electron_phi_pre"+hname,"#phi^{#mu}", 20, 0, 3.15);
  helectrons["reliso_pre"] = new TH1F("electron_reliso_pre"+hname,"I_{rel}", 40, 0, 0.2);
  helectrons["relisocorr_pre"] = new TH1F("electron_relisocorr_pre"+hname,"I_{rel}^{corr}", 40, 0, 0.2);
  // jets
  hjets["N_pre"] = new TH1F("jet_N_pre"+hname,"Number of jets",8,-0.5,7.5);
  hjets["pt_pre"] = new TH1F("jet_pt_pre"+hname,"jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt1_pre"] = new TH1F("jet_pt1_pre"+hname,"1st jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt2_pre"] = new TH1F("jet_pt2_pre"+hname,"2nd jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt3_pre"] = new TH1F("jet_pt3_pre"+hname,"3rd jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt4_pre"] = new TH1F("jet_pt4_pre"+hname,"4th jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["deltaR_lepton_pre"] = new TH1F("deltaR_lepton_pre"+hname,"#DeltaR(jet,lepton)", 50, 0, 7);
  hjets["N"] = new TH1F("jet_N"+hname,"Number of jets",8,-0.5,7.5);
  hjets["pt"] = new TH1F("jet_pt"+hname,"jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt1"] = new TH1F("jet_pt1"+hname,"1st jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt2"] = new TH1F("jet_pt2"+hname,"2nd jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt3"] = new TH1F("jet_pt3"+hname,"3rd jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["pt4"] = new TH1F("jet_pt4"+hname,"4th jet p_{T} [GeV/c]", 50, 30, 800);
  hjets["Nbtags_CSVM"] = new TH1F("Nbjets_CSVM"+hname,"Tagged b-jets",3,-0.5,2.5);
  hjets["btag_CSV"] = new TH1F("btag_CSV"+hname,"btag CSV discriminator",50,0,1);
  hjets["pt_top"]  = new TH1F("pt_top"+hname,"top p_{T} [GeV]",50,0,1500);
  // btagging
  hbtag["jet_pt_b"] = new  TH1F("jet_pt_b","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["jet_pt_c"] = new  TH1F("jet_pt_c","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["jet_pt_udsg"] = new  TH1F("jet_pt_udsg","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["jet_pt_CSVM_b"] = new  TH1F("jet_pt_CSVM_b","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["jet_pt_CSVM_c"] = new  TH1F("jet_pt_CSVM_c","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["jet_pt_CSVM_udsg"] = new  TH1F("jet_pt_CSVM_udsg","Jet p_{T} [GeV/c]", 50, 0, 200);
  hbtag["discriminator_CSV"] = new TH1F("discriminator_CSV", "Discriminator CSV", 1500, -20, 20);
  hbtag["discriminator_CSV_b"] = new TH1F("discriminator_CSV_b", "Discriminator CSV", 1500, -20, 20);
  hbtag["discriminator_CSV_c"] = new TH1F("discriminator_CSV_c", "Discriminator CSV", 1500, -20, 20);
  hbtag["discriminator_CSV_udsg"] = new TH1F("discriminator_CSV_udsg", "Discriminator CSV", 1500, -20, 20);
  // MET
  hMET["MET_pre"] = new TH1F("MET_pre"+hname,"Missing Transverse Energy [GeV]", 50, 0, 300);
  hMET["MET"] = new TH1F("MET"+hname,"Missing Transverse Energy [GeV]", 50, 0, 300);
  // photons
  hphotons["mc_momID"] = new TH1F("photon_MCmomID","MC photon mother ID",51,-25.5,25.5);
  hphotons["mc_momIDoverlap"] = new TH1F("photon_MCmomIDoverlap","MC photon mother ID",51,-25.5,25.5);
  hphotons["N"] = new TH1F("photon_N"+hname,"Number of photons",8,-0.5,7.5);
  hphotons["cut0_et"] = new TH1F("photon_cut0_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  hphotons["cut0_eta"] = new TH1F("photon_cut0_eta"+hname,"#eta^{#gamma}", 20, -2.1, 2.1);
  hphotons["cut0_phi"] = new TH1F("photon_cut0_phi"+hname,"#phi^{#gamma}", 20, 0, 3.15);
  //hphotons["cut1_pt"] = new TH1F("photon_cut1_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  //hphotons["cut1_eta"] = new TH1F("photon_cut1_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  //hphotons["cut1_phi"] = new TH1F("photon_cut1_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  hphotons["cut1_eVeto"] = new TH1F("photon_cut1_eVeto"+hname,"electron veto",2,-0.5,1.5);
  hphotons["cut1_hasConvTrk"] = new TH1F("photon_cut1_hasConvTrk"+hname,"HasConversionTrk",2,-0.5,1.5);
  hphotons["cut2_HoverE"] = new TH1F("photon_cut2_HoverE"+hname,"H/E",20,0,1);
  hphotons["cut3_sigmaietaieta"] = new TH1F("photon_cut3_sigmaietaieta"+hname,"#sigma_{i#etai#eta}",50,0,0.03);
  hphotons["cut4_chHadIso"] = new TH1F("photon_cut4_chHadIso"+hname,"Charged Hadron Isolation",40,0,10);
  hphotons["cut5_ntHadIso"] = new TH1F("photon_cut5_ntHadIso"+hname,"Neutral Hadron Isolation",40,0,10);
  hphotons["cut6_pfIso"] = new TH1F("photon_cut6_pfIso"+hname,"PF photon Isolation",40,0,10);
  hphotons["sigmaietaieta"] = new TH1F("photon_sigmaietaieta"+hname,"#sigma_{i#etai#eta}",50,0,0.03);
  hphotons["chHadIso"] = new TH1F("photon_chHadIso"+hname,"Charged Hadron Isolation",40,0,10);
  hphotons["Ngood"] = new TH1F("photon_Ngood"+hname,"Number of photons",8,-0.5,7.5);
  hphotons["et"] = new TH1F("photon_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20,500);
  hphotons["eta"] = new TH1F("photon_eta"+hname,"#eta^{#gamma}", 20, -2.1, 2.1);
  hphotons["phi"] = new TH1F("photon_phi"+hname,"#phi^{#gamma}", 20, 0, 3.15);
  hphotons["pfIso"] = new TH1F("photon_pfIso"+hname,"PF photon Isolation",40,0,10);
  hphotons["ntHadIso"] = new TH1F("photon_ntHadIso"+hname,"Neutral Hadron Isolation",40,0,10);
  hphotons["SCFRChIso"] = new TH1F("photon_SCFRChIso"+hname,"SCFR Charged Hadron Isolation",100,-0.5,20);
  hphotons["RandIso"] = new TH1F("photon_RandIso"+hname,"Random Cone Isolation",100,-0.5,20);
  hphotons["temp_sigmaietaieta"] = new TH1F("photon_temp_sigmaietaieta"+hname,"#sigma_{i#etai#eta}",50,0,0.03);
  hphotons["temp_SCFRChIso"] = new TH1F("photon_temp_SCFRChIso"+hname,"SCFR Charged Hadron Isolation",100,-0.5,20);
  //hphotons["temp2D_sigmaVsSCFRCh"] = new TH2F("photon_temp_sigmaVsSCFRCh"+hname,"#sigma_{i#etai#eta} vs ",50,0,0.03,100,-0.5,20); 
  // mass
  hM["M3"] = new TH1F("M3"+hname,"M3 [GeV/c^{2}]", 100, 0, 1000); 
  hM["WMt"] = new TH1F("Mt"+hname,"M_{T}(W) [GeV/c^{2}]", 50, 0, 300);
  // MC
  hMC["PID"] = new TH1F("MC_ID","MC ID",51,-25.5,25.5);

  map<string,TH1* > allhistos = hmuons;
  allhistos.insert( helectrons.begin(), helectrons.end() );
  allhistos.insert( hphotons.begin(), hphotons.end() );
  allhistos.insert( hmuons.begin(), hmuons.end() );
  allhistos.insert( hjets.begin(), hjets.end() );
  allhistos.insert( hbtag.begin(), hbtag.end() );
  allhistos.insert( hMET.begin(), hMET.end() );
  allhistos.insert( hM.begin(), hM.end() );
  allhistos.insert( hMC.begin(), hMC.end() );

  for ( map<string,TH1* >::const_iterator imap = allhistos.begin(); imap!=allhistos.end(); ++imap )
    {
      TH1 *temp = imap->second;
      temp->Sumw2();
      temp->SetXTitle( temp->GetTitle() );
    }
  
   h1test = new TH1F("h1test","Electron p_{T}",100,10.,400);
   h1test->SetDirectory(fFile);   

   // cut flow
   if (fChannel==1)
     { //muon +jets
       fCutLabels.push_back("Processed");
       fCutLabels.push_back("Cleaning");
       fCutLabels.push_back("HLT");
       fCutLabels.push_back("GoodPV");
       fCutLabels.push_back("OneIsoMuon");
       fCutLabels.push_back("LooseMuVeto");
       fCutLabels.push_back("LooseEleVeto");
       fCutLabels.push_back("1Jets");
       fCutLabels.push_back("2Jets");
       fCutLabels.push_back("3Jets");
       fCutLabels.push_back("4Jets");
       fCutLabels.push_back("Onebtag");
       fCutLabels.push_back("PhotonFiducial");
       fCutLabels.push_back("PhotonEleVeto");
       fCutLabels.push_back("PhotonHoverE");
       fCutLabels.push_back("PhotonSigmaieta");
       fCutLabels.push_back("PhotonchHadIso");
       fCutLabels.push_back("PhotonntHadIso");
       fCutLabels.push_back("PhotonIso");
       fCutLabels.push_back("PhotonDeltaRLepton");
       fCutLabels.push_back("PhotonDeltaRJet");
       fCutLabels.push_back("OnePhoton");
     }
   else if (fChannel==2)
     { //electron+jets
       fCutLabels.push_back("Processed");
       fCutLabels.push_back("Cleaning");
       fCutLabels.push_back("HLT");
       fCutLabels.push_back("GoodPV");
       fCutLabels.push_back("OneIsoEle");
       fCutLabels.push_back("LooseMuVeto");
       fCutLabels.push_back("DileptonVeto");
       fCutLabels.push_back("ConvRejection");
       fCutLabels.push_back("1Jets");
       fCutLabels.push_back("2Jets");
       fCutLabels.push_back("3Jets");
       fCutLabels.push_back("4Jets");
       fCutLabels.push_back("Onebtag");
       fCutLabels.push_back("PhotonFiducial");
       fCutLabels.push_back("PhotonEleVeto");
       fCutLabels.push_back("PhotonHoverE");
       fCutLabels.push_back("PhotonSigmaieta");
       fCutLabels.push_back("PhotonchHadIso");
       fCutLabels.push_back("PhotonntHadIso");
       fCutLabels.push_back("PhotonIso");
       fCutLabels.push_back("PhotonDeltaRLepton");
       fCutLabels.push_back("PhotonDeltaRJet");
       fCutLabels.push_back("OnePhoton");
     }

   hcutflow = new TH1F("cutflow","cut flow", fCutLabels.size(), 0.5, fCutLabels.size() +0.5 );

   for ( vector<string>::const_iterator ivec= fCutLabels.begin(); ivec!=fCutLabels.end(); ++ivec)
     {
       cutmap[ *ivec ] = 0;
     }

   //SKIM
   fReader->SetDirectory( fFile );
}

Bool_t ttgamma3::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ttgamma3::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  double EvtWeight = 1.;
  TLorentzVector p4muon, p4ele, p4lepton, p4MET;
  TLorentzVector p4photon;
  TLorentzVector p4Nu, p4OtherNu;
  vector< TLorentzVector > p4jets;
  vector< bool > vec_btags;
  
  if (fVerbose) cout << "Processing entry: " << entry << "  ====================="<<endl;
  else if ( entry%1000 == 0 )
    Info("Process","Entries read: %i", entry);

  //fChain->GetEntry(entry);
  fReader->GetEntry(entry);

  // MC plots for ttbar samples
  if ( fSample == "ttg" || fSample == "ttgWhizard" || fSample == "ttjets_0l" || fSample == "ttjets_1l" || fSample == "ttjets_2l" || fSample == "ttjets_0l_g" || fSample == "ttjets_1l_g" || fSample == "ttjets_2l_g")
    {
      for(int imc = 0; imc < fReader->nMC; ++imc)             //Loop over gen particles
        {
          hMC["PID"]->Fill( fReader->mcPID->at(imc) );
        
          ///////////////////////////////////
          // check photon mom in ttg sample
          ///////////////////////////////////
          if ( fSample == "ttg" || fSample == "ttgWhizard" ) 
            {
              // photons
              if ( fReader->mcPID->at(imc) == 22 && fReader->mcStatus->at(imc) == 1)
                {
                  hphotons["mc_momID"]->Fill( fReader->mcMomPID->at(imc) );
                }
            }
        }
    }

  ////////////////////////////////////////////////
  // Remove overlapping photon events in ttjets MG
  ////////////////////////////////////////////////
  if ( fSample == "ttjets_1l" || fSample == "ttjets_2l" || fSample == "ttjets_0l" || fSample == "ttjets_0l_g" || fSample == "ttjets_1l_g" || fSample == "ttjets_2l_g" )
    {
      if (fVerbose) cout << "Checking if event has overlapping photons" << endl;
 
      int mcNphotons = 0;
      int mcNphotons_overlap = 0;
      
      TLorentzVector p4bQuark;
      TLorentzVector p4bbarQuark;
      vector < TLorentzVector > p4MCPhotonsVec;
      vector < int > momPhotonIDVec;
      TLorentzVector tmpp4;

      for(int imc = 0; imc < fReader->nMC; ++imc)             //Loop over gen particles
        {
          // photons
          if ( fReader->mcPID->at(imc) == 22 && fReader->mcStatus->at(imc) == 1)
            {
              tmpp4.SetPtEtaPhiE( fReader->mcPt->at(imc), fReader->mcEta->at(imc), fReader->mcPhi->at(imc), fReader->mcE->at(imc) );
              hphotons["mc_momID"]->Fill( fReader->mcMomPID->at(imc) );
              
              //if ( fabs(fReader->mcMomPID->at(imc)) == 5 ||
              //     fabs(fReader->mcMomPID->at(imc)) == 24 )
              if (fReader->mcParentage->at(imc)==2 || (fReader->mcParentage->at(imc)==10 && fabs(fReader->mcMomPID->at(imc))==24) )
                // photon mom can only be from b or W
                {
                  p4MCPhotonsVec.push_back( tmpp4 );
                  momPhotonIDVec.push_back( fReader->mcMomPID->at(imc) );
                  mcNphotons++;
                }
            }
          if ( fabs(fReader->mcPID->at(imc)) == 5 && fabs(fReader->mcMomPID->at(imc)) == 6 ) // b from tops
            {
              if (fReader->mcPID->at(imc) == 5) p4bQuark.SetPtEtaPhiE( fReader->mcPt->at(imc), fReader->mcEta->at(imc), fReader->mcPhi->at(imc), fReader->mcE->at(imc) );
              if (fReader->mcPID->at(imc) ==-5) p4bbarQuark.SetPtEtaPhiE( fReader->mcPt->at(imc), fReader->mcEta->at(imc), fReader->mcPhi->at(imc), fReader->mcE->at(imc) );
            }
        }
      
      int ipho = 0;
      for ( vector<TLorentzVector>::const_iterator ivec= p4MCPhotonsVec.begin(); ivec != p4MCPhotonsVec.end(); ++ivec )
        {
          tmpp4 = *ivec;
          //hphotons["mc_momID"]->Fill( momPhotonIDVec[ipho] );
          if ( tmpp4.Et() > 20 && tmpp4.DeltaR( p4bQuark ) > 0.1 && tmpp4.DeltaR( p4bbarQuark ) > 0.1 )
            {
              mcNphotons_overlap++;
              hphotons["mc_momIDoverlap"]->Fill( momPhotonIDVec[ipho] );
            }
          //else 
          //  {
          //    hphotons["mc_momIDpass"]->Fill( momPhotonIDVec[ipho] );
          //  }
          ipho++;
        }
      p4MCPhotonsVec.clear();

      if ( mcNphotons_overlap > 0 && fKeepOnlyPhotons == false )
        {
          if (fVerbose) cout << "overlap photon, skip event" << endl;
          fN_tt_filter++;
          return kTRUE;
        }
      if ( mcNphotons_overlap == 0 && fKeepOnlyPhotons == true )
        {
          if (fVerbose) cout << "no photon found, skip event" << endl;
          return kTRUE;
        }
          
      if (fVerbose) cout << "event kept." << endl;
    }

  ////////////////////////////////////////
  // top pt reweighting
  ///////////////////////////////////////
  float gen_toppt=0.0;
  float gen_antitoppt=0.0;
  float toppt_weight = 1.0;
  for(int imc = 0; imc < fReader->nMC; ++imc)             //Loop over gen particles
    {
      if ( fReader->mcPID->at(imc) == 6 ) gen_toppt = fReader->mcPt->at(imc);
      if ( fReader->mcPID->at(imc) == -6 ) gen_antitoppt = fReader->mcPt->at(imc);
    }
  
  if( gen_toppt > 0.001 && gen_antitoppt > 0.001)
    toppt_weight = sqrt( SFtop(gen_toppt) * SFtop(gen_antitoppt) );

  if (fdoTOPPTup && fdoTOPPT) { EvtWeight = EvtWeight*toppt_weight*toppt_weight; }
  else if ( fdoTOPPT && fdoTOPPTup==false && fdoTOPPTdown==false ) { EvtWeight = EvtWeight*toppt_weight; }
  // fdoTOPPTdown is when no top pt reweighting is applied!
  

  ////////////////////////////////////////
  // get PU weight for MC samples
  ///////////////////////////////////////
  if (fPUreweighting) {
    // for BX 0 only
    float numInter = fReader->puTrue->at(1);
    EvtWeight = EvtWeight*hPU_weights->GetBinContent( hPU_weights->FindBin( numInter )  );
    if (fVerbose) cout << " PU weight = " << EvtWeight << endl;
  }

  cutmap["Processed"] += EvtWeight;

  /////////////////////////////////////////
  // Event Clean up
  ////////////////////////////////////////

  EvtCleaning EventFilter = EvtCleaning();
  EventFilter.Init(fReader, fVerbose );

  if ( EventFilter.Pass() )
    cutmap["Cleaning"] += EvtWeight;
  else
    return kTRUE;

  /////////////////////////////////////////
  //  HLT
  /////////////////////////////////////////
  
  if ( fdoHLT )
    {

      if ( fChannel==1 && fReader->HLTIndex[18] == 0 ) //HLT_IsoMu24_eta2p1_v
        { return kTRUE; }
      if ( fChannel==2 && fReader->HLTIndex[17] == 0 )
        { return kTRUE; } // 17 --> HLT_Ele27_WP80_v Trigger
    }

  cutmap["HLT"] += EvtWeight;
  if (fVerbose) cout << "Pass HLT" << endl;
  ////////////////////////////////////
  // PRIMARY VERTICES
  ///////////////////////////////////
  // check if first vertex is a good vtx
  // required NDF >= 4 && fabs(PVz) <= 24 && fabs(rho) <= 2 
  if ( fReader->IsVtxGood != 0 ) return kTRUE;
  //if ( fReader->vtxNTrk->at(0) >= 4 &&
  //fabs( fReader->vtx_z->at(0) ) <= 24 &&
  //     fabs( fReader->vtxD0->at(0) ) <= 2 )
  //  {
  cutmap["GoodPV"] += EvtWeight;
  //hPVs["N"]->Fill( fReader->nVtx , EvtWeight );
  
  if (fVerbose) cout << "Pass good vertex" << endl;

  //////////////////////////////////
  // ELECTRONS
  //////////////////////////////////
  if (fVerbose) cout << "Total number of electrons nEle = "<< fReader->nEle << endl;

  int Ngood_Ele = 0;
  int Nloose_Ele = 0;
  float relIso_Ele = -1;
  float relIsocorr_Ele = -1;
  bool passconversionveto = false;

  ElectronSelector eSelector = ElectronSelector();
  eSelector.Init(fReader, fVerbose, fIsMC );

  for(int ie = 0; ie < fReader->nEle; ++ie)		//Loop over the electrons in a event
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( fReader->elePt->at(ie),
                          fReader->eleSCEta->at(ie),
                          fReader->elePhi->at(ie),
                          fReader->eleEcalEn->at(ie) );

      h1test->Fill( tmpp4.Pt() );

      // loose electrons
      if ( eSelector.PassLoose( ie ) )
        {
          float relIso = eSelector.get_relIso();
          float relIsocorr = eSelector.get_relIsocorr();
          Nloose_Ele++;

          // good tight electrons
          if ( eSelector.PassTight(ie) )
            {
              // leading electron
              if ( Ngood_Ele == 0 )
                {
                  p4ele.SetPtEtaPhiE( tmpp4.Pt(), tmpp4.Eta(), tmpp4.Phi(), tmpp4.E() );
                  relIso_Ele = relIso;
                  relIsocorr_Ele = relIsocorr;
                  passconversionveto = (fReader->eleConvVtxFit->at(ie) == 0);
                  p4lepton = p4ele;
                }
              Ngood_Ele++;
          //helectrons["pt"]->Fill( tmpp4.Pt() );
            }
        } //loose

    } // end electrons
    
  //////////////////////////
  // MUONS
  /////////////////////////
  if (fVerbose) cout << "Total number of muons nMu = "<< fReader->nMu << endl;

  int Ngood_Mu = 0;
  int Nloose_Mu = 0;
  float relIso_Mu = -1;
  float relIsocorr_Mu = -1;
  MuonSelector muSelector = MuonSelector();
  muSelector.Init(fReader, fVerbose, fIsMC );

  for(int imu = 0; imu < fReader->nMu; ++imu)             //Loop over the muons in a event
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( fReader->muPt->at(imu),
                          fReader->muEta->at(imu),
                          fReader->muPhi->at(imu),
                          sqrt(fReader->muPt->at(imu)*fReader->muPt->at(imu) + fReader->muPz->at(imu)*fReader->muPz->at(imu)) );

      // loose muon
      if ( muSelector.PassLoose(imu) ) 
        {

          float relIso = muSelector.get_relIso();
          float relIsocorr = eSelector.get_relIsocorr(); 

          Nloose_Mu++;

          // tight muons
          if ( muSelector.PassTight(imu) )
            {
              // leading muon
              if ( Ngood_Mu == 0 )
                {
                  p4muon.SetPtEtaPhiE( tmpp4.Pt(), tmpp4.Eta(), tmpp4.Phi(), tmpp4.E() );
                  relIso_Mu = relIso;
                  relIsocorr_Mu = relIsocorr;
                  p4lepton = p4muon;
                }
              
              Ngood_Mu++;
            }
        }
    } // end muons

  ////////////////////////////
  // muon+jets selection
  ///////////////////////////
  if ( fChannel==1 ) 
    {

      if ( Ngood_Mu == 1 )
        {
          // Get SF for muon
          MuonScaleFactor muSF = MuonScaleFactor();
          muSF.Init();
          float muon_sf = 1.0;
          muon_sf = muSF.GetSF( p4lepton.Eta() );
          if (fVerbose) cout << "muon SF = " << muon_sf << endl;
          EvtWeight *= muon_sf;

          cutmap["OneIsoMuon"] += EvtWeight;
        }
      else 
        return kTRUE;

      if ( Nloose_Mu >= 2 )
        return kTRUE;
      else
        cutmap["LooseMuVeto"] += EvtWeight;
        
      if ( Nloose_Ele >= 1 )
        return kTRUE;
      else
        cutmap["LooseEleVeto"] += EvtWeight;
    }
  ////////////////////////////
  // electron+jets selection
  ///////////////////////////
  if ( fChannel==2 )
    {

      if ( Ngood_Ele == 1 )
        cutmap["OneIsoEle"] += EvtWeight;
      else
        return kTRUE;

      if ( Nloose_Mu >= 1 )
        return kTRUE;
      else
        cutmap["LooseMuVeto"] += EvtWeight;

      if ( Nloose_Ele >= 2 )
        return kTRUE;
      else
        cutmap["DileptonVeto"] += EvtWeight;

      if ( passconversionveto )
        cutmap["ConvRejection"] += EvtWeight;
      else
        return kTRUE;
    }


  ////////////////////////////////
  // JETS
  ///////////////////////////////
  if (fVerbose) cout << "Total number of jets nJet = "<< fReader->nJet << endl;

  int Ngood_Jets = 0;
  float met_x = 0.0;
  float met_y = 0.0;
  float CSVM = 0.679;
  int Nbtags = 0;
  bool pass1stJet, pass2ndJet, pass3rdJet, pass4thJet;
  pass1stJet = pass2ndJet = pass3rdJet = pass4thJet = false;

  double the_btag_weight = 1.0;
  BTagScaleFactor btagSF = BTagScaleFactor();

  for (int ij= 0; ij < fReader->nJet; ++ij)
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( fReader->jetPt->at(ij), fReader->jetEta->at(ij), fReader->jetPhi->at(ij), fReader->jetEn->at(ij) );

      float uncorr_jet_pt = fReader->jetRawPt->at(ij);
      met_x += uncorr_jet_pt*TMath::Cos( fReader->jetPhi->at(ij) );
      met_y += uncorr_jet_pt*TMath::Sin( fReader->jetPhi->at(ij) );

      double jet_pt = fReader->jetPt->at(ij);
      double jet_eta = fReader->jetEta->at(ij);
      
      // Apply JER smearing in MC
      if ( fIsMC && fdoJER && fReader->jetPt->at(ij) > 10 )
        {
          float gen_pt = fReader->jetGenJetPt->at(ij);
          if ( gen_pt < 0 ) continue;        
          // factor is (c - 1), where c are the eta-dependent scale factors, to be taken from the official twiki
          float c = 0;
          // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          // nominal
          if ( fabs(jet_eta) > 0 && fabs(jet_eta) < 0.5 )  c = 1.052;
          if ( fabs(jet_eta) > 0.5 && fabs(jet_eta) < 1.1) c = 1.057;
          if ( fabs(jet_eta) > 1.1 && fabs(jet_eta) < 1.7) c = 1.096;
          if ( fabs(jet_eta) > 1.7 && fabs(jet_eta) < 2.3) c = 1.134;
          if ( fabs(jet_eta) > 2.3 && fabs(jet_eta) < 5.0) c = 1.288;
          // minus sigma
          if ( fdoJERdown) {
            if ( fabs(jet_eta) > 0 && fabs(jet_eta) < 0.5 )  c = 0.990;
            if ( fabs(jet_eta) > 0.5 && fabs(jet_eta) < 1.1) c = 1.001;
            if ( fabs(jet_eta) > 1.1 && fabs(jet_eta) < 1.7) c = 1.032;
            if ( fabs(jet_eta) > 1.7 && fabs(jet_eta) < 2.3) c = 1.042;
            if ( fabs(jet_eta) > 2.3 && fabs(jet_eta) < 5.0) c = 1.089;
          }
          // plus sigma
          else if (fdoJERup) {
            if ( fabs(jet_eta) > 0 && fabs(jet_eta) < 0.5 )  c = 1.115;
            if ( fabs(jet_eta) > 0.5 && fabs(jet_eta) < 1.1) c = 1.114;
            if ( fabs(jet_eta) > 1.1 && fabs(jet_eta) < 1.7) c = 1.161;
            if ( fabs(jet_eta) > 1.7 && fabs(jet_eta) < 2.3) c = 1.228;
            if ( fabs(jet_eta) > 2.3 && fabs(jet_eta) < 5.0) c = 1.488;
          }
          float factor = (c-1);
          float deltapt = (jet_pt - gen_pt) * factor;
          float ptscale = fmax(0.0, (jet_pt + deltapt) / jet_pt);
          tmpp4 *= ptscale;
          met_x -= uncorr_jet_pt*TMath::Cos( fReader->jetPhi->at(ij) ) * ptscale;
          met_y -= uncorr_jet_pt*TMath::Sin( fReader->jetPhi->at(ij) ) * ptscale;
        }
      
      if ( tmpp4.Pt() > 30 &&
           fabs( tmpp4.Eta() ) < 2.4 &&
           fReader->jetPFLooseId->at(ij) == true )
        {
          
          float tmpDeltaR = tmpp4.DeltaR(p4lepton);
          
          if (fVerbose) cout << "Ngood_Jets = " << Ngood_Jets << " Pt = " << tmpp4.Pt() << " DeltaR(lepton)= " << tmpDeltaR << endl;

          if ( tmpDeltaR < 0.3 ) continue;
          //if ( aDeltaR <= jminDeltaR ) jminDeltaR = aDeltaR;

          float bdiscriminator = fReader->jetCombinedSecondaryVtxBJetTags->at(ij);
          //if (fVerbose) cout << "bdiscriminator = " << bdiscriminator << endl;
          //if (fVerbose) cout <<"jetPartonID = " << jetPartonID << endl;

          if (fIsMC && bdiscriminator > CSVM)
            {
              int jetPartonID = fReader->jetPartonID->at(ij);
              float SFb = 0.0;
              if ( fabs(jetPartonID) == 5 ) SFb = btagSF.bJetSF(jet_pt);
              else if ( fabs(jetPartonID) == 4) SFb = btagSF.cJetSF(jet_pt);
              //else if ( fabs(jetPartonID) == 1 ||
              //          fabs(jetPartonID) == 2 ||
              //          fabs(jetPartonID) == 3 ||
              //          fabs(jetPartonID) == 21 ) SFb = btagSF.lfJetSF(jet_pt, jet_eta);
              else SFb = btagSF.lfJetSF(jet_pt, jet_eta);

              the_btag_weight *= 1.0 - SFb;
            }

          // store jets and btags
          p4jets.push_back( tmpp4 );
          if ( bdiscriminator > CSVM )
            {
              //vec_btags.push_back( true );
              Nbtags++;
            }
          //else
          //vec_btags.push_back( false );

          if (Ngood_Jets == 0 && tmpp4.Pt() > 55.0 )
            {
              pass1stJet = true;

              hjets["deltaR_lepton_pre"]->Fill(tmpDeltaR);
            }

          if (Ngood_Jets == 1 && tmpp4.Pt() > 45.0 )
            {
              pass2ndJet = true;
            }

          if (Ngood_Jets == 2 && tmpp4.Pt() > 35.0 )
            {
              pass3rdJet = true;
            }

          if (Ngood_Jets == 3 && tmpp4.Pt() > 20.0 )// mim jet pt is 30, different from top ref. sel.
            {
              pass4thJet = true;
            }
          
          Ngood_Jets++;
        }

    } // end jets

  if (fVerbose) cout << "Number of good jets = " << Ngood_Jets << endl;
  if (fVerbose) cout << "Pass jets = " << pass1stJet << " " << pass2ndJet << " " << pass3rdJet << " " << pass4thJet << endl;

  
  // for the top ref. selection:
  /*
  if ( pass1stJet  )
    {
      cutmap["1Jets"] += EvtWeight;
      if ( pass2ndJet )
        {
          cutmap["2Jets"] += EvtWeight;
          if ( pass3rdJet )
            {
              cutmap["3Jets"] += EvtWeight;
              if ( pass4thJet )
                {
                cutmap["4Jets"] += EvtWeight;
                } else return kTRUE;
            } else return kTRUE;
        } else return kTRUE;
    } else return kTRUE;
  */

  //h1test->Fill( p4lepton.Pt() ); // leading pT lepton

  ////////////////////////
  // Jet selection
  if ( Ngood_Jets > 3 )
    {
      cutmap["4Jets"] += EvtWeight;
    } else return kTRUE;

  ////////////////////////
  // b-tagging selection
  ///////////////////////
  if (fVerbose) cout << "b-tagging weight ( 0 = no btags ) = " << (1. - the_btag_weight) << endl;
  if ( fdoBTAG && (1. - the_btag_weight) > 0 ) EvtWeight = EvtWeight*(1. - the_btag_weight);

  if ( Nbtags >= 1 )
    {
      cutmap["Onebtag"] += EvtWeight;
    }
  //else
  //  return kTRUE;

  ////////////////////////
  // SKIM Output
  ////////////////////////

  if (fdoSkim) 
    {
      // skim with at least 4 good jets
      //if ( pass1stJet && pass2ndJet && pass3rdJet && pass4thJet )
      if ( Ngood_Jets > 3 )
        {
          fReader->FillSkim( fFile );
          if (entry%1000 == 1) fReader->AutoSave();
        }
      return kTRUE;
    }

  ////////////////////////
  // pre-selection plots
  ///////////////////////
  bool passPreSel = false;
  //if ( pass1stJet && pass2ndJet && pass3rdJet && pass4thJet &&
  if ( Ngood_Jets > 3 &&
       Nbtags >= 1 )
    passPreSel = true;

  if ( ! passPreSel ) return kTRUE;
  if (fVerbose) cout << "pass pre-selection" << endl;

  // Compute M3
  float maxPt = -1.0;
  float  M3 = -1.0 ;
  TLorentzVector M3p4;

  for ( int ij = 0; ij < p4jets.size(); ++ij)
    {
      TLorentzVector tmpp4_1 = p4jets[ij];
      for ( int kj= ij+1; kj != p4jets.size(); ++kj )
        {
          TLorentzVector tmpp4_2 = p4jets[kj];
          for ( int nj= kj+1; nj != p4jets.size(); ++nj )
            {
              TLorentzVector tmpp4_3 = p4jets[nj];
              
              M3p4 = tmpp4_1 + tmpp4_2 + tmpp4_3;

              float sumPt = M3p4.Pt();
              
              if (sumPt > maxPt )
                {
                  maxPt = sumPt;
                  M3 = M3p4.M();
                }
            }
        }
    }

  //hPVs["N_pre"]->Fill( fReader->nVtx , EvtWeight );
  if (fChannel==1)
    {
      hmuons["pt_pre"]->Fill( p4lepton.Pt(), EvtWeight );
      hmuons["eta_pre"]->Fill( p4lepton.Eta(), EvtWeight );
      hmuons["phi_pre"]->Fill( p4lepton.Phi(), EvtWeight );
      hmuons["relisocorr_pre"]->Fill( relIsocorr_Mu, EvtWeight );
    }
  else
    {
      helectrons["pt_pre"]->Fill( p4lepton.Pt(), EvtWeight );
      helectrons["eta_pre"]->Fill( p4lepton.Eta(), EvtWeight );
      helectrons["phi_pre"]->Fill( p4lepton.Phi(), EvtWeight );
      helectrons["relisocorr_pre"]->Fill( relIsocorr_Ele, EvtWeight );
    }
  hPVs["N_pre"]->Fill( fReader->nVtx , EvtWeight );
  hjets["N_pre"]->Fill( Ngood_Jets, EvtWeight );
  hjets["pt1_pre"]->Fill( p4jets[0].Pt(), EvtWeight );
  hMET["MET_pre"]->Fill( fReader->pfMET, EvtWeight );
  hM["M3"]->Fill( M3, EvtWeight);

  if (fVerbose) cout << " pre-selection done." << endl;

  ///////////////////////////////////
  // PHOTONS
  ///////////////////////////////////
  float Ngamma[8] = {0.};
  int Ngood_gamma = 0;
  
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
  // photon ID is not going to be changed every time this code runs
  // barrel/endcap, Loose/Medium/Tight
  //int    photonID_IsConv[2][3]                = { {0, 0, 0} , {0, 0, 0} };
  double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} , {0.05, 0.05, 0.05} };
  double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
  double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} , {2.3, 1.2, 0.5} };
  double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} , {2.9, 1.5, 1.5} };
  double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} , {0.04, 0.04, 0.04} };
  double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} , {999, 1.0, 1.0} };
  double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };
  // we use loose selection
  int photon_ID = 0;

  for(int ip = 0; ip < fReader->nPho; ++ip)	 //Loop over the photons in a event
    {
      Ngamma[0] += EvtWeight;

      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( fReader->phoEt->at(ip),
                          fReader->phoEta->at(ip),
                          fReader->phoPhi->at(ip),
                          fReader->phoE->at(ip) );

      int region = 0; //barrel
      if( fabs( tmpp4.Eta() )>1.5) region = 1; //endcap

      hphotons["cut0_et"]->Fill( tmpp4.Et(), EvtWeight );
      hphotons["cut0_eta"]->Fill( tmpp4.Eta(), EvtWeight );
      hphotons["cut0_phi"]->Fill( tmpp4.Phi(), EvtWeight );

      if (fVerbose) cout << "photon: et= "<< fReader->phoEt->at(ip) << " eta= "<< fReader->phoEta->at(ip) << endl;

      if( tmpp4.Et() > 25 &&
          fabs( tmpp4.Eta() ) < 1.4442 )
      {
        Ngamma[1]+= EvtWeight;

        bool passelectronveto = (fReader->phoEleVeto->at(ip) == 0);        
        hphotons["cut1_eVeto"]->Fill( int(passelectronveto), EvtWeight );
        if (fVerbose) cout << "photon pass e veto: " << passelectronveto << endl;

        if ( passelectronveto )
          {
            Ngamma[2]+= EvtWeight;
            
            // check conversion tracks
            hphotons["cut1_hasConvTrk"]->Fill( int(fReader->phoIsConv->at(ip)) , EvtWeight );

            float HoverE = fReader->phoHoverE->at(ip);
            hphotons["cut2_HoverE"]->Fill( HoverE, EvtWeight );
            if (fVerbose) cout << "photon H/E= " << HoverE << endl;

            if ( HoverE < photonID_HoverE[region][photon_ID] )
              {
                Ngamma[3]+= EvtWeight;

                float ietaieta = fReader->phoSigmaIEtaIEta->at(ip);
                hphotons["cut3_sigmaietaieta"]->Fill( ietaieta, EvtWeight );
                if (fVerbose) cout << "photon: ietaieta = " << ietaieta << endl;

                float chHadIso = fReader->phoPFChIso->at(ip) - fReader->rho2012 * phoEffArea03ChHad( tmpp4.Eta() );
                float SCFRChIso = fReader->phoSCRChIso->at(ip) - fReader->rho2012 * phoEffArea03ChHad( tmpp4.Eta() );
                float ntHadIso = fReader->phoPFNeuIso->at(ip) - fReader->rho2012 * phoEffArea03NeuHad( tmpp4.Eta() );
                float rand_chHadIso = fReader->phoRandConeChIso->at(ip) - fReader->rho2012 * phoEffArea03ChHad( tmpp4.Eta() );
                float pfIso = fReader->phoPFPhoIso->at(ip)  - fReader->rho2012 * phoEffArea03Pho( tmpp4.Eta() );

                // templates
                if ( ntHadIso < (photonID_RhoCorrR03NeuHadIso_0[region][photon_ID] + tmpp4.Et() * photonID_RhoCorrR03NeuHadIso_1[region][photon_ID]) )
                  {
                    hphotons["temp_sigmaietaieta"]->Fill( ietaieta, EvtWeight );
                    hphotons["temp_SCFRChIso"]->Fill( SCFRChIso, EvtWeight );
                    //hphotons["temp2D_sigmaVsSCFRCh"]->Fill( (double)ietaieta, (double)SCFRChIso, EvtWeight);
                  }

                if ( ietaieta < photonID_SigmaIEtaIEta[region][photon_ID] )
                  {
                    Ngamma[4]+= EvtWeight;

                    hphotons["cut4_chHadIso"]->Fill( chHadIso, EvtWeight );


                    //if ( chHadIso < (photonID_RhoCorrR03ChHadIso[region][photon_ID] + tmpp4.Et() * photonID_RhoCorrR03ChHadIso[region][photon_ID]) )
                    //{
                    //  Ngamma[5]+=EvtWeight;
                    hphotons["cut5_ntHadIso"]->Fill( ntHadIso, EvtWeight );

                    if ( ntHadIso < (photonID_RhoCorrR03NeuHadIso_0[region][photon_ID] + tmpp4.Et() * photonID_RhoCorrR03NeuHadIso_1[region][photon_ID]) )
                      {
                        Ngamma[6]+= EvtWeight;
                        hphotons["cut6_pfIso"]->Fill( chHadIso, EvtWeight );
                        if ( pfIso < ( photonID_RhoCorrR03PhoIso_0[region][photon_ID] + tmpp4.Et() * photonID_RhoCorrR03PhoIso_1[region][photon_ID] ) ) 
                          {
                            Ngamma[7]+= EvtWeight;
                            if ( Ngood_gamma == 0 )
                              {
                                p4photon = tmpp4;
                                hphotons["sigmaietaieta"]->Fill( ietaieta, EvtWeight );
                                hphotons["chHadIso"]->Fill( chHadIso, EvtWeight );
                                hphotons["SCFRChIso"]->Fill( SCFRChIso, EvtWeight );
                                hphotons["RandIso"]->Fill( rand_chHadIso, EvtWeight );
                              }
                            Ngood_gamma++;
                          }
                      }
                      
                  }
              } // HoverE
          } //eVeto
      } //fidutial

    } // end photons

  for ( int ibin=0; ibin < 8; ibin++)
    {
      hphotons["N"]->SetBinContent( ibin, Ngamma[ibin] );
      if ( Ngamma[ibin] > 0 )
        {
          if (ibin==1) cutmap["PhotonFiducial"] += EvtWeight;
          if (ibin==2) cutmap["PhotonEleVeto"] += EvtWeight;
          if (ibin==3) cutmap["PhotonHoverE"] += EvtWeight;
          if (ibin==4) cutmap["PhotonSigmaieta"] += EvtWeight;
          if (ibin==5) cutmap["PhotonchHadIso"] += EvtWeight;
          if (ibin==6) cutmap["PhotonntHadIso"] += EvtWeight;
          if (ibin==7) cutmap["PhotonIso"] += EvtWeight;
          
        }
    }
  hphotons["Ngood"]->Fill( Ngood_gamma );

  if ( Ngood_gamma >0 )
    {

      cutmap["OnePhoton"] += EvtWeight;

      if (fChannel==1)
        {
          hmuons["pt"]->Fill( p4lepton.Pt(), EvtWeight );
          hmuons["eta"]->Fill( p4lepton.Eta(), EvtWeight );
          hmuons["phi"]->Fill( p4lepton.Phi(), EvtWeight );
          hmuons["relisocorr"]->Fill( relIsocorr_Mu, EvtWeight );
        }
      else
        {
          helectrons["pt"]->Fill( p4lepton.Pt(), EvtWeight );
          helectrons["eta"]->Fill( p4lepton.Eta(), EvtWeight );
          helectrons["phi"]->Fill( p4lepton.Phi(), EvtWeight );
          helectrons["relisocorr"]->Fill( relIsocorr_Ele, EvtWeight );
        }
      hPVs["N"]->Fill( fReader->nVtx , EvtWeight );
      hjets["N"]->Fill( Ngood_Jets, EvtWeight );
      hjets["pt1"]->Fill( p4jets[0].Pt(), EvtWeight );
      hMET["MET"]->Fill( fReader->pfMET, EvtWeight );
    }

  // clear vectors
  p4jets.clear();
  vec_btags.clear();

  if (fVerbose) cout << "analysis done." << endl;
   return kTRUE;
}

void ttgamma3::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  Info("SlaveTerminate","begin");

  // fill cutflow histogram
  cout << "Cutflow table" << endl;
  int ibin = 1;
  for ( vector<string>::const_iterator ivec= fCutLabels.begin(); ivec != fCutLabels.end(); ++ivec )
    //  for ( map<string, int >::const_iterator imap=cutmap.begin(); imap!=cutmap.end(); ++imap )
    {
      hcutflow->SetBinContent( ibin, cutmap[ *ivec ] );
      cout << ibin << " cut: " << *ivec << " entries: " << cutmap[ *ivec ] << endl;
      ibin++;
    }

  if ( fN_tt_filter > 0 )
    cout << "Number of ttjets events filtered = " << fN_tt_filter << endl;

  // Write the ntuple to the file
  if (fFile) {
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    if (fVerbose) cout << "number of entries in h1test histogram = " << h1test->GetEntries() << endl;

    if(h1test->GetEntries() > 0){
      fFile->cd();
      h1test->Write();
      h1test->SetDirectory(0);
      hcutflow->Write();
      hcutflow->SetDirectory(0);

      if (! fdoSkim )
        {
          WriteHistograms("electrons", helectrons);
          WriteHistograms("muons", hmuons);
          WriteHistograms("jets", hjets);
          WriteHistograms("btag", hbtag);
          WriteHistograms("photons", hphotons);
          WriteHistograms("MET", hMET);
          WriteHistograms("PV", hPVs);
          WriteHistograms("MC", hMC);
          WriteHistograms("mass",hM);
        }
      else
        {
          fReader->AutoSave();
        }

      fFile->cd();
      
      //// SKIM
      //if (fdoSkim)
      //  {
          //fFile->mkdir("ggNtuplizer");
          //fFile->cd("ggNtuplizer");
          //TTree *ChainSkim = fReader->GetSkim();
          //ChainSkim->Print();
          //ChainSkim->AutoSave();
      //}
      //fProofFile->Print();
      //fOutput->Add(fProofFile);

    } else {
      cleanup = kTRUE;
    }
        
    Info("SlaveTerminate", "Everything saved. Closing file.");

    gDirectory = savedir;
    fFile->Close();
    // Cleanup, if needed
    if (cleanup) {
      Info("SlaveTerminate", "nothing to save: just cleanup everything ...");
      TUrl uf(*(fFile->GetEndpointUrl()));
      SafeDelete(fFile);
      gSystem->Unlink(uf.GetFile());
      SafeDelete(fProofFile);
    } else {

      //Info("SlaveTerminate", "objects saved into '%s%s': sending related TProofOutputFile ...",
      //    fProofFile->GetFileName(), fProofFile->GetOptionsAnchor());

    }

  }
  
  Info("SlaveTerminate", "done.");
}

void ttgamma3::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  Info("Terminate","Analyzer done.");
}
