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

  if (fIsMC) fPUreweighting = true;
  else fPUreweighting = false;

  if (fMyOpt.Contains("sync"))
    {
      fdoSync = true;
      fdoJER = false;
      fdoHLT = false;
      fPUreweighting = false;
    }
  
  if (fMyOpt.Contains("skim"))
    {
      fdoSkim = true;
      Info("Begin","Running in SKIM mode");
      fPUreweighting = false;
      fdoJER = false;
    }
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
    TFile *PU_data_file     = new TFile( "HelperFiles/MyDataPileupHistogram.root" );
    if ( PU_data_file->IsZombie() ) {
      Info("SlaveBegin","Error opening file: %s", PU_data_file->GetName() );
    }
    hPU_weights = new TH1F( *(static_cast<TH1F*>(PU_data_file->Get( "pileup" )->Clone() )) );
    TFile *PU_MC_file = new TFile( "HelperFiles/MyMCPileupHistogram.root");
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
  fFile = new TFile( tmpfilename, "RECREATE");
  //Info("SlaveBegin", "Output filename: %s", tmpfilename );

  TString hname = "_"+fSample;
  // vertices
  hPVs["N"] = new TH1F("NPV"+hname,"Number of PVs",35,-0.5,34.5);
  hPVs["N_fin"] = new TH1F("NPV_fin"+hname,"Number of PVs",35,-0.5,34.5);
  // muons
  hmuons["pt"] = new TH1F("muon_pt"+hname,"p_{T}^{#mu} [GeV/c]", 50, 0, 500);
  hmuons["eta"] = new TH1F("muon_eta"+hname,"#eta^{#mu}", 20, -2.1, 2.1);
  hmuons["phi"] = new TH1F("muon_phi"+hname,"#phi^{#mu}", 30, -3.15, 3.15);
  hmuons["reliso"] = new TH1F("muon_reliso"+hname,"I_{rel}", 40, 0, 0.2);
  hmuons["relisocorr"] = new TH1F("muon_relisocorr"+hname,"I_{rel}^{corr}", 40, 0, 0.2);
  // electrons
  helectrons["pt"] = new TH1F("electron_pt"+hname,"p_{T}^{e} [GeV/c]", 50, 20, 500);
  helectrons["eta"] = new TH1F("electron_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  helectrons["phi"] = new TH1F("electron_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  helectrons["reliso"] = new TH1F("electron_reliso"+hname,"I_{rel}", 40, 0, 0.2);
  helectrons["relisocorr"] = new TH1F("electron_relisocorr"+hname,"I_{rel}^{corr}", 40, 0, 0.2);
  // jets
  hjets["pt"] = new TH1F("jet_pt"+hname,"jet p_{T} [GeV/c]", 60, 30, 800);
  hjets["pt1"] = new TH1F("jet_pt1"+hname,"1st jet p_{T} [GeV/c]", 60, 30, 800);
  hjets["pt2"] = new TH1F("jet_pt2"+hname,"2nd jet p_{T} [GeV/c]", 60, 30, 800);
  hjets["pt3"] = new TH1F("jet_pt3"+hname,"3rd jet p_{T} [GeV/c]", 60, 30, 800);
  hjets["pt4"] = new TH1F("jet_pt4"+hname,"4th jet p_{T} [GeV/c]", 60, 30, 800);
  hjets["deltaR_lepton"] = new TH1F("deltaR_lepton"+hname,"#DeltaR(jet,lepton)", 50, 0, 7);
  // MET
  hMET["MET"] = new TH1F("MET"+hname,"Missing Transverse Energy [GeV]", 50, 0, 300);
  // photons
  hphotons["cut0_pt"] = new TH1F("photon_cut0_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  hphotons["cut0_eta"] = new TH1F("photon_cut0_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  hphotons["cut0_phi"] = new TH1F("photon_cut0_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  hphotons["cut1_pt"] = new TH1F("photon_cut1_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  hphotons["cut1_eta"] = new TH1F("photon_cut1_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  hphotons["cut1_phi"] = new TH1F("photon_cut1_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  hphotons["cut2_pt"] = new TH1F("photon_cut2_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  hphotons["cut2_eta"] = new TH1F("photon_cut2_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  hphotons["cut2_phi"] = new TH1F("photon_cut2_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  hphotons["HoverE"] = new TH1F("photon_HoverE"+hname,"H/E",20,0,1);
  hphotons["chHadIso"] = new TH1F("photon_chHadIso"+hname,"Charged Hadron Isolation",20,0,10);

  map<string,TH1* > allhistos = hmuons;
  allhistos.insert( helectrons.begin(), helectrons.end() );
  allhistos.insert( hphotons.begin(), hphotons.end() );
  allhistos.insert( hmuons.begin(), hmuons.end() );
  allhistos.insert( hjets.begin(), hjets.end() );
  allhistos.insert( hMET.begin(), hMET.end() );

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
       fCutLabels.push_back("Processeed");
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
        cutmap["OneIsoMuon"] += EvtWeight;
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

  for (int ij= 0; ij < fReader->nJet; ++ij)
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( fReader->jetPt->at(ij), fReader->jetEta->at(ij), fReader->jetPhi->at(ij), fReader->jetEn->at(ij) );

      float uncorr_jet_pt = fReader->jetRawPt->at(ij);
      met_x += uncorr_jet_pt*TMath::Cos( fReader->jetPhi->at(ij) );
      met_y += uncorr_jet_pt*TMath::Sin( fReader->jetPhi->at(ij) );

      // Apply JER smearing in MC
      if ( fIsMC && fdoJER && fReader->jetPt->at(ij) > 10 )
        {
          float gen_pt = fReader->jetGenJetPt->at(ij);
          if ( gen_pt < 0 ) continue;
          float jet_pt = fReader->jetPt->at(ij);
          float jet_eta = fReader->jetEta->at(ij);
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

      if ( tmpp4.Pt() > 10 &&
           fabs( tmpp4.Eta() ) < 2.5 &&
           fReader->jetPFLooseId->at(ij) == true )
        {
          
          float tmpDeltaR = tmpp4.DeltaR(p4lepton);
          
          if (fVerbose) cout << "Ngood_Jets = " << Ngood_Jets << " Pt = " << tmpp4.Pt() << " DeltaR(lepton)= " << tmpDeltaR << endl;

          if ( tmpDeltaR < 0.3 ) continue;
          //if ( aDeltaR <= jminDeltaR ) jminDeltaR = aDeltaR;

          if (Ngood_Jets == 0 && tmpp4.Pt() > 55.0 )
            {
              p4jets.push_back( tmpp4 );
              pass1stJet = true;
              if ( fReader->jetCombinedSecondaryVtxBJetTags->at(ij) > CSVM )
                {
                  vec_btags.push_back( true );
                  Nbtags++;
                }
              else
                vec_btags.push_back( false );

              hjets["deltaR_lepton"]->Fill(tmpDeltaR);
            }

          if (Ngood_Jets == 1 && tmpp4.Pt() > 45.0 )
            {
              p4jets.push_back( tmpp4 );
              pass2ndJet = true;
              if ( fReader->jetCombinedSecondaryVtxBJetTags->at(ij) > CSVM )
                {
                  vec_btags.push_back( true );
                  Nbtags++;
                }
              else
                vec_btags.push_back( false );
            }

          if (Ngood_Jets == 2 && tmpp4.Pt() > 35.0 )
            {
              p4jets.push_back( tmpp4 );
              pass3rdJet = true;
              if ( fReader->jetCombinedSecondaryVtxBJetTags->at(ij) > CSVM )
                {
                  vec_btags.push_back( true );
                  Nbtags++;
                }
              else
                vec_btags.push_back( false );
            }

          if (Ngood_Jets == 3 && tmpp4.Pt() > 20.0 )
            {
              p4jets.push_back( tmpp4 );
              pass4thJet = true;
              if ( fReader->jetCombinedSecondaryVtxBJetTags->at(ij) > CSVM )
                {
                  vec_btags.push_back( true );
                  Nbtags++;
                }
              else
                vec_btags.push_back( false );
            }
          
          Ngood_Jets++;
        }

    } // end jets

  if (fVerbose) cout << "Pass jets = " << pass1stJet << " " << pass2ndJet << " " << pass3rdJet << " " << pass4thJet << endl;

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

  //hPVs["N"]->Fill( fReader->nVtx , EvtWeight );
  h1test->Fill( p4lepton.Pt() ); // leading pT lepton

  ////////////////////////
  // SKIM Output
  ////////////////////////

  if (fdoSkim) 
    {
      fReader->FillSkim( fFile );
      if (entry%1000 == 1) fReader->AutoSave();
    }
  ////////////////////////
  // b-tagging selection
  ///////////////////////

  if ( Nbtags >= 1 )
    {
      cutmap["Onebtag"] += EvtWeight;
    }
  else
    return kTRUE;

  ///////////////////////////////////
  // PHOTONS
  ///////////////////////////////////
  int Nloose_gamma = 0;
  int Ngood_gamma = 0;

  for(int ip = 0; ip < fReader->nPho; ++ip)		//Loop over the photons in a event
    {
  
      if( fReader->phoEt->at(ip) > 25 &&
          fabs( fReader->phoEta->at(ip) ) < 1.44 )
      {
        hphotons["cut0_pt"]->Fill( fReader->phoEt->at(ip) );
        hphotons["cut0_eta"]->Fill( fReader->phoEta->at(ip) );
        
        bool passelectronveto = fReader->phoEleVeto->at(ip);
        
        if ( passelectronveto )
          {

          }

      }
    } // end photons
  
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
          WriteHistograms("photons", hphotons);
          WriteHistograms("MET", hMET);
          WriteHistograms("PV", hPVs);
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

}

void ttgamma3::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  Info("Terminate","Analyzer done.");
}
