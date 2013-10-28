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
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "../../EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

void ttgamma3::ParseInput()
{

  if (fMyOpt.Contains("muon"))
    {
      fChannel = 1;
    }
  if (fMyOpt.Contains("electron"))
    {
      fChannel = 2;
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
          Info("Begin","This sample is treated as MC");
        }
    }

  if (fIsMC) fPUreweighting = true;
  else fPUreweighting = false;

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
   ParseInput();
   //fVerbose = true;
   Info("Begin", "starting with process option: %s", option.Data());
   if (fPUreweighting) Info("Begin", "Apply PU reweighting");
}

void ttgamma3::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  fMyOpt = option;
  ParseInput();

  //initialize the Tree branch addresses
  Init(tree);

  // get PU weights
  if (fPUreweighting) {
    TFile *PU_data_file     = new TFile( "MyDataPileupHistogram.root" );
    hPU_weights = new TH1F( *(static_cast<TH1F*>(PU_data_file->Get( "pileup" )->Clone() )) );
    TFile *PU_MC_file = new TFile( "MyMCPileupHistogram.root");
    TString hmc_name = "hPUTrue_";
    hmc_name += fSample;
    TH1F* hPU_mc = new TH1F( *(static_cast<TH1F*>(PU_MC_file->Get( hmc_name )->Clone() )) );
    hPU_weights->Scale(1./(hPU_weights->Integral()));
    hPU_mc->Scale(1./(hPU_mc->Integral()));
    hPU_weights->Divide( hPU_mc );

  }

  // Create Output file
  TString tmpfilename = "results";
  if ( fSample != "" ) tmpfilename = fOutdir+"results_"+fSample+".root";
  else tmpfilename = fOutdir + "results_data.root";
  fFile = new TFile( tmpfilename, "RECREATE");
  //Info("SlaveBegin", "Output filename: %s", tmpfilename );

  TString hname = "_"+fSample;
  // vertices
  hPVs["N"] = new TH1F("NPV"+hname,"Number of PVs",35,-0.5,34.5);
  // muons
  hmuons["pt"] = new TH1F("muon_pt"+hname,"p_{T}^{#mu} [GeV/c]", 50, 0, 500);
  hmuons["eta"] = new TH1F("muon_eta"+hname,"#eta^{#mu}", 20, -2.1, 2.1);
  hmuons["phi"] = new TH1F("muon_phi"+hname,"#phi^{#mu}", 30, -3.15, 3.15);
  // electrons
  helectrons["pt"] = new TH1F("electron_pt"+hname,"p_{T}^{e} [GeV/c]", 50, 20, 500);
  helectrons["eta"] = new TH1F("electron_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  helectrons["phi"] = new TH1F("electron_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);
  helectrons["reliso"] = new TH1F("electron_reliso"+hname,"Relative Isolation", 40, 0, 0.2);
  // jets
  hjets["pt"] = new TH1F("jet_pt"+hname,"jet p_{T} [GeV/c]", 60, 30, 800);
  // MET
  hMET["MET"] = new TH1F("MET"+hname,"Missing Transverse Energy [GeV]", 50, 0, 300);
  // photons
  hphotons["pt"] = new TH1F("photon_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);

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
       fCutLabels.push_back("Skim");
       fCutLabels.push_back("HLT");
       fCutLabels.push_back("GoodPV");
       fCutLabels.push_back("OneIsoMuon");
       fCutLabels.push_back("LooseMuVeto");
       fCutLabels.push_back("LooseEleVeto");
       fCutLabels.push_back("3Jets");
       fCutLabels.push_back("4Jets");
       fCutLabels.push_back("1bJet");
       fCutLabels.push_back("1Photon");
     }
   else
     { //electron+jets
       fCutLabels.push_back("Skim");
       fCutLabels.push_back("HLT");
       fCutLabels.push_back("GoodPV");
       fCutLabels.push_back("OneIsoEle");
       fCutLabels.push_back("LooseMuVeto");
       fCutLabels.push_back("LooseEleVeto");
       fCutLabels.push_back("3Jets");
       fCutLabels.push_back("4Jets");
       fCutLabels.push_back("1bJet");
       fCutLabels.push_back("1Photon");
     }

   hcutflow = new TH1F("cutflow","cut flow", fCutLabels.size(), 0.5, fCutLabels.size() +0.5 );

   for ( vector<string>::const_iterator ivec= fCutLabels.begin(); ivec!=fCutLabels.end(); ++ivec)
     {
       cutmap[ *ivec ] = 0;
     }

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

  if (fVerbose) cout << "Processing entry: " << entry << endl;

  fChain->GetEntry(entry);

  ////////////////////////////////////////
  // get PU weight for MC samples
  ///////////////////////////////////////
  if (fPUreweighting) {
    // for BX 0 only
    float numInter = puTrue[1];
    EvtWeight = EvtWeight*hPU_weights->GetBinContent( hPU_weights->FindBin( numInter )  );
    if (fVerbose) cout << " PU weight = " << EvtWeight << endl;
  }

  cutmap["Skim"] += EvtWeight;

  /////////////////////////////////////////
  //  HLT
  /////////////////////////////////////////
  if ( HLTIndex[17] == 0 )
    { return kTRUE; } // 17 --> HLT_Ele27_WP80_v Trigger
  
  cutmap["HLT"] += EvtWeight;
  if (fVerbose) cout << "Pass HLT" << endl;
  ////////////////////////////////////
  // PRIMARY VERTICES
  ///////////////////////////////////
  // check if first vertex is a good vtx
  if ( IsVtxGood != 0 ) return kTRUE;
  cutmap["GoodPV"] += EvtWeight;
  hPVs["N"]->Fill( nVtx , EvtWeight );

  if (fVerbose) cout << "Pass good vertex" << endl;
  //////////////////////////////////
  // ELECTRONS
  //////////////////////////////////
  if (fVerbose) cout << "Total number of electrons nEle = "<< nEle << endl;

  for(int ie = 0; ie < nEle; ++ie)		//Loop over the electrons in a event
    {
      p4ele.SetPtEtaPhiE( elePt[ie],
                          eleSCEta[ie],
                          elePhi[ie],
                          eleEcalEn[ie] );

      h1test->Fill( p4ele.Pt() );

      float AEff03 = 0.00;

      if (!fIsMC) {
        AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eleSCEta[ie], ElectronEffectiveArea::kEleEAData2011);
      } 
      else {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eleSCEta[ie], ElectronEffectiveArea::kEleEAFall11MC);
      }

      double relIso = (elePFChIso03[ie] + elePFNeuIso03[ie] + elePFPhoIso03[ie])/ p4ele.Pt();
      double relIsorho = ( elePFChIso03[ie] + fmax(0.0, elePFNeuIso03[ie] + elePFPhoIso03[ie] - rho25_elePFiso*AEff03) )/ p4ele.Pt();

      if ( p4ele.Pt() > 35.0 &&
           fabs(p4ele.Eta()) < 2.5 &&
           !( (1.4442 < fabs( eleSCEta[ie] )) && fabs( eleSCEta[ie]) < 1.5660) &&
           eleD0[ie] < 0.02 &&
           eleIDMVATrig[ie] > 0.5 &&
           eleMissHits[ie] <= 0 )
        {


          helectrons["pt"]->Fill( p4ele.Pt() );
        }
    } // electrons
    
  /*  
    for(int ip = 0; ip < nPho; ++ip)		//Loop over the photons in a event
    {  
      if(nPho >= 1 && phoEt[ip] > 25)
      {
	double highet = phoEt[0];
    
	double gammaphi = phoPhi[0];
    
	double gammaeta = phoEta[0];

      }
    }
  
    if(nJet >= 3)		//Cut for Jets
    {
      double leadjet = jetPt[0];
    
      double leadjetphi = jetPhi[0];
    
      double leadjeteta = jetEta[0];
    
    }
  */
  
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
  int ibin = 1;
  for ( vector<string>::const_iterator ivec= fCutLabels.begin(); ivec != fCutLabels.end(); ++ivec )
    //  for ( map<string, int >::const_iterator imap=cutmap.begin(); imap!=cutmap.end(); ++imap )
    {
      hcutflow->SetBinContent( ibin, cutmap[ *ivec ] );
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

      WriteHistograms("electrons", helectrons);
      WriteHistograms("muons", hmuons);
      WriteHistograms("jets", hjets);
      WriteHistograms("photons", hphotons);
      WriteHistograms("MET", hMET);
      WriteHistograms("PV", hPVs);

      fFile->cd();
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
