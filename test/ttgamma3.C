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
          fIsMC = true;
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
   if (fPUreweighting) Info("Begin", "PU reweighting is ON");
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
  // MET
  hMET["MET"] = new TH1F("MET"+hname,"Missing Transverse Energy [GeV]", 50, 0, 300);
  // photons
  hphotons["pt"] = new TH1F("photon_pt"+hname,"p_{T}^{#gamma} [GeV/c]", 50, 20, 500);
  hphotons["eta"] = new TH1F("photon_eta"+hname,"#eta^{e}", 20, -2.1, 2.1);
  hphotons["phi"] = new TH1F("photon_phi"+hname,"#phi^{#mu}", 20, 0, 3.15);

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
       fCutLabels.push_back("DileptonVeto");
       fCutLabels.push_back("1Jets");
       fCutLabels.push_back("2Jets");
       fCutLabels.push_back("3Jets");
       fCutLabels.push_back("4Jets");
       fCutLabels.push_back("Onebtag");
       fCutLabels.push_back("OnePhoton");
     }
   else
     { //electron+jets
       fCutLabels.push_back("Skim");
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
  // required NDF > 4 && fabs(PVz)<= 24 && rho <= 2 
  if ( IsVtxGood != 0 ) return kTRUE;
  cutmap["GoodPV"] += EvtWeight;
  hPVs["N"]->Fill( nVtx , EvtWeight );

  if (fVerbose) cout << "Pass good vertex" << endl;

  //////////////////////////////////
  // ELECTRONS
  //////////////////////////////////
  if (fVerbose) cout << "Total number of electrons nEle = "<< nEle << endl;

  int Ngood_Ele = 0;
  int Nloose_Ele = 0;
  float relIso_Ele = -1;
  float relIsocorr_Ele = -1;
  bool passconversionveto = false;

  for(int ie = 0; ie < nEle; ++ie)		//Loop over the electrons in a event
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( elePt[ie],
                          eleSCEta[ie],
                          elePhi[ie],
                          eleEcalEn[ie] );

      h1test->Fill( tmpp4.Pt() );

      float AEff03 = 0.00;

      if (!fIsMC) {
        AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eleSCEta[ie], ElectronEffectiveArea::kEleEAData2011);
      } 
      else {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eleSCEta[ie], ElectronEffectiveArea::kEleEAFall11MC);
      }

      float relIso = (elePFChIso03[ie] + elePFNeuIso03[ie] + elePFPhoIso03[ie])/ tmpp4.Pt();
      float relIsocorr = ( elePFChIso03[ie] + fmax(0.0, elePFNeuIso03[ie] + elePFPhoIso03[ie] - rho25_elePFiso*AEff03) )/ tmpp4.Pt();

      // loose electrons
      if ( tmpp4.Pt() > 20.0 &&
           fabs(tmpp4.Eta()) < 2.5 &&
           eleIDMVATrig[ie] > 0.0 && eleIDMVATrig[ie] < 1.0 &&
           relIsocorr < 0.15 )
        {
          Nloose_Ele++;

          // good tight electrons
          if ( tmpp4.Pt() > 30.0 &&
               fabs(tmpp4.Eta()) < 2.5 &&
               !( (1.4442 < fabs( eleSCEta[ie] )) && fabs( eleSCEta[ie]) < 1.5660) &&
               eleD0[ie] < 0.02 &&
               eleIDMVATrig[ie] > 0.5 && eleIDMVATrig[ie] < 1.0 &&
               eleMissHits[ie] <= 0 &&
               relIsocorr < 0.1 )
            {
              // leading electron
              if ( Ngood_Ele == 0 )
                {
                  p4ele.SetPtEtaPhiE( tmpp4.Pt(), tmpp4.Eta(), tmpp4.Phi(), tmpp4.E() );
                  relIso_Ele = relIso;
                  relIsocorr_Ele = relIsocorr;
                  passconversionveto = (eleConvVtxFit[ie] == 0);
                }
              Ngood_Ele++;
          //helectrons["pt"]->Fill( tmpp4.Pt() );
            }
        } //loose

    } // end electrons
    
  //////////////////////////
  // MUONS
  /////////////////////////
  if (fVerbose) cout << "Total number of muons nMu = "<< nMu << endl;

  int Ngood_Mu = 0;
  int Nloose_Mu = 0;
  float relIso_Mu = -1;
  float relIsocorr_Mu = -1;

  for(int imu = 0; imu < nMu; ++imu)             //Loop over the muons in a event
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( muPt[imu],
                          muEta[imu],
                          muPhi[imu],
                          sqrt(muPt[imu]*muPt[imu] + muPz[imu]*muPz[imu]) );

      float relIso = (muPFIsoR04_CH[imu] + muPFIsoR04_NH[imu])/ tmpp4.Pt();
      // delta beta corrections
      //  I = [sumChargedHadronPt+ max(0.,sumNeutralHadronPt+sumPhotonPt-0.5sumPUPt]/pt
      float relIsocorr = ( muPFIsoR04_CH[imu] + fmax(0.0, muPFIsoR04_NH[imu] + muPFIsoR04_Pho[imu] -0.5*muPFIsoR04_PU[imu]) )/ tmpp4.Pt();
      // loose muon selection
      if ( tmpp4.Pt() > 10.0 &&
           fabs(tmpp4.Eta()) < 2.5 &&
           relIsocorr < 0.2 )
        {
          Nloose_Mu++;
          // good muons
          if ( tmpp4.Pt() > 26.0 &&
               fabs(tmpp4.Eta()) < 2.1 &&
               muChi2NDF[imu] < 10 &&
               muNumberOfValidTrkLayers[imu] > 5 &&
               muNumberOfValidMuonHits[imu] > 0 &&
               muD0[imu] < 0.2 &&
               fabs( muDz[imu] ) < 0.5 && //check this
               muNumberOfValidPixelHits[imu] > 0 &&
               muStations[imu] > 1 &&
               relIsocorr < 0.12 )
            {
              // leading muon
              if ( Ngood_Mu == 0 )
                {
                  p4muon.SetPtEtaPhiE( tmpp4.Pt(), tmpp4.Eta(), tmpp4.Phi(), tmpp4.E() );
                  relIso_Mu = relIso;
                  relIsocorr_Mu = relIsocorr;
                }
              
              Ngood_Mu++;
            }
        }
    } // end muons

  ////////////////////////////////
  // JETS
  ///////////////////////////////
  if (fVerbose) cout << "Total number of jets nJet = "<< nJet << endl;

  int Ngood_Jets = 0;
  float met_x = 0.0;
  float met_y = 0.0;

  for (int ij= 0; ij < nJet; ++ij)
    {
      TLorentzVector tmpp4;
      tmpp4.SetPtEtaPhiE( jetPt[ij], jetEta[ij], jetPhi[ij], jetEn[ij] );

      float uncorr_jet_pt = jetRawPt[ij];
      met_x += uncorr_jet_pt*TMath::Cos( jetPhi[ij] );
      met_y += uncorr_jet_pt*TMath::Sin( jetPhi[ij] );

      // Apply JER smearing in MC
      if ( fIsMC && jetPt[ij] > 10 )
        {
          float gen_pt = jetGenJetPt[ij];
          if ( gen_pt < 0 ) continue;
          float jet_pt = jetPt[ij];
          float jet_eta = jetEta[ij];
          // factor is (c - 1), where c are the eta-dependent scale factors, to be taken from the official twiki
          float c = 0;
          // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          if ( fabs(jet_eta) > 0 && fabs(jet_eta) < 0.5 )  c = 1.052;
          if ( fabs(jet_eta) > 0.5 && fabs(jet_eta) < 1.1) c = 1.057;
          if ( fabs(jet_eta) > 1.1 && fabs(jet_eta) < 1.7) c = 1.096;
          if ( fabs(jet_eta) > 1.7 && fabs(jet_eta) < 2.3) c = 1.134;
          if ( fabs(jet_eta) > 2.3 && fabs(jet_eta) < 5.0) c = 1.288;
          float factor = (c-1);
          float deltapt = (jet_pt - gen_pt) * factor;
          float ptscale = fmax(0.0, (jet_pt + deltapt) / jet_pt);
          tmpp4 *= ptscale;
          met_x -= uncorr_jet_pt*TMath::Cos( jetPhi[ij] ) * ptscale;
          met_y -= uncorr_jet_pt*TMath::Sin( jetPhi[ij] ) * ptscale;
        }

      if ( tmpp4.Pt() > 30 &&
           fabs( tmpp4.Eta() ) < 2.5 &&
           jetPFLooseId[ij] == true )
        {
          if (Ngood_Jets == 0 && tmpp4.Pt() > 55 )
            {
              p4jets.push_back( tmpp4 );
            }
          if (Ngood_Jets == 1 && tmpp4.Pt() > 45 )
            {
              p4jets.push_back( tmpp4 );
            }
          if (Ngood_Jets == 2 && tmpp4.Pt() > 35 )
            {
              p4jets.push_back( tmpp4 );
            }
          if (Ngood_Jets == 3 && tmpp4.Pt() > 20 )
            {
              p4jets.push_back( tmpp4 );
            }
          Ngood_Jets++;
        }

    } // end jets

  if ( p4jets.size() >= 1 )
    cutmap["1Jets"] += EvtWeight;
  if ( p4jets.size() >= 2 )
    cutmap["2Jets"] += EvtWeight;
  if ( p4jets.size() >=3 )
    cutmap["3Jets"] += EvtWeight;
  if ( p4jets.size() >=4 )
    cutmap["4Jets"] += EvtWeight;

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
