#define ttgamma_cxx
// The class definition in ttgamma.h has been generated automatically
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
// Root > T->Process("ttgamma.C")
// Root > T->Process("ttgamma.C","some options")
// Root > T->Process("ttgamma.C+")
//

#include "ttgamma.h"
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

void ttgamma::ParseInput()
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
    fOutdir = tmp+"/";
    Info("Begin","output files will be written to directory: %s", fOutdir.Data());
  }
  if (fMyOpt.Contains("sample"))
    {
      TString tmp = fMyOpt;
      tmp = tmp.Remove(0,fMyOpt.Index("sample")+7);
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
}

void ttgamma::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void ttgamma::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   h1test = new TH1F("h1test","Electron p_{T}",100,10.,400);
   h2test = new TH1F("h2test","Leading Electron P_{T}",100,10.,400);
   h3test = new TH1F("h3test","Leading Electron eta",100,-4.,4);
   h4test = new TH1F("h4test","Leading Electron phi",100,-4.,4);
   h5test = new TH1F("h5test","Number of Electrons",10,0.,10);   
   
   h6test = new TH1F("h6test","Leading Photon E_{T}",100,10.,400);
   h7test = new TH1F("h7test","Leading Photon eta",100,-4.,4);
   h8test = new TH1F("h8test","Leading Photon phi",100,-4.,4);
   h9test = new TH1F("h9test","Number of Photons",10,0.,10);      
   
   h10test = new TH1F("h10test","Leading Jet P_{T}",100,10.,400);
   h11test = new TH1F("h11test","Leading Jet eta",100,-4.,4);
   h12test = new TH1F("h12test","Leading Jet phi",100,-4.,4);
   h13test = new TH1F("h13test","Number of Jets",70,0.,70);  
   
   h14test = new TH1F("h14test","Number of Primary Vertex",50,0.,50); 
   
   h1test->SetDirectory(fFile);   
   h2test->SetDirectory(fFile);
   h3test->SetDirectory(fFile);   
   h4test->SetDirectory(fFile);   
   h5test->SetDirectory(fFile);   
   
   h6test->SetDirectory(fFile);
   h7test->SetDirectory(fFile);   
   h8test->SetDirectory(fFile);   
   h9test->SetDirectory(fFile);     
   
   h10test->SetDirectory(fFile);
   h11test->SetDirectory(fFile);   
   h12test->SetDirectory(fFile);   
   h13test->SetDirectory(fFile);  
   
   h14test->SetDirectory(fFile);  

   //initialize the Tree branch addresses
   Init(tree);

   // We may be creating a dataset or a merge file: check it
   TNamed *nm = dynamic_cast<TNamed *>(fInput->FindObject("SimpleNtuple.root"));
   if (nm) {
     // Just create the object
     
     UInt_t opt = TProofOutputFile::kRegister | TProofOutputFile::kOverwrite | TProofOutputFile::kVerify;
     fProofFile = new TProofOutputFile( "SimpleNtuple.root",TProofOutputFile::kDataset, opt, nm->GetTitle());
   } else {
     // For the ntuple, we use the automatic file merging facility
     // Check if an output URL has been given
     TNamed *out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE_LOCATION");
     Info("SlaveBegin", "PROOF_OUTPUTFILE_LOCATION: %s", (out ? out->GetTitle() : "undef"));
     TString tmpfilename = "results";
     if ( fSample != "" ) tmpfilename += "_"+fSample+".root";
     else tmpfilename = "results_ttjets1.root";
     fProofFile = new TProofOutputFile(tmpfilename, (out ? out->GetTitle() : "M"));
     out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE");
     if (out) fProofFile->SetOutputFileName(fOutdir + out->GetTitle());
   }

   if (!(fFile = fProofFile->OpenFile("RECREATE"))) {
     Warning("SlaveBegin", "problems opening file: %s/%s",
             fProofFile->GetDir(), fProofFile->GetFileName());
   }

}

Bool_t ttgamma::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ttgamma::GetEntry() or TBranch::GetEntry()
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

  fChain->GetEntry(entry);
  
  if(HLTIndex[17] == 0){return kTRUE;} // 17 --> HLT_Ele27_WP80_v Trigger
  //{
        
    h14test->Fill( nVtx ); 
    
    for(int ie = 0; ie < nEle; ++ie)		//Loop over the electrons in a event
    {
      if(nEle >= 1 && elePt[ie] > 35)
      {
	double pt = elePt[ie];
	h1test->Fill( pt );      
      
	double leadpt = elePt[0];
	h2test->Fill( leadpt ); 
    
	double leadphi = elePhi[0];
	h3test->Fill( leadphi );
    
	double leadeta = eleSCEta[0];
	h4test->Fill( leadeta );
    
	h5test->Fill( nEle );  
      }
    }
  
    for(int ip = 0; ip < nPho; ++ip)		//Loop over the photons in a event
    {  
      if(nPho >= 1 && phoEt[ip] > 25)
      {
	double highet = phoEt[0];
	h6test->Fill( highet ); 
    
	double gammaphi = phoPhi[0];
	h7test->Fill( gammaphi );
    
	double gammaeta = phoEta[0];
	h8test->Fill( gammaeta );
    
	h9test->Fill( nPho );  
      }
    }
  
    if(nJet > 3)		//Cut for Jets
    {
      double leadjet = jetPt[0];
      h10test->Fill( leadjet ); 
    
      double leadjetphi = jetPhi[0];
      h11test->Fill( leadjetphi );
    
      double leadjeteta = jetEta[0];
      h12test->Fill( leadjeteta );
    
      h13test->Fill( nJet );  
    }

  //}
  
  if (fVerbose) cout << "analysis done." << endl;
   return kTRUE;
}

void ttgamma::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  // Write the ntuple to the file
  if (fFile) {
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    if (fVerbose) cout << "number of entries in h1test histogram = " << h1test->GetEntries() << endl;

    if(h1test->GetEntries() > 0){
      fFile->cd();
      h1test->Write();
      h1test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }
    
    if(h2test->GetEntries() > 0){
      fFile->cd();
      h2test->Write();
      h2test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }
    
    if(h3test->GetEntries() > 0){
      fFile->cd();
      h3test->Write();
      h3test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }    

    if(h4test->GetEntries() > 0){
      fFile->cd();
      h4test->Write();
      h4test->SetDirectory(0);
    } else {
      cleanup = kTRUE;    } 
    
    if(h5test->GetEntries() > 0){
      fFile->cd();
      h5test->Write();
      h5test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }    

    if(h6test->GetEntries() > 0){
      fFile->cd();
      h6test->Write();
      h6test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }
    
    if(h7test->GetEntries() > 0){
      fFile->cd();
      h7test->Write();
      h7test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }    

    if(h8test->GetEntries() > 0){
      fFile->cd();
      h8test->Write();
      h8test->SetDirectory(0);
    } else {
      cleanup = kTRUE;    } 
    
    if(h9test->GetEntries() > 0){
      fFile->cd();
      h9test->Write();
      h9test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }  
    
    if(h10test->GetEntries() > 0){
      fFile->cd();
      h10test->Write();
      h10test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }
    
    if(h11test->GetEntries() > 0){
      fFile->cd();
      h11test->Write();
      h11test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }    

    if(h12test->GetEntries() > 0){
      fFile->cd();
      h12test->Write();
      h12test->SetDirectory(0);
    } else {
      cleanup = kTRUE;    } 
    
    if(h13test->GetEntries() > 0){
      fFile->cd();
      h13test->Write();
      h13test->SetDirectory(0);
    } else {
      cleanup = kTRUE;
    }    
    
    if(h14test->GetEntries() > 0){
      fFile->cd();
      h14test->Write();
      h14test->SetDirectory(0);
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

      Info("SlaveTerminate", "objects saved into '%s%s': sending related TProofOutputFile ...",
           fProofFile->GetFileName(), fProofFile->GetOptionsAnchor());
      fProofFile->Print();
      fOutput->Add(fProofFile);
    }
  }

}

void ttgamma::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  Info("Terminate","Analyzer done.");
}
