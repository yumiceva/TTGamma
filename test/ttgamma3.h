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

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxeleESEffSigmaRR = 1;
const Int_t kMaxphoESEffSigmaRR = 1;
const Int_t kMaxnPFEle = 6;
const Int_t kMaxPFElePt = 1;
const Int_t kMaxPFEleEta = 1;
const Int_t kMaxPFElePhi = 1;
const Int_t kMaxPFEleEn = 1;
const Int_t kMaxPFEleCharge = 1;

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
  map< string, TH1*> hPVs;
  map< string, TH1*> hMET;
  void               WriteHistograms(const char* name, map<string, TH1*> hcontainer);
  TFile          *fFile;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[443];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[63][3];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nVtxBS;
   Float_t         vtxbs[63][3];   //[nVtxBS]
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   Int_t           mcPID[12];   //[nMC]
   Float_t         mcVtx[12][3];   //[nMC]
   Float_t         mcPt[12];   //[nMC]
   Float_t         mcMass[12];   //[nMC]
   Float_t         mcEta[12];   //[nMC]
   Float_t         mcPhi[12];   //[nMC]
   Float_t         mcE[12];   //[nMC]
   Float_t         mcEt[12];   //[nMC]
   Int_t           mcGMomPID[12];   //[nMC]
   Int_t           mcMomPID[12];   //[nMC]
   Float_t         mcMomPt[12];   //[nMC]
   Float_t         mcMomMass[12];   //[nMC]
   Float_t         mcMomEta[12];   //[nMC]
   Float_t         mcMomPhi[12];   //[nMC]
   Int_t           mcIndex[12];   //[nMC]
   Int_t           mcDecayType[12];   //[nMC]
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           nPUInfo;
   Int_t           nPU[3];   //[nPUInfo]
   Int_t           puBX[3];   //[nPUInfo]
   Float_t         puTrue[3];   //[nPUInfo]
   Float_t         MET;
   Float_t         METPhi;
   Float_t         METsumEt;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Int_t           metFilters[10];
   Int_t           nEle;
   Int_t           eleTrg[8][16];   //[nEle]
   Int_t           eleClass[8];   //[nEle]
   Int_t           eleIsEcalDriven[8];   //[nEle]
   Int_t           eleCharge[8];   //[nEle]
   Float_t         eleEn[8];   //[nEle]
   Float_t         eleEcalEn[8];   //[nEle]
   Float_t         eleSCRawEn[8];   //[nEle]
   Float_t         eleSCEn[8];   //[nEle]
   Float_t         eleESEn[8];   //[nEle]
   Float_t         elePt[8];   //[nEle]
   Float_t         eleEta[8];   //[nEle]
   Float_t         elePhi[8];   //[nEle]
   Float_t         eleEtaVtx[8][100];   //[nEle]
   Float_t         elePhiVtx[8][100];   //[nEle]
   Float_t         eleEtVtx[8][100];   //[nEle]
   Float_t         eleSCEta[8];   //[nEle]
   Float_t         eleSCPhi[8];   //[nEle]
   Float_t         eleSCEtaWidth[8];   //[nEle]
   Float_t         eleSCPhiWidth[8];   //[nEle]
   Float_t         eleVtx[8][3];   //[nEle]
   Float_t         eleD0[8];   //[nEle]
   Float_t         eleDz[8];   //[nEle]
   Float_t         eleD0GV[8];   //[nEle]
   Float_t         eleDzGV[8];   //[nEle]
   Float_t         eleD0Vtx[8][100];   //[nEle]
   Float_t         eleDzVtx[8][100];   //[nEle]
   Float_t         eleHoverE[8];   //[nEle]
   Float_t         eleHoverE12[8];   //[nEle]
   Float_t         eleEoverP[8];   //[nEle]
   Float_t         elePin[8];   //[nEle]
   Float_t         elePout[8];   //[nEle]
   Float_t         eleTrkMomErr[8];   //[nEle]
   Float_t         eleBrem[8];   //[nEle]
   Float_t         eledEtaAtVtx[8];   //[nEle]
   Float_t         eledPhiAtVtx[8];   //[nEle]
   Float_t         eleSigmaIEtaIEta[8];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[8];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[8];   //[nEle]
   Float_t         eleEmax[8];   //[nEle]
   Float_t         eleE1x5[8];   //[nEle]
   Float_t         eleE3x3[8];   //[nEle]
   Float_t         eleE5x5[8];   //[nEle]
   Float_t         eleE2x5Right[8];   //[nEle]
   Float_t         eleE2x5Left[8];   //[nEle]
   Float_t         eleE2x5Top[8];   //[nEle]
   Float_t         eleE2x5Bottom[8];   //[nEle]
   Float_t         eleRegrE[8];   //[nEle]
   Float_t         eleRegrEerr[8];   //[nEle]
   Float_t         elePhoRegrE[8];   //[nEle]
   Float_t         elePhoRegrEerr[8];   //[nEle]
   Float_t         eleSeedTime[8];   //[nEle]
   Int_t           eleRecoFlag[8];   //[nEle]
   Int_t           elePos[8];   //[nEle]
   Int_t           eleGenIndex[8];   //[nEle]
   Int_t           eleGenGMomPID[8];   //[nEle]
   Int_t           eleGenMomPID[8];   //[nEle]
   Float_t         eleGenMomPt[8];   //[nEle]
   Float_t         eleIsoTrkDR03[8];   //[nEle]
   Float_t         eleIsoEcalDR03[8];   //[nEle]
   Float_t         eleIsoHcalDR03[8];   //[nEle]
   Float_t         eleIsoHcalDR0312[8];   //[nEle]
   Float_t         eleIsoTrkDR04[8];   //[nEle]
   Float_t         eleIsoEcalDR04[8];   //[nEle]
   Float_t         eleIsoHcalDR04[8];   //[nEle]
   Float_t         eleIsoHcalDR0412[8];   //[nEle]
   Int_t           eleMissHits[8];   //[nEle]
   Float_t         eleConvDist[8];   //[nEle]
   Float_t         eleConvDcot[8];   //[nEle]
   Int_t           eleConvVtxFit[8];   //[nEle]
   Float_t         eleIP3D[8];   //[nEle]
   Float_t         eleIP3DErr[8];   //[nEle]
   Float_t         eleIDMVANonTrig[8];   //[nEle]
   Float_t         eleIDMVATrig[8];   //[nEle]
   Float_t         eleIDMVATrigIDIso[8];   //[nEle]
   Float_t         elePFChIso03[8];   //[nEle]
   Float_t         elePFPhoIso03[8];   //[nEle]
   Float_t         elePFNeuIso03[8];   //[nEle]
   Float_t         elePFChIso04[8];   //[nEle]
   Float_t         elePFPhoIso04[8];   //[nEle]
   Float_t         elePFNeuIso04[8];   //[nEle]
   Float_t         eleESEffSigmaRR[8][3];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[10][8];   //[nPho]
   Int_t           phoTrgFilter[10][50];   //[nPho]
   Bool_t          phoIsPhoton[10];   //[nPho]
   Float_t         phoSCPos[10][3];   //[nPho]
   Float_t         phoCaloPos[10][3];   //[nPho]
   Float_t         phoE[10];   //[nPho]
   Float_t         phoEt[10];   //[nPho]
   Float_t         phoEta[10];   //[nPho]
   Float_t         phoVtx[10][3];   //[nPho]
   Float_t         phoPhi[10];   //[nPho]
   Float_t         phoEtVtx[10][100];   //[nPho]
   Float_t         phoEtaVtx[10][100];   //[nPho]
   Float_t         phoPhiVtx[10][100];   //[nPho]
   Float_t         phoR9[10];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[10];   //[nPho]
   Float_t         phoEcalIsoDR03[10];   //[nPho]
   Float_t         phoHcalIsoDR03[10];   //[nPho]
   Float_t         phoHcalIsoDR0312[10];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[10];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[10][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[10][100];   //[nPho]
   Float_t         phoCiCdRtoTrk[10];   //[nPho]
   Float_t         phoEcalIsoDR04[10];   //[nPho]
   Float_t         phoHcalIsoDR04[10];   //[nPho]
   Float_t         phoHcalIsoDR0412[10];   //[nPho]
   Float_t         phoHoverE[10];   //[nPho]
   Float_t         phoHoverE12[10];   //[nPho]
   Int_t           phoEleVeto[10];   //[nPho]
   Float_t         phoSigmaIEtaIEta[10];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[10];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[10];   //[nPho]
   Float_t         phoCiCPF4phopfIso03[10];   //[nPho]
   Float_t         phoCiCPF4phopfIso04[10];   //[nPho]
   Float_t         phoCiCPF4chgpfIso02[10][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso03[10][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso04[10][100];   //[nPho]
   Float_t         phoEmax[10];   //[nPho]
   Float_t         phoEtop[10];   //[nPho]
   Float_t         phoEbottom[10];   //[nPho]
   Float_t         phoEleft[10];   //[nPho]
   Float_t         phoEright[10];   //[nPho]
   Float_t         phoE3x3[10];   //[nPho]
   Float_t         phoE3x1[10];   //[nPho]
   Float_t         phoE1x3[10];   //[nPho]
   Float_t         phoE5x5[10];   //[nPho]
   Float_t         phoE1x5[10];   //[nPho]
   Float_t         phoE2x2[10];   //[nPho]
   Float_t         phoE2x5Max[10];   //[nPho]
   Float_t         phoE2x5Right[10];   //[nPho]
   Float_t         phoE2x5Left[10];   //[nPho]
   Float_t         phoE2x5Top[10];   //[nPho]
   Float_t         phoE2x5Bottom[10];   //[nPho]
   Float_t         phoPFChIso[10];   //[nPho]
   Float_t         phoPFPhoIso[10];   //[nPho]
   Float_t         phoPFNeuIso[10];   //[nPho]
   Float_t         phoRegrE[10];   //[nPho]
   Float_t         phoRegrEerr[10];   //[nPho]
   Float_t         phoSeedTime[10];   //[nPho]
   Int_t           phoSeedDetId1[10];   //[nPho]
   Int_t           phoSeedDetId2[10];   //[nPho]
   Int_t           phoRecoFlag[10];   //[nPho]
   Int_t           phoPos[10];   //[nPho]
   Int_t           phoGenIndex[10];   //[nPho]
   Int_t           phoGenGMomPID[10];   //[nPho]
   Int_t           phoGenMomPID[10];   //[nPho]
   Float_t         phoGenMomPt[10];   //[nPho]
   Float_t         phoSCE[10];   //[nPho]
   Float_t         phoSCRawE[10];   //[nPho]
   Float_t         phoESEn[10];   //[nPho]
   Float_t         phoSCEt[10];   //[nPho]
   Float_t         phoSCEta[10];   //[nPho]
   Float_t         phoSCPhi[10];   //[nPho]
   Float_t         phoSCEtaWidth[10];   //[nPho]
   Float_t         phoSCPhiWidth[10];   //[nPho]
   Float_t         phoSCBrem[10];   //[nPho]
   Int_t           phoOverlap[10];   //[nPho]
   Int_t           phohasPixelSeed[10];   //[nPho]
   Int_t           pho_hasConvPf[10];   //[nPho]
   Int_t           pho_hasSLConvPf[10];   //[nPho]
   Float_t         pho_pfconvVtxZ[10];   //[nPho]
   Float_t         pho_pfconvVtxZErr[10];   //[nPho]
   Int_t           pho_nSLConv[10];   //[nPho]
   Float_t         pho_pfSLConvPos[10][3];   //[nPho]
   Float_t         pho_pfSLConvVtxZ[10][20];   //[nPho]
   Int_t           phoIsConv[10];   //[nPho]
   Int_t           phoNConv[10];   //[nPho]
   Float_t         phoConvInvMass[10];   //[nPho]
   Float_t         phoConvCotTheta[10];   //[nPho]
   Float_t         phoConvEoverP[10];   //[nPho]
   Float_t         phoConvZofPVfromTrks[10];   //[nPho]
   Float_t         phoConvMinDist[10];   //[nPho]
   Float_t         phoConvdPhiAtVtx[10];   //[nPho]
   Float_t         phoConvdPhiAtCalo[10];   //[nPho]
   Float_t         phoConvdEtaAtCalo[10];   //[nPho]
   Float_t         phoConvTrkd0[10][2];   //[nPho]
   Float_t         phoConvTrkPin[10][2];   //[nPho]
   Float_t         phoConvTrkPout[10][2];   //[nPho]
   Float_t         phoConvTrkdz[10][2];   //[nPho]
   Float_t         phoConvTrkdzErr[10][2];   //[nPho]
   Float_t         phoConvChi2[10];   //[nPho]
   Float_t         phoConvChi2Prob[10];   //[nPho]
   Int_t           phoConvNTrks[10];   //[nPho]
   Float_t         phoConvCharge[10][2];   //[nPho]
   Float_t         phoConvValidVtx[10];   //[nPho]
   Float_t         phoConvLikeLihood[10];   //[nPho]
   Float_t         phoConvP4[10][4];   //[nPho]
   Float_t         phoConvVtx[10][3];   //[nPho]
   Float_t         phoConvVtxErr[10][3];   //[nPho]
   Float_t         phoConvPairMomentum[10][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[10][3];   //[nPho]
   Int_t           SingleLegConv[10];   //[nPho]
   Float_t         phoPFConvVtx[10][3];   //[nPho]
   Float_t         phoPFConvMom[10][3];   //[nPho]
   Float_t         phoESEffSigmaRR[10][3];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[16][10];   //[nMu]
   Float_t         muEta[16];   //[nMu]
   Float_t         muPhi[16];   //[nMu]
   Int_t           muCharge[16];   //[nMu]
   Float_t         muPt[16];   //[nMu]
   Float_t         muPz[16];   //[nMu]
   Float_t         muVtx[16][3];   //[nMu]
   Float_t         muVtxGlb[16][3];   //[nMu]
   Int_t           muGenIndex[16];   //[nMu]
   Float_t         mucktPt_[16];   //[nMu]
   Float_t         mucktPtErr_[16];   //[nMu]
   Float_t         mucktdxy_[16];   //[nMu]
   Float_t         mucktdz_[16];   //[nMu]
   Float_t         muIsoTrk[16];   //[nMu]
   Float_t         muIsoCalo[16];   //[nMu]
   Float_t         muIsoEcal[16];   //[nMu]
   Float_t         muIsoHcal[16];   //[nMu]
   Float_t         muChi2NDF[16];   //[nMu]
   Float_t         muInnerChi2NDF[16];   //[nMu]
   Float_t         muPFIsoR04_CH[16];   //[nMu]
   Float_t         muPFIsoR04_NH[16];   //[nMu]
   Float_t         muPFIsoR04_Pho[16];   //[nMu]
   Float_t         muPFIsoR04_PU[16];   //[nMu]
   Float_t         muPFIsoR04_CPart[16];   //[nMu]
   Float_t         muPFIsoR04_NHHT[16];   //[nMu]
   Float_t         muPFIsoR04_PhoHT[16];   //[nMu]
   Float_t         muPFIsoR03_CH[16];   //[nMu]
   Float_t         muPFIsoR03_NH[16];   //[nMu]
   Float_t         muPFIsoR03_Pho[16];   //[nMu]
   Float_t         muPFIsoR03_PU[16];   //[nMu]
   Float_t         muPFIsoR03_CPart[16];   //[nMu]
   Float_t         muPFIsoR03_NHHT[16];   //[nMu]
   Float_t         muPFIsoR03_PhoHT[16];   //[nMu]
   Int_t           muType[16];   //[nMu]
   Bool_t          muID[16][6];   //[nMu]
   Float_t         muD0[16];   //[nMu]
   Float_t         muDz[16];   //[nMu]
   Float_t         muD0GV[16];   //[nMu]
   Float_t         muDzGV[16];   //[nMu]
   Float_t         muD0Vtx[16][100];   //[nMu]
   Float_t         muDzVtx[16][100];   //[nMu]
   Float_t         muInnerD0[16];   //[nMu]
   Float_t         muInnerDz[16];   //[nMu]
   Float_t         muInnerD0GV[16];   //[nMu]
   Float_t         muInnerDzGV[16];   //[nMu]
   Int_t           muNumberOfValidTrkLayers[16];   //[nMu]
   Int_t           muNumberOfValidTrkHits[16];   //[nMu]
   Int_t           muNumberOfValidPixelLayers[16];   //[nMu]
   Int_t           muNumberOfValidPixelHits[16];   //[nMu]
   Int_t           muNumberOfValidMuonHits[16];   //[nMu]
   Int_t           muStations[16];   //[nMu]
   Int_t           muChambers[16];   //[nMu]
   Float_t         muIP3D[16];   //[nMu]
   Float_t         muIP3DErr[16];   //[nMu]
   Int_t           nPFEle_;
   Float_t         PFElePt_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFEleEta_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFElePhi_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFEleEn_[kMaxnPFEle];   //[nPFEle_]
   Int_t           PFEleCharge[kMaxnPFEle];   //[nPFEle_]
   Float_t         rho25;
   Float_t         rho25_neu;
   Float_t         rho25_muPFiso;
   Float_t         rho25_elePFiso;
   Float_t         rho2011;
   Float_t         rho2012;
   Int_t           nJet;
   Int_t           jetTrg[92][14];   //[nJet]
   Float_t         jetEn[92];   //[nJet]
   Float_t         jetPt[92];   //[nJet]
   Float_t         jetEta[92];   //[nJet]
   Float_t         jetPhi[92];   //[nJet]
   Float_t         jetCharge[92];   //[nJet]
   Float_t         jetEt[92];   //[nJet]
   Float_t         jetRawPt[92];   //[nJet]
   Float_t         jetRawEn[92];   //[nJet]
   Float_t         jetArea[92];   //[nJet]
   Float_t         jetCHF[92];   //[nJet]
   Float_t         jetNHF[92];   //[nJet]
   Float_t         jetCEF[92];   //[nJet]
   Float_t         jetNEF[92];   //[nJet]
   Int_t           jetNCH[92];   //[nJet]
   Float_t         jetHFHAE[92];   //[nJet]
   Float_t         jetHFEME[92];   //[nJet]
   Int_t           jetNConstituents[92];   //[nJet]
   Float_t         jetCombinedSecondaryVtxBJetTags[92];   //[nJet]
   Float_t         jetCombinedSecondaryVtxMVABJetTags[92];   //[nJet]
   Float_t         jetJetProbabilityBJetTags[92];   //[nJet]
   Float_t         jetJetBProbabilityBJetTags[92];   //[nJet]
   Float_t         jetTrackCountingHighPurBJetTags[92];   //[nJet]
   Float_t         jetBetaStar[92][100];   //[nJet]
   Bool_t          jetPFLooseId[92];   //[nJet]
   Float_t         jetDRMean[92];   //[nJet]
   Float_t         jetDR2Mean[92];   //[nJet]
   Float_t         jetDZ[92];   //[nJet]
   Float_t         jetFrac01[92];   //[nJet]
   Float_t         jetFrac02[92];   //[nJet]
   Float_t         jetFrac03[92];   //[nJet]
   Float_t         jetFrac04[92];   //[nJet]
   Float_t         jetFrac05[92];   //[nJet]
   Float_t         jetFrac06[92];   //[nJet]
   Float_t         jetFrac07[92];   //[nJet]
   Float_t         jetBeta[92];   //[nJet]
   Float_t         jetBetaStarCMG[92];   //[nJet]
   Float_t         jetBetaStarClassic[92];   //[nJet]
   Float_t         jetBetaExt[92][100];   //[nJet]
   Float_t         jetBetaStarCMGExt[92][100];   //[nJet]
   Float_t         jetBetaStarClassicExt[92][100];   //[nJet]
   Float_t         jetNNeutrals[92];   //[nJet]
   Float_t         jetNCharged[92];   //[nJet]
   Float_t         jetMVAs[92][4];   //[nJet]
   Int_t           jetWPLevels[92][4];   //[nJet]
   Float_t         jetMVAsExt[92][4][100];   //[nJet]
   Int_t           jetWPLevelsExt[92][4][100];   //[nJet]
   Int_t           jetPartonID[92];   //[nJet]
   Int_t           jetGenJetIndex[92];   //[nJet]
   Float_t         jetGenJetEn[92];   //[nJet]
   Float_t         jetGenJetPt[92];   //[nJet]
   Float_t         jetGenJetEta[92];   //[nJet]
   Float_t         jetGenJetPhi[92];   //[nJet]
   Int_t           jetGenPartonID[92];   //[nJet]
   Float_t         jetGenEn[92];   //[nJet]
   Float_t         jetGenPt[92];   //[nJet]
   Float_t         jetGenEta[92];   //[nJet]
   Float_t         jetGenPhi[92];   //[nJet]
   Int_t           nLowPtJet;
   Float_t         jetLowPtEn[61];   //[nLowPtJet]
   Float_t         jetLowPtPt[61];   //[nLowPtJet]
   Float_t         jetLowPtEta[61];   //[nLowPtJet]
   Float_t         jetLowPtPhi[61];   //[nLowPtJet]
   Float_t         jetLowPtCharge[61];   //[nLowPtJet]
   Float_t         jetLowPtEt[61];   //[nLowPtJet]
   Float_t         jetLowPtRawPt[61];   //[nLowPtJet]
   Float_t         jetLowPtRawEn[61];   //[nLowPtJet]
   Float_t         jetLowPtArea[61];   //[nLowPtJet]
   Int_t           jetLowPtPartonID[61];   //[nLowPtJet]
   Float_t         jetLowPtGenJetEn[61];   //[nLowPtJet]
   Float_t         jetLowPtGenJetPt[61];   //[nLowPtJet]
   Float_t         jetLowPtGenJetEta[61];   //[nLowPtJet]
   Float_t         jetLowPtGenJetPhi[61];   //[nLowPtJet]
   Int_t           jetLowPtGenPartonID[61];   //[nLowPtJet]
   Float_t         jetLowPtGenEn[61];   //[nLowPtJet]
   Float_t         jetLowPtGenPt[61];   //[nLowPtJet]
   Float_t         jetLowPtGenEta[61];   //[nLowPtJet]
   Float_t         jetLowPtGenPhi[61];   //[nLowPtJet]
   Int_t           nConv;
   Float_t         convP4[147][4];   //[nConv]
   Float_t         convVtx[147][3];   //[nConv]
   Float_t         convVtxErr[147][3];   //[nConv]
   Float_t         convPairMomentum[147][3];   //[nConv]
   Float_t         convRefittedMomentum[147][3];   //[nConv]
   Int_t           convNTracks[147];   //[nConv]
   Float_t         convPairInvMass[147];   //[nConv]
   Float_t         convPairCotThetaSep[147];   //[nConv]
   Float_t         convEoverP[147];   //[nConv]
   Float_t         convDistOfMinApproach[147];   //[nConv]
   Float_t         convDPhiTrksAtVtx[147];   //[nConv]
   Float_t         convDPhiTrksAtEcal[147];   //[nConv]
   Float_t         convDEtaTrksAtEcal[147];   //[nConv]
   Float_t         convDxy[147];   //[nConv]
   Float_t         convDz[147];   //[nConv]
   Float_t         convLxy[147];   //[nConv]
   Float_t         convLz[147];   //[nConv]
   Float_t         convZofPrimVtxFromTrks[147];   //[nConv]
   Int_t           convNHitsBeforeVtx[147][2];   //[nConv]
   Int_t           convNSharedHits[147];   //[nConv]
   Int_t           convValidVtx[147];   //[nConv]
   Float_t         convMVALikelihood[147];   //[nConv]
   Float_t         convChi2[147];   //[nConv]
   Float_t         convChi2Probability[147];   //[nConv]
   Float_t         convTk1Dz[147];   //[nConv]
   Float_t         convTk2Dz[147];   //[nConv]
   Float_t         convTk1DzErr[147];   //[nConv]
   Float_t         convTk2DzErr[147];   //[nConv]
   Int_t           convCh1Ch2[147];   //[nConv]
   Float_t         convTk1D0[147];   //[nConv]
   Float_t         convTk1Pout[147];   //[nConv]
   Float_t         convTk1Pin[147];   //[nConv]
   Float_t         convTk2D0[147];   //[nConv]
   Float_t         convTk2Pout[147];   //[nConv]
   Float_t         convTk2Pin[147];   //[nConv]

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
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METsumEt;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleIsEcalDriven;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleEtVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleD0GV;   //!
   TBranch        *b_eleDzGV;   //!
   TBranch        *b_eleD0Vtx;   //!
   TBranch        *b_eleDzVtx;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleHoverE12;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleTrkMomErr;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleRegrE;   //!
   TBranch        *b_eleRegrEerr;   //!
   TBranch        *b_elePhoRegrE;   //!
   TBranch        *b_elePhoRegrEerr;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_elePos;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleGenGMomPID;   //!
   TBranch        *b_eleGenMomPID;   //!
   TBranch        *b_eleGenMomPt;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalDR0312;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalDR0412;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleConvVtxFit;   //!
   TBranch        *b_eleIP3D;   //!
   TBranch        *b_eleIP3DErr;   //!
   TBranch        *b_eleIDMVANonTrig;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_eleIDMVATrigIDIso;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_elePFChIso04;   //!
   TBranch        *b_elePFPhoIso04;   //!
   TBranch        *b_elePFNeuIso04;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoTrgFilter;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoSCPos;   //!
   TBranch        *b_phoCaloPos;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR0312;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoCiCTrkIsoDR03;   //!
   TBranch        *b_phoCiCTrkIsoDR04;   //!
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR0412;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoHoverE12;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoCiCPF4phopfIso03;   //!
   TBranch        *b_phoCiCPF4phopfIso04;   //!
   TBranch        *b_phoCiCPF4chgpfIso02;   //!
   TBranch        *b_phoCiCPF4chgpfIso03;   //!
   TBranch        *b_phoCiCPF4chgpfIso04;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoEtop;   //!
   TBranch        *b_phoEbottom;   //!
   TBranch        *b_phoEleft;   //!
   TBranch        *b_phoEright;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE3x1;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrEerr;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoGenIndex;   //!
   TBranch        *b_phoGenGMomPID;   //!
   TBranch        *b_phoGenMomPID;   //!
   TBranch        *b_phoGenMomPt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_pho_hasConvPf;   //!
   TBranch        *b_pho_hasSLConvPf;   //!
   TBranch        *b_pho_pfconvVtxZ;   //!
   TBranch        *b_pho_pfconvVtxZErr;   //!
   TBranch        *b_pho_nSLConv;   //!
   TBranch        *b_pho_pfSLConvPos;   //!
   TBranch        *b_pho_pfSLConvVtxZ;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0;   //!
   TBranch        *b_phoConvTrkPin;   //!
   TBranch        *b_phoConvTrkPout;   //!
   TBranch        *b_phoConvTrkdz;   //!
   TBranch        *b_phoConvTrkdzErr;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvNTrks;   //!
   TBranch        *b_phoConvCharge;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4;   //!
   TBranch        *b_phoConvVtx;   //!
   TBranch        *b_phoConvVtxErr;   //!
   TBranch        *b_phoConvPairMomentum;   //!
   TBranch        *b_phoConvRefittedMomentum;   //!
   TBranch        *b_SingleLegConv;   //!
   TBranch        *b_phoPFConvVtx;   //!
   TBranch        *b_phoPFConvMom;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muVtx;   //!
   TBranch        *b_muVtxGlb;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_mucktPt_;   //!
   TBranch        *b_mucktPtErr_;   //!
   TBranch        *b_mucktdxy_;   //!
   TBranch        *b_mucktdz_;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerChi2NDF;   //!
   TBranch        *b_muPFIsoR04_CH;   //!
   TBranch        *b_muPFIsoR04_NH;   //!
   TBranch        *b_muPFIsoR04_Pho;   //!
   TBranch        *b_muPFIsoR04_PU;   //!
   TBranch        *b_muPFIsoR04_CPart;   //!
   TBranch        *b_muPFIsoR04_NHHT;   //!
   TBranch        *b_muPFIsoR04_PhoHT;   //!
   TBranch        *b_muPFIsoR03_CH;   //!
   TBranch        *b_muPFIsoR03_NH;   //!
   TBranch        *b_muPFIsoR03_Pho;   //!
   TBranch        *b_muPFIsoR03_PU;   //!
   TBranch        *b_muPFIsoR03_CPart;   //!
   TBranch        *b_muPFIsoR03_NHHT;   //!
   TBranch        *b_muPFIsoR03_PhoHT;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muD0GV;   //!
   TBranch        *b_muDzGV;   //!
   TBranch        *b_muD0Vtx;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muInnerD0GV;   //!
   TBranch        *b_muInnerDzGV;   //!
   TBranch        *b_muNumberOfValidTrkLayers;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelLayers;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_muIP3D;   //!
   TBranch        *b_muIP3DErr;   //!
   TBranch        *b_nPFEle_;   //!
   TBranch        *b_PFElePt_;   //!
   TBranch        *b_PFEleEta_;   //!
   TBranch        *b_PFElePhi_;   //!
   TBranch        *b_PFEleEn_;   //!
   TBranch        *b_PFEleCharge;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho25_neu;   //!
   TBranch        *b_rho25_muPFiso;   //!
   TBranch        *b_rho25_elePFiso;   //!
   TBranch        *b_rho2011;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetCombinedSecondaryVtxBJetTags;   //!
   TBranch        *b_jetCombinedSecondaryVtxMVABJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetJetBProbabilityBJetTags;   //!
   TBranch        *b_jetTrackCountingHighPurBJetTags;   //!
   TBranch        *b_jetBetaStar;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetDRMean;   //!
   TBranch        *b_jetDR2Mean;   //!
   TBranch        *b_jetDZ;   //!
   TBranch        *b_jetFrac01;   //!
   TBranch        *b_jetFrac02;   //!
   TBranch        *b_jetFrac03;   //!
   TBranch        *b_jetFrac04;   //!
   TBranch        *b_jetFrac05;   //!
   TBranch        *b_jetFrac06;   //!
   TBranch        *b_jetFrac07;   //!
   TBranch        *b_jetBeta;   //!
   TBranch        *b_jetBetaStarCMG;   //!
   TBranch        *b_jetBetaStarClassic;   //!
   TBranch        *b_jetBetaExt;   //!
   TBranch        *b_jetBetaStarCMGExt;   //!
   TBranch        *b_jetBetaStarClassicExt;   //!
   TBranch        *b_jetNNeutrals;   //!
   TBranch        *b_jetNCharged;   //!
   TBranch        *b_jetMVAs;   //!
   TBranch        *b_jetWPLevels;   //!
   TBranch        *b_jetMVAsExt;   //!
   TBranch        *b_jetWPLevelsExt;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_nLowPtJet;   //!
   TBranch        *b_jetLowPtEn;   //!
   TBranch        *b_jetLowPtPt;   //!
   TBranch        *b_jetLowPtEta;   //!
   TBranch        *b_jetLowPtPhi;   //!
   TBranch        *b_jetLowPtCharge;   //!
   TBranch        *b_jetLowPtEt;   //!
   TBranch        *b_jetLowPtRawPt;   //!
   TBranch        *b_jetLowPtRawEn;   //!
   TBranch        *b_jetLowPtArea;   //!
   TBranch        *b_jetLowPtPartonID;   //!
   TBranch        *b_jetLowPtGenJetEn;   //!
   TBranch        *b_jetLowPtGenJetPt;   //!
   TBranch        *b_jetLowPtGenJetEta;   //!
   TBranch        *b_jetLowPtGenJetPhi;   //!
   TBranch        *b_jetLowPtGenPartonID;   //!
   TBranch        *b_jetLowPtGenEn;   //!
   TBranch        *b_jetLowPtGenPt;   //!
   TBranch        *b_jetLowPtGenEta;   //!
   TBranch        *b_jetLowPtGenPhi;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4;   //!
   TBranch        *b_convVtx;   //!
   TBranch        *b_convVtxErr;   //!
   TBranch        *b_convPairMomentum;   //!
   TBranch        *b_convRefittedMomentum;   //!
   TBranch        *b_convNTracks;   //!
   TBranch        *b_convPairInvMass;   //!
   TBranch        *b_convPairCotThetaSep;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convDistOfMinApproach;   //!
   TBranch        *b_convDPhiTrksAtVtx;   //!
   TBranch        *b_convDPhiTrksAtEcal;   //!
   TBranch        *b_convDEtaTrksAtEcal;   //!
   TBranch        *b_convDxy;   //!
   TBranch        *b_convDz;   //!
   TBranch        *b_convLxy;   //!
   TBranch        *b_convLz;   //!
   TBranch        *b_convZofPrimVtxFromTrks;   //!
   TBranch        *b_convNHitsBeforeVtx;   //!
   TBranch        *b_convNSharedHits;   //!
   TBranch        *b_convValidVtx;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convChi2;   //!
   TBranch        *b_convChi2Probability;   //!
   TBranch        *b_convTk1Dz;   //!
   TBranch        *b_convTk2Dz;   //!
   TBranch        *b_convTk1DzErr;   //!
   TBranch        *b_convTk2DzErr;   //!
   TBranch        *b_convCh1Ch2;   //!
   TBranch        *b_convTk1D0;   //!
   TBranch        *b_convTk1Pout;   //!
   TBranch        *b_convTk1Pin;   //!
   TBranch        *b_convTk2D0;   //!
   TBranch        *b_convTk2Pout;   //!
   TBranch        *b_convTk2Pin;   //!

  ttgamma3(TTree * /*tree*/ =0) : fProofFile(0),h1test(0),hPU_weights(0),fFile(0),fChain(0) 
  { 
    fChannel = 2; //default e+jets
    fVerbose = false;
    fIsMC    = false;
    fPUreweighting = false;
    fSample  = "";
    fOutdir  = "";
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
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs", vtxbs, &b_vtxbs);
   fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcDecayType", mcDecayType, &b_mcDecayType);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", puTrue, &b_puTrue);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METsumEt", &METsumEt, &b_METsumEt);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("metFilters", metFilters, &b_metFilters);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleIsEcalDriven", eleIsEcalDriven, &b_eleIsEcalDriven);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleEcalEn", eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleEtaVtx", eleEtaVtx, &b_eleEtaVtx);
   fChain->SetBranchAddress("elePhiVtx", elePhiVtx, &b_elePhiVtx);
   fChain->SetBranchAddress("eleEtVtx", eleEtVtx, &b_eleEtVtx);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
   fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleD0GV", eleD0GV, &b_eleD0GV);
   fChain->SetBranchAddress("eleDzGV", eleDzGV, &b_eleDzGV);
   fChain->SetBranchAddress("eleD0Vtx", eleD0Vtx, &b_eleD0Vtx);
   fChain->SetBranchAddress("eleDzVtx", eleDzVtx, &b_eleDzVtx);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleHoverE12", eleHoverE12, &b_eleHoverE12);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", elePout, &b_elePout);
   fChain->SetBranchAddress("eleTrkMomErr", eleTrkMomErr, &b_eleTrkMomErr);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleEmax", eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE1x5", eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Right", eleE2x5Right, &b_eleE2x5Right);
   fChain->SetBranchAddress("eleE2x5Left", eleE2x5Left, &b_eleE2x5Left);
   fChain->SetBranchAddress("eleE2x5Top", eleE2x5Top, &b_eleE2x5Top);
   fChain->SetBranchAddress("eleE2x5Bottom", eleE2x5Bottom, &b_eleE2x5Bottom);
   fChain->SetBranchAddress("eleRegrE", eleRegrE, &b_eleRegrE);
   fChain->SetBranchAddress("eleRegrEerr", eleRegrEerr, &b_eleRegrEerr);
   fChain->SetBranchAddress("elePhoRegrE", elePhoRegrE, &b_elePhoRegrE);
   fChain->SetBranchAddress("elePhoRegrEerr", elePhoRegrEerr, &b_elePhoRegrEerr);
   fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
   fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("elePos", elePos, &b_elePos);
   fChain->SetBranchAddress("eleGenIndex", eleGenIndex, &b_eleGenIndex);
   fChain->SetBranchAddress("eleGenGMomPID", eleGenGMomPID, &b_eleGenGMomPID);
   fChain->SetBranchAddress("eleGenMomPID", eleGenMomPID, &b_eleGenMomPID);
   fChain->SetBranchAddress("eleGenMomPt", eleGenMomPt, &b_eleGenMomPt);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR0312", eleIsoHcalDR0312, &b_eleIsoHcalDR0312);
   fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR0412", eleIsoHcalDR0412, &b_eleIsoHcalDR0412);
   fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleConvVtxFit", eleConvVtxFit, &b_eleConvVtxFit);
   fChain->SetBranchAddress("eleIP3D", eleIP3D, &b_eleIP3D);
   fChain->SetBranchAddress("eleIP3DErr", eleIP3DErr, &b_eleIP3DErr);
   fChain->SetBranchAddress("eleIDMVANonTrig", eleIDMVANonTrig, &b_eleIDMVANonTrig);
   fChain->SetBranchAddress("eleIDMVATrig", eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("eleIDMVATrigIDIso", eleIDMVATrigIDIso, &b_eleIDMVATrigIDIso);
   fChain->SetBranchAddress("elePFChIso03", elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("elePFChIso04", elePFChIso04, &b_elePFChIso04);
   fChain->SetBranchAddress("elePFPhoIso04", elePFPhoIso04, &b_elePFPhoIso04);
   fChain->SetBranchAddress("elePFNeuIso04", elePFNeuIso04, &b_elePFNeuIso04);
   fChain->SetBranchAddress("eleESEffSigmaRR", eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoTrgFilter", phoTrgFilter, &b_phoTrgFilter);
   fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoSCPos", phoSCPos, &b_phoSCPos);
   fChain->SetBranchAddress("phoCaloPos", phoCaloPos, &b_phoCaloPos);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx", phoVtx, &b_phoVtx);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR0312", phoHcalIsoDR0312, &b_phoHcalIsoDR0312);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03, &b_phoCiCTrkIsoDR03);
   fChain->SetBranchAddress("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04, &b_phoCiCTrkIsoDR04);
   fChain->SetBranchAddress("phoCiCdRtoTrk", phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR0412", phoHcalIsoDR0412, &b_phoHcalIsoDR0412);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoHoverE12", phoHoverE12, &b_phoHoverE12);
   fChain->SetBranchAddress("phoEleVeto", phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoCiCPF4phopfIso03", phoCiCPF4phopfIso03, &b_phoCiCPF4phopfIso03);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04", phoCiCPF4phopfIso04, &b_phoCiCPF4phopfIso04);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso02", phoCiCPF4chgpfIso02, &b_phoCiCPF4chgpfIso02);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso03", phoCiCPF4chgpfIso03, &b_phoCiCPF4chgpfIso03);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04", phoCiCPF4chgpfIso04, &b_phoCiCPF4chgpfIso04);
   fChain->SetBranchAddress("phoEmax", phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoEtop", phoEtop, &b_phoEtop);
   fChain->SetBranchAddress("phoEbottom", phoEbottom, &b_phoEbottom);
   fChain->SetBranchAddress("phoEleft", phoEleft, &b_phoEleft);
   fChain->SetBranchAddress("phoEright", phoEright, &b_phoEright);
   fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE3x1", phoE3x1, &b_phoE3x1);
   fChain->SetBranchAddress("phoE1x3", phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE5x5", phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE1x5", phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE2x5Right", phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoE2x5Left", phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Top", phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoPFChIso", phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoRegrE", phoRegrE, &b_phoRegrE);
   fChain->SetBranchAddress("phoRegrEerr", phoRegrEerr, &b_phoRegrEerr);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoGenIndex", phoGenIndex, &b_phoGenIndex);
   fChain->SetBranchAddress("phoGenGMomPID", phoGenGMomPID, &b_phoGenGMomPID);
   fChain->SetBranchAddress("phoGenMomPID", phoGenMomPID, &b_phoGenMomPID);
   fChain->SetBranchAddress("phoGenMomPt", phoGenMomPt, &b_phoGenMomPt);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("pho_hasConvPf", pho_hasConvPf, &b_pho_hasConvPf);
   fChain->SetBranchAddress("pho_hasSLConvPf", pho_hasSLConvPf, &b_pho_hasSLConvPf);
   fChain->SetBranchAddress("pho_pfconvVtxZ", pho_pfconvVtxZ, &b_pho_pfconvVtxZ);
   fChain->SetBranchAddress("pho_pfconvVtxZErr", pho_pfconvVtxZErr, &b_pho_pfconvVtxZErr);
   fChain->SetBranchAddress("pho_nSLConv", pho_nSLConv, &b_pho_nSLConv);
   fChain->SetBranchAddress("pho_pfSLConvPos", pho_pfSLConvPos, &b_pho_pfSLConvPos);
   fChain->SetBranchAddress("pho_pfSLConvVtxZ", pho_pfSLConvVtxZ, &b_pho_pfSLConvVtxZ);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0", phoConvTrkd0, &b_phoConvTrkd0);
   fChain->SetBranchAddress("phoConvTrkPin", phoConvTrkPin, &b_phoConvTrkPin);
   fChain->SetBranchAddress("phoConvTrkPout", phoConvTrkPout, &b_phoConvTrkPout);
   fChain->SetBranchAddress("phoConvTrkdz", phoConvTrkdz, &b_phoConvTrkdz);
   fChain->SetBranchAddress("phoConvTrkdzErr", phoConvTrkdzErr, &b_phoConvTrkdzErr);
   fChain->SetBranchAddress("phoConvChi2", phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvNTrks", phoConvNTrks, &b_phoConvNTrks);
   fChain->SetBranchAddress("phoConvCharge", phoConvCharge, &b_phoConvCharge);
   fChain->SetBranchAddress("phoConvValidVtx", phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4", phoConvP4, &b_phoConvP4);
   fChain->SetBranchAddress("phoConvVtx", phoConvVtx, &b_phoConvVtx);
   fChain->SetBranchAddress("phoConvVtxErr", phoConvVtxErr, &b_phoConvVtxErr);
   fChain->SetBranchAddress("phoConvPairMomentum", phoConvPairMomentum, &b_phoConvPairMomentum);
   fChain->SetBranchAddress("phoConvRefittedMomentum", phoConvRefittedMomentum, &b_phoConvRefittedMomentum);
   fChain->SetBranchAddress("SingleLegConv", SingleLegConv, &b_SingleLegConv);
   fChain->SetBranchAddress("phoPFConvVtx", phoPFConvVtx, &b_phoPFConvVtx);
   fChain->SetBranchAddress("phoPFConvMom", phoPFConvMom, &b_phoPFConvMom);
   fChain->SetBranchAddress("phoESEffSigmaRR", phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", muPz, &b_muPz);
   fChain->SetBranchAddress("muVtx", muVtx, &b_muVtx);
   fChain->SetBranchAddress("muVtxGlb", muVtxGlb, &b_muVtxGlb);
   fChain->SetBranchAddress("muGenIndex", muGenIndex, &b_muGenIndex);
   fChain->SetBranchAddress("mucktPt_", mucktPt_, &b_mucktPt_);
   fChain->SetBranchAddress("mucktPtErr_", mucktPtErr_, &b_mucktPtErr_);
   fChain->SetBranchAddress("mucktdxy_", mucktdxy_, &b_mucktdxy_);
   fChain->SetBranchAddress("mucktdz_", mucktdz_, &b_mucktdz_);
   fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerChi2NDF", muInnerChi2NDF, &b_muInnerChi2NDF);
   fChain->SetBranchAddress("muPFIsoR04_CH", muPFIsoR04_CH, &b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH", muPFIsoR04_NH, &b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho", muPFIsoR04_Pho, &b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("muPFIsoR04_PU", muPFIsoR04_PU, &b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR04_CPart", muPFIsoR04_CPart, &b_muPFIsoR04_CPart);
   fChain->SetBranchAddress("muPFIsoR04_NHHT", muPFIsoR04_NHHT, &b_muPFIsoR04_NHHT);
   fChain->SetBranchAddress("muPFIsoR04_PhoHT", muPFIsoR04_PhoHT, &b_muPFIsoR04_PhoHT);
   fChain->SetBranchAddress("muPFIsoR03_CH", muPFIsoR03_CH, &b_muPFIsoR03_CH);
   fChain->SetBranchAddress("muPFIsoR03_NH", muPFIsoR03_NH, &b_muPFIsoR03_NH);
   fChain->SetBranchAddress("muPFIsoR03_Pho", muPFIsoR03_Pho, &b_muPFIsoR03_Pho);
   fChain->SetBranchAddress("muPFIsoR03_PU", muPFIsoR03_PU, &b_muPFIsoR03_PU);
   fChain->SetBranchAddress("muPFIsoR03_CPart", muPFIsoR03_CPart, &b_muPFIsoR03_CPart);
   fChain->SetBranchAddress("muPFIsoR03_NHHT", muPFIsoR03_NHHT, &b_muPFIsoR03_NHHT);
   fChain->SetBranchAddress("muPFIsoR03_PhoHT", muPFIsoR03_PhoHT, &b_muPFIsoR03_PhoHT);
   fChain->SetBranchAddress("muType", muType, &b_muType);
   fChain->SetBranchAddress("muID", muID, &b_muID);
   fChain->SetBranchAddress("muD0", muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", muDz, &b_muDz);
   fChain->SetBranchAddress("muD0GV", muD0GV, &b_muD0GV);
   fChain->SetBranchAddress("muDzGV", muDzGV, &b_muDzGV);
   fChain->SetBranchAddress("muD0Vtx", muD0Vtx, &b_muD0Vtx);
   fChain->SetBranchAddress("muDzVtx", muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muInnerD0", muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muInnerD0GV", muInnerD0GV, &b_muInnerD0GV);
   fChain->SetBranchAddress("muInnerDzGV", muInnerDzGV, &b_muInnerDzGV);
   fChain->SetBranchAddress("muNumberOfValidTrkLayers", muNumberOfValidTrkLayers, &b_muNumberOfValidTrkLayers);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelLayers", muNumberOfValidPixelLayers, &b_muNumberOfValidPixelLayers);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers);
   fChain->SetBranchAddress("muIP3D", muIP3D, &b_muIP3D);
   fChain->SetBranchAddress("muIP3DErr", muIP3DErr, &b_muIP3DErr);
   fChain->SetBranchAddress("nPFEle_", &nPFEle_, &b_nPFEle_);
   fChain->SetBranchAddress("PFElePt_", PFElePt_, &b_PFElePt_);
   fChain->SetBranchAddress("PFEleEta_", PFEleEta_, &b_PFEleEta_);
   fChain->SetBranchAddress("PFElePhi_", PFElePhi_, &b_PFElePhi_);
   fChain->SetBranchAddress("PFEleEn_", PFEleEn_, &b_PFEleEn_);
   fChain->SetBranchAddress("PFEleCharge", PFEleCharge, &b_PFEleCharge);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho25_neu", &rho25_neu, &b_rho25_neu);
   fChain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso, &b_rho25_muPFiso);
   fChain->SetBranchAddress("rho25_elePFiso", &rho25_elePFiso, &b_rho25_elePFiso);
   fChain->SetBranchAddress("rho2011", &rho2011, &b_rho2011);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", jetCombinedSecondaryVtxMVABJetTags, &b_jetCombinedSecondaryVtxMVABJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetJetBProbabilityBJetTags", jetJetBProbabilityBJetTags, &b_jetJetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetTrackCountingHighPurBJetTags", jetTrackCountingHighPurBJetTags, &b_jetTrackCountingHighPurBJetTags);
   fChain->SetBranchAddress("jetBetaStar", jetBetaStar, &b_jetBetaStar);
   fChain->SetBranchAddress("jetPFLooseId", jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetDRMean", jetDRMean, &b_jetDRMean);
   fChain->SetBranchAddress("jetDR2Mean", jetDR2Mean, &b_jetDR2Mean);
   fChain->SetBranchAddress("jetDZ", jetDZ, &b_jetDZ);
   fChain->SetBranchAddress("jetFrac01", jetFrac01, &b_jetFrac01);
   fChain->SetBranchAddress("jetFrac02", jetFrac02, &b_jetFrac02);
   fChain->SetBranchAddress("jetFrac03", jetFrac03, &b_jetFrac03);
   fChain->SetBranchAddress("jetFrac04", jetFrac04, &b_jetFrac04);
   fChain->SetBranchAddress("jetFrac05", jetFrac05, &b_jetFrac05);
   fChain->SetBranchAddress("jetFrac06", jetFrac06, &b_jetFrac06);
   fChain->SetBranchAddress("jetFrac07", jetFrac07, &b_jetFrac07);
   fChain->SetBranchAddress("jetBeta", jetBeta, &b_jetBeta);
   fChain->SetBranchAddress("jetBetaStarCMG", jetBetaStarCMG, &b_jetBetaStarCMG);
   fChain->SetBranchAddress("jetBetaStarClassic", jetBetaStarClassic, &b_jetBetaStarClassic);
   fChain->SetBranchAddress("jetBetaExt", jetBetaExt, &b_jetBetaExt);
   fChain->SetBranchAddress("jetBetaStarCMGExt", jetBetaStarCMGExt, &b_jetBetaStarCMGExt);
   fChain->SetBranchAddress("jetBetaStarClassicExt", jetBetaStarClassicExt, &b_jetBetaStarClassicExt);
   fChain->SetBranchAddress("jetNNeutrals", jetNNeutrals, &b_jetNNeutrals);
   fChain->SetBranchAddress("jetNCharged", jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetMVAs", jetMVAs, &b_jetMVAs);
   fChain->SetBranchAddress("jetWPLevels", jetWPLevels, &b_jetWPLevels);
   fChain->SetBranchAddress("jetMVAsExt", jetMVAsExt, &b_jetMVAsExt);
   fChain->SetBranchAddress("jetWPLevelsExt", jetWPLevelsExt, &b_jetWPLevelsExt);
   fChain->SetBranchAddress("jetPartonID", jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("nLowPtJet", &nLowPtJet, &b_nLowPtJet);
   fChain->SetBranchAddress("jetLowPtEn", jetLowPtEn, &b_jetLowPtEn);
   fChain->SetBranchAddress("jetLowPtPt", jetLowPtPt, &b_jetLowPtPt);
   fChain->SetBranchAddress("jetLowPtEta", jetLowPtEta, &b_jetLowPtEta);
   fChain->SetBranchAddress("jetLowPtPhi", jetLowPtPhi, &b_jetLowPtPhi);
   fChain->SetBranchAddress("jetLowPtCharge", jetLowPtCharge, &b_jetLowPtCharge);
   fChain->SetBranchAddress("jetLowPtEt", jetLowPtEt, &b_jetLowPtEt);
   fChain->SetBranchAddress("jetLowPtRawPt", jetLowPtRawPt, &b_jetLowPtRawPt);
   fChain->SetBranchAddress("jetLowPtRawEn", jetLowPtRawEn, &b_jetLowPtRawEn);
   fChain->SetBranchAddress("jetLowPtArea", jetLowPtArea, &b_jetLowPtArea);
   fChain->SetBranchAddress("jetLowPtPartonID", jetLowPtPartonID, &b_jetLowPtPartonID);
   fChain->SetBranchAddress("jetLowPtGenJetEn", jetLowPtGenJetEn, &b_jetLowPtGenJetEn);
   fChain->SetBranchAddress("jetLowPtGenJetPt", jetLowPtGenJetPt, &b_jetLowPtGenJetPt);
   fChain->SetBranchAddress("jetLowPtGenJetEta", jetLowPtGenJetEta, &b_jetLowPtGenJetEta);
   fChain->SetBranchAddress("jetLowPtGenJetPhi", jetLowPtGenJetPhi, &b_jetLowPtGenJetPhi);
   fChain->SetBranchAddress("jetLowPtGenPartonID", jetLowPtGenPartonID, &b_jetLowPtGenPartonID);
   fChain->SetBranchAddress("jetLowPtGenEn", jetLowPtGenEn, &b_jetLowPtGenEn);
   fChain->SetBranchAddress("jetLowPtGenPt", jetLowPtGenPt, &b_jetLowPtGenPt);
   fChain->SetBranchAddress("jetLowPtGenEta", jetLowPtGenEta, &b_jetLowPtGenEta);
   fChain->SetBranchAddress("jetLowPtGenPhi", jetLowPtGenPhi, &b_jetLowPtGenPhi);
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("convP4", convP4, &b_convP4);
   fChain->SetBranchAddress("convVtx", convVtx, &b_convVtx);
   fChain->SetBranchAddress("convVtxErr", convVtxErr, &b_convVtxErr);
   fChain->SetBranchAddress("convPairMomentum", convPairMomentum, &b_convPairMomentum);
   fChain->SetBranchAddress("convRefittedMomentum", convRefittedMomentum, &b_convRefittedMomentum);
   fChain->SetBranchAddress("convNTracks", convNTracks, &b_convNTracks);
   fChain->SetBranchAddress("convPairInvMass", convPairInvMass, &b_convPairInvMass);
   fChain->SetBranchAddress("convPairCotThetaSep", convPairCotThetaSep, &b_convPairCotThetaSep);
   fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convDistOfMinApproach", convDistOfMinApproach, &b_convDistOfMinApproach);
   fChain->SetBranchAddress("convDPhiTrksAtVtx", convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
   fChain->SetBranchAddress("convDPhiTrksAtEcal", convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
   fChain->SetBranchAddress("convDEtaTrksAtEcal", convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
   fChain->SetBranchAddress("convDxy", convDxy, &b_convDxy);
   fChain->SetBranchAddress("convDz", convDz, &b_convDz);
   fChain->SetBranchAddress("convLxy", convLxy, &b_convLxy);
   fChain->SetBranchAddress("convLz", convLz, &b_convLz);
   fChain->SetBranchAddress("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
   fChain->SetBranchAddress("convNHitsBeforeVtx", convNHitsBeforeVtx, &b_convNHitsBeforeVtx);
   fChain->SetBranchAddress("convNSharedHits", convNSharedHits, &b_convNSharedHits);
   fChain->SetBranchAddress("convValidVtx", convValidVtx, &b_convValidVtx);
   fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convChi2", convChi2, &b_convChi2);
   fChain->SetBranchAddress("convChi2Probability", convChi2Probability, &b_convChi2Probability);
   fChain->SetBranchAddress("convTk1Dz", convTk1Dz, &b_convTk1Dz);
   fChain->SetBranchAddress("convTk2Dz", convTk2Dz, &b_convTk2Dz);
   fChain->SetBranchAddress("convTk1DzErr", convTk1DzErr, &b_convTk1DzErr);
   fChain->SetBranchAddress("convTk2DzErr", convTk2DzErr, &b_convTk2DzErr);
   fChain->SetBranchAddress("convCh1Ch2", convCh1Ch2, &b_convCh1Ch2);
   fChain->SetBranchAddress("convTk1D0", convTk1D0, &b_convTk1D0);
   fChain->SetBranchAddress("convTk1Pout", convTk1Pout, &b_convTk1Pout);
   fChain->SetBranchAddress("convTk1Pin", convTk1Pin, &b_convTk1Pin);
   fChain->SetBranchAddress("convTk2D0", convTk2D0, &b_convTk2D0);
   fChain->SetBranchAddress("convTk2Pout", convTk2Pout, &b_convTk2Pout);
   fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin);
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
