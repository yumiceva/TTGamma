//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 19 14:24:38 2013 by ROOT version 5.32/00
// from TTree EventTree/Event data
// found on file: ggtree_mc.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *fChainSkim;
   void            InitSkim();
   void            SetDirectory( TFile *fSkimFile);
   void            FillSkim( TFile *fSkimFile);
   void            AutoSave();
   void            Write();
   TTree*          GetSkim() { return fChainSkim;}
   Int_t           fCurrent; //!current Tree number in a TChain

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
   vector<float>   *vtx_x;
   vector<float>   *vtx_y;
   vector<float>   *vtx_z;
   Int_t           IsVtxGood;
   Int_t           nGoodVtx;
   Int_t           nVtxBS;
   vector<float>   *vtxbs_x;
   vector<float>   *vtxbs_y;
   vector<float>   *vtxbs_z;
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx_x;
   vector<float>   *mcVtx_y;
   vector<float>   *mcVtx_z;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<int>     *mcIndex;
   vector<int>     *mcDecayType;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfType01MET;
   Float_t         pfType01METPhi;
   Float_t         pfType01METsumEt;
   Float_t         pfType01METmEtSig;
   Float_t         pfType01METSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Float_t         trkMETxPV;
   Float_t         trkMETyPV;
   Float_t         trkMETPhiPV;
   Float_t         trkMETPV;
   vector<float>   *trkMETx;
   vector<float>   *trkMETy;
   vector<float>   *trkMETPhi;
   vector<float>   *trkMET;
   Int_t           metFilters[10];
   Int_t           nEle;
   vector<unsigned long> *eleTrg;
   vector<int>     *eleClass;
   vector<int>     *eleIsEcalDriven;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleEcalEn;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleESEn;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<vector<float> > *eleEtaVtx;
   vector<vector<float> > *elePhiVtx;
   vector<vector<float> > *eleEtVtx;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleVtx_x;
   vector<float>   *eleVtx_y;
   vector<float>   *eleVtx_z;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleD0GV;
   vector<float>   *eleDzGV;
   vector<vector<float> > *eleD0Vtx;
   vector<vector<float> > *eleDzVtx;
   vector<float>   *eleHoverE;
   vector<float>   *eleHoverE12;
   vector<float>   *eleEoverP;
   vector<float>   *elePin;
   vector<float>   *elePout;
   vector<float>   *eleTrkMomErr;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eleSigmaIEtaIEta;
   vector<float>   *eleSigmaIEtaIPhi;
   vector<float>   *eleSigmaIPhiIPhi;
   vector<float>   *eleEmax;
   vector<float>   *eleE2ndMax;
   vector<float>   *eleETop;
   vector<float>   *eleEBottom;
   vector<float>   *eleELeft;
   vector<float>   *eleERight;
   vector<float>   *eleE1x5;
   vector<float>   *eleE3x3;
   vector<float>   *eleE5x5;
   vector<float>   *eleE2x5Max;
   vector<float>   *eleE2x5Top;
   vector<float>   *eleE2x5Bottom;
   vector<float>   *eleE2x5Left;
   vector<float>   *eleE2x5Right;
   vector<float>   *eleSeedEta;
   vector<float>   *eleSeedE;
   vector<float>   *eleSeedPhi;
   vector<float>   *eleCrysEta;
   vector<float>   *eleCrysPhi;
   vector<int>     *eleCrysIEta;
   vector<int>     *eleCrysIPhi;
   vector<float>   *eleRegrE;
   vector<float>   *eleRegrEerr;
   vector<float>   *elePhoRegrE;
   vector<float>   *elePhoRegrEerr;
   vector<float>   *eleSeedTime;
  //vector<vector<float> > *eleGSFPt;
  //vector<vector<float> > *eleGSFEta;
  //vector<vector<float> > *eleGSFPhi;
  //vector<vector<float> > *eleGSFCharge;
  //vector<vector<float> > *eleGSFChi2NDF;
   vector<int>     *eleRecoFlag;
   vector<int>     *elePos;
   vector<int>     *eleGenIndex;
   vector<int>     *eleGenGMomPID;
   vector<int>     *eleGenMomPID;
   vector<float>   *eleGenMomPt;
   vector<float>   *eleIsoTrkDR03;
   vector<float>   *eleIsoEcalDR03;
   vector<float>   *eleIsoHcalDR03;
   vector<float>   *eleIsoHcalDR0312;
   vector<float>   *eleIsoTrkDR04;
   vector<float>   *eleIsoEcalDR04;
   vector<float>   *eleIsoHcalDR04;
   vector<float>   *eleIsoHcalDR0412;
   vector<float>   *eleModIsoTrk;
   vector<float>   *eleModIsoEcal;
   vector<float>   *eleModIsoHcal;
   vector<int>     *eleMissHits;
   vector<float>   *eleConvDist;
   vector<float>   *eleConvDcot;
   vector<int>     *eleConvVtxFit;
   vector<float>   *eleIP3D;
   vector<float>   *eleIP3DErr;
   vector<float>   *eleIDMVANonTrig;
   vector<float>   *eleIDMVATrig;
   vector<float>   *elePFChIso03;
   vector<float>   *elePFPhoIso03;
   vector<float>   *elePFNeuIso03;
   vector<float>   *elePFChIso04;
   vector<float>   *elePFPhoIso04;
   vector<float>   *elePFNeuIso04;
   vector<float>   *eleESEffSigmaRR_x;
   vector<float>   *eleESEffSigmaRR_y;
   vector<float>   *eleESEffSigmaRR_z;
   Int_t           nPho;
   vector<unsigned long> *phoTrg;
   vector<unsigned long> *phoTrgFilter;
   vector<bool>    *phoIsPhoton;
   vector<float>   *phoSCPos_x;
   vector<float>   *phoSCPos_y;
   vector<float>   *phoSCPos_z;
   vector<float>   *phoCaloPos_x;
   vector<float>   *phoCaloPos_y;
   vector<float>   *phoCaloPos_z;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoVtx_x;
   vector<float>   *phoVtx_y;
   vector<float>   *phoVtx_z;
   vector<float>   *phoPhi;
   vector<vector<float> > *phoEtVtx;
   vector<vector<float> > *phoEtaVtx;
   vector<vector<float> > *phoPhiVtx;
   vector<float>   *phoR9;
   vector<int>     *phoNClus;
   vector<float>   *phoTrkIsoHollowDR03;
   vector<float>   *phoEcalIsoDR03;
   vector<float>   *phoHcalIsoDR03;
   vector<float>   *phoHcalIsoDR0312;
   vector<float>   *phoTrkIsoHollowDR04;
   vector<float>   *phoCiCdRtoTrk;
   vector<float>   *phoEcalIsoDR04;
   vector<float>   *phoHcalIsoDR04;
   vector<float>   *phoHcalIsoDR0412;
   vector<float>   *phoHoverE;
   vector<float>   *phoHoverE12;
   vector<int>     *phoEleVeto;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoSigmaIEtaIPhi;
   vector<float>   *phoSigmaIPhiIPhi;
   vector<float>   *phoCiCPF4phopfIso03;
   vector<float>   *phoCiCPF4phopfIso04;
   vector<vector<float> > *phoCiCPF4chgpfIso02;
   vector<vector<float> > *phoCiCPF4chgpfIso03;
   vector<vector<float> > *phoCiCPF4chgpfIso04;
   vector<float>   *phoEmax;
   vector<float>   *phoETop;
   vector<float>   *phoEBottom;
   vector<float>   *phoELeft;
   vector<float>   *phoERight;
   vector<float>   *phoE2ndMax;
   vector<float>   *phoE3x3;
   vector<float>   *phoE3x1;
   vector<float>   *phoE1x3;
   vector<float>   *phoE5x5;
   vector<float>   *phoE1x5;
   vector<float>   *phoE2x2;
   vector<float>   *phoE2x5Max;
   vector<float>   *phoE2x5Top;
   vector<float>   *phoE2x5Bottom;
   vector<float>   *phoE2x5Left;
   vector<float>   *phoE2x5Right;
   vector<float>   *phoSeedE;
   vector<float>   *phoSeedEta;
   vector<float>   *phoSeedPhi;
   vector<float>   *phoCrysEta;
   vector<float>   *phoCrysPhi;
   vector<int>     *phoCrysIEta;
   vector<int>     *phoCrysIPhi;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoSCRChIso;
   vector<float>   *phoSCRPhoIso;
   vector<float>   *phoSCRNeuIso;
   vector<float>   *phoSCRChIso04;
   vector<float>   *phoSCRPhoIso04;
   vector<float>   *phoSCRNeuIso04;
   vector<float>   *phoRandConeChIso;
   vector<float>   *phoRandConePhoIso;
   vector<float>   *phoRandConeNeuIso;
   vector<float>   *phoRandConeChIso04;
   vector<float>   *phoRandConePhoIso04;
   vector<float>   *phoRandConeNeuIso04;
   vector<float>   *phoRegrE;
   vector<float>   *phoRegrEerr;
   vector<float>   *phoSeedTime;
   vector<int>     *phoSeedDetId1;
   vector<int>     *phoSeedDetId2;
   vector<float>   *phoLICTD;
   vector<int>     *phoRecoFlag;
   vector<int>     *phoPos;
   vector<int>     *phoGenIndex;
   vector<int>     *phoGenGMomPID;
   vector<int>     *phoGenMomPID;
   vector<float>   *phoGenMomPt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEn;
   vector<float>   *phoSCEt;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phoOverlap;
   vector<int>     *phohasPixelSeed;
   vector<int>     *pho_hasConvPf;
   vector<int>     *pho_hasSLConvPf;
   vector<float>   *pho_pfconvVtxZ;
   vector<float>   *pho_pfconvVtxZErr;
   vector<int>     *pho_nSLConv;
   vector<vector<float> > *pho_pfSLConvPos_x;
   vector<vector<float> > *pho_pfSLConvPos_y;
   vector<vector<float> > *pho_pfSLConvPos_z;
   vector<vector<float> > *pho_pfSLConvVtxZ;
   vector<int>     *phoIsConv;
   vector<int>     *phoNConv;
   vector<float>   *phoConvInvMass;
   vector<float>   *phoConvCotTheta;
   vector<float>   *phoConvEoverP;
   vector<float>   *phoConvZofPVfromTrks;
   vector<float>   *phoConvMinDist;
   vector<float>   *phoConvdPhiAtVtx;
   vector<float>   *phoConvdPhiAtCalo;
   vector<float>   *phoConvdEtaAtCalo;
   vector<float>   *phoConvTrkd0_x;
   vector<float>   *phoConvTrkd0_y;
   vector<float>   *phoConvTrkPin_x;
   vector<float>   *phoConvTrkPin_y;
   vector<float>   *phoConvTrkPout_x;
   vector<float>   *phoConvTrkPout_y;
   vector<float>   *phoConvTrkdz_x;
   vector<float>   *phoConvTrkdz_y;
   vector<float>   *phoConvTrkdzErr_x;
   vector<float>   *phoConvTrkdzErr_y;
   vector<float>   *phoConvChi2;
   vector<float>   *phoConvChi2Prob;
   vector<int>     *phoConvNTrks;
   vector<float>   *phoConvCharge1;
   vector<float>   *phoConvCharge2;
   vector<int>     *phoConvValidVtx;
   vector<float>   *phoConvLikeLihood;
   vector<float>   *phoConvP4_0;
   vector<float>   *phoConvP4_1;
   vector<float>   *phoConvP4_2;
   vector<float>   *phoConvP4_3;
   vector<float>   *phoConvVtx_x;
   vector<float>   *phoConvVtx_y;
   vector<float>   *phoConvVtx_z;
   vector<float>   *phoConvVtxErr_x;
   vector<float>   *phoConvVtxErr_y;
   vector<float>   *phoConvVtxErr_z;
   vector<float>   *phoConvPairMomentum_x;
   vector<float>   *phoConvPairMomentum_y;
   vector<float>   *phoConvPairMomentum_z;
   vector<float>   *phoConvRefittedMomentum_x;
   vector<float>   *phoConvRefittedMomentum_y;
   vector<float>   *phoConvRefittedMomentum_z;
   vector<int>     *SingleLegConv;
   vector<vector<float> > *phoPFConvVtx_x;
   vector<vector<float> > *phoPFConvVtx_y;
   vector<vector<float> > *phoPFConvVtx_z;
   vector<vector<float> > *phoPFConvMom_x;
   vector<vector<float> > *phoPFConvMom_y;
   vector<vector<float> > *phoPFConvMom_z;
   vector<float>   *phoESEffSigmaRR_x;
   vector<float>   *phoESEffSigmaRR_y;
   vector<float>   *phoESEffSigmaRR_z;
   Int_t           nMu;
   vector<unsigned long> *muTrg;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<float>   *muPt;
   vector<float>   *muPz;
   vector<float>   *muVtx_x;
   vector<float>   *muVtx_y;
   vector<float>   *muVtx_z;
   vector<float>   *muVtxGlb_x;
   vector<float>   *muVtxGlb_y;
   vector<float>   *muVtxGlb_z;
   vector<int>     *muGenIndex;
   vector<float>   *mucktPt;
   vector<float>   *mucktPtErr;
   vector<float>   *mucktEta;
   vector<float>   *mucktPhi;
   vector<float>   *mucktdxy;
   vector<float>   *mucktdz;
   vector<float>   *muIsoTrk;
   vector<float>   *muIsoCalo;
   vector<float>   *muIsoEcal;
   vector<float>   *muIsoHcal;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerChi2NDF;
   vector<float>   *muPFIsoR04_CH;
   vector<float>   *muPFIsoR04_NH;
   vector<float>   *muPFIsoR04_Pho;
   vector<float>   *muPFIsoR04_PU;
   vector<float>   *muPFIsoR04_CPart;
   vector<float>   *muPFIsoR04_NHHT;
   vector<float>   *muPFIsoR04_PhoHT;
   vector<float>   *muPFIsoR03_CH;
   vector<float>   *muPFIsoR03_NH;
   vector<float>   *muPFIsoR03_Pho;
   vector<float>   *muPFIsoR03_PU;
   vector<float>   *muPFIsoR03_CPart;
   vector<float>   *muPFIsoR03_NHHT;
   vector<float>   *muPFIsoR03_PhoHT;
   vector<int>     *muType;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muD0GV;
   vector<float>   *muDzGV;
   vector<vector<float> > *muD0Vtx;
   vector<vector<float> > *muDzVtx;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<float>   *muInnerD0GV;
   vector<float>   *muInnerDzGV;
   vector<float>   *muInnerPt;
   vector<float>   *muInnerPtErr;
   vector<int>     *muNumberOfValidTrkLayers;
   vector<int>     *muNumberOfValidTrkHits;
   vector<int>     *muNumberOfValidPixelLayers;
   vector<int>     *muNumberOfValidPixelHits;
   vector<int>     *muNumberOfValidMuonHits;
   vector<int>     *muStations;
   vector<int>     *muChambers;
   vector<float>   *muIP3D;
   vector<float>   *muIP3DErr;
   Int_t           nTau;
   vector<bool>    *tauDecayModeFinding;
   vector<bool>    *tauAgainstElectronLooseMVA3;
   vector<bool>    *tauAgainstElectronMediumMVA3;
   vector<bool>    *tauAgainstElectronTightMVA3;
   vector<bool>    *tauAgainstElectronVTightMVA3;
   vector<bool>    *tauAgainstElectronDeadECAL;
   vector<bool>    *tauAgainstMuonLoose2;
   vector<bool>    *tauAgainstMuonMedium2;
   vector<bool>    *tauAgainstMuonTight2;
   vector<bool>    *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<bool>    *tauLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *tauMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *tauTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<float>   *tauEta;
   vector<float>   *tauPhi;
   vector<float>   *tauPt;
   vector<float>   *tauEt;
   vector<float>   *tauCharge;
   vector<int>     *tauDecayMode;
   vector<float>   *tauEMFraction;
   vector<float>   *tauHCAL3x3OverPLead;
   vector<float>   *tauHCALMaxOverPLead;
   vector<float>   *tauHCALTotOverPLead;
   vector<float>   *tauIsolationPFChargedHadrCandsPtSum;
   vector<float>   *tauIsolationPFGammaCandsEtSum;
   vector<float>   *tauLeadPFChargedHadrCandsignedSipt;
   vector<bool>    *tauLeadChargedHadronExists;
   vector<float>   *tauLeadChargedHadronEta;
   vector<float>   *tauLeadChargedHadronPhi;
   vector<float>   *tauLeadChargedHadronPt;
   Float_t         rho25;
   Float_t         rho25_neu;
   Float_t         rho25_muPFiso;
   Float_t         rho25_elePFiso;
   Float_t         rho2011;
   Float_t         rho2012;
   Float_t         QGTag_MLP;
   Float_t         QGTag_likelihood;
   Int_t           nCA8Jet;
   vector<float>   *CA8JetPt;
   vector<float>   *CA8JetEta;
   vector<float>   *CA8JetPhi;
   vector<float>   *CA8JetMass;
   vector<float>   *CA8JetArea;
   vector<float>   *CA8Jet_tau1;
   vector<float>   *CA8Jet_tau2;
   vector<float>   *CA8Jet_tau3;
   vector<float>   *CA8JetCHF;
   vector<float>   *CA8JetNHF;
   vector<float>   *CA8JetCEF;
   vector<float>   *CA8JetNEF;
   vector<int>     *CA8JetNCH;
   vector<int>     *CA8Jetnconstituents;
   vector<float>   *CA8prunedJetMass;
   vector<int>     *CA8prunedJet_nSubJets;
   vector<vector<float> > *CA8prunedJet_SubjetPt;
   vector<vector<float> > *CA8prunedJet_SubjetEta;
   vector<vector<float> > *CA8prunedJet_SubjetPhi;
   vector<vector<float> > *CA8prunedJet_SubjetMass;
   Int_t           nJet;
   vector<unsigned long> *jetTrg;
   vector<float>   *jetEn;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetCharge;
   vector<float>   *jetEt;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetArea;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<int>     *jetNCH;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConstituents;
   vector<float>   *jetCombinedSecondaryVtxBJetTags;
   vector<float>   *jetCombinedSecondaryVtxMVABJetTags;
   vector<float>   *jetJetProbabilityBJetTags;
   vector<float>   *jetJetBProbabilityBJetTags;
   vector<vector<float> > *jetBetaStar;
   vector<bool>    *jetPFLooseId;
   vector<float>   *jetDRMean;
   vector<float>   *jetDR2Mean;
   vector<float>   *jetDZ;
   vector<float>   *jetFrac01;
   vector<float>   *jetFrac02;
   vector<float>   *jetFrac03;
   vector<float>   *jetFrac04;
   vector<float>   *jetFrac05;
   vector<float>   *jetFrac06;
   vector<float>   *jetFrac07;
   vector<float>   *jetBeta;
   vector<float>   *jetBetaStarCMG;
   vector<float>   *jetBetaStarClassic;
   vector<vector<float> > *jetBetaExt;
   vector<vector<float> > *jetBetaStarCMGExt;
   vector<vector<float> > *jetBetaStarClassicExt;
   vector<float>   *jetNNeutrals;
   vector<float>   *jetNCharged;
   vector<vector<float> > *jetMVAs;
   vector<vector<int> > *jetWPLevels;
   vector<vector<float> > *jetMVAsExt_simple;
   vector<vector<int> > *jetWPLevelsExt_simple;
   vector<vector<float> > *jetMVAsExt_full;
   vector<vector<int> > *jetWPLevelsExt_full;
   vector<vector<float> > *jetMVAsExt_cutBased;
   vector<vector<int> > *jetWPLevelsExt_cutBased;
   vector<vector<float> > *jetMVAsExt_philv1;
   vector<vector<int> > *jetWPLevelsExt_philv1;
   vector<float>   *jetMt;
   vector<float>   *jetJECUnc;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtx3dL;
   vector<float>   *jetVtx3deL;
   vector<float>   *jetSoftLeptPt;
   vector<float>   *jetSoftLeptPtRel;
   vector<float>   *jetSoftLeptdR;
   vector<float>   *jetSoftLeptIdlooseMu;
   vector<float>   *jetSoftLeptIdEle95;
   vector<float>   *jetDPhiMETJet;
   vector<float>   *jetPuJetIdL;
   vector<float>   *jetPuJetIdM;
   vector<float>   *jetPuJetIdT;
   vector<int>     *jetPartonID;
   vector<int>     *jetGenJetIndex;
   vector<float>   *jetGenJetEn;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenEn;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   Int_t           nConv;
   vector<float>   *convP4_x;
   vector<float>   *convP4_y;
   vector<float>   *convP4_z;
   vector<float>   *convP4_E;
   vector<float>   *convVtx_x;
   vector<float>   *convVtx_y;
   vector<float>   *convVtx_z;
   vector<float>   *convVtxErr_x;
   vector<float>   *convVtxErr_y;
   vector<float>   *convVtxErr_z;
   vector<float>   *convPairMomentum_x;
   vector<float>   *convPairMomentum_y;
   vector<float>   *convPairMomentum_z;
   vector<float>   *convRefittedMomentum_x;
   vector<float>   *convRefittedMomentum_y;
   vector<float>   *convRefittedMomentum_z;
   vector<int>     *convNTracks;
   vector<float>   *convPairInvMass;
   vector<float>   *convPairCotThetaSep;
   vector<float>   *convEoverP;
   vector<float>   *convDistOfMinApproach;
   vector<float>   *convDPhiTrksAtVtx;
   vector<float>   *convDPhiTrksAtEcal;
   vector<float>   *convDEtaTrksAtEcal;
   vector<float>   *convDxy;
   vector<float>   *convDz;
   vector<float>   *convLxy;
   vector<float>   *convLz;
   vector<float>   *convZofPrimVtxFromTrks;
   vector<int>     *convNHitsBeforeVtx_0;
   vector<int>     *convNHitsBeforeVtx_1;
   vector<int>     *convNSharedHits;
   vector<int>     *convValidVtx;
   vector<float>   *convMVALikelihood;
   vector<float>   *convChi2;
   vector<float>   *convChi2Probability;
   vector<float>   *convTk1Dz;
   vector<float>   *convTk2Dz;
   vector<float>   *convTk1DzErr;
   vector<float>   *convTk2DzErr;
   vector<int>     *convCh1Ch2;
   vector<float>   *convTk1D0;
   vector<float>   *convTk1Pout;
   vector<float>   *convTk1Pin;
   vector<float>   *convTk2D0;
   vector<float>   *convTk2Pout;
   vector<float>   *convTk2Pin;

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
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs_x;   //!
   TBranch        *b_vtxbs_y;   //!
   TBranch        *b_vtxbs_z;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
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
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfType01MET;   //!
   TBranch        *b_pfType01METPhi;   //!
   TBranch        *b_pfType01METsumEt;   //!
   TBranch        *b_pfType01METmEtSig;   //!
   TBranch        *b_pfType01METSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_trkMETxPV;   //!
   TBranch        *b_trkMETyPV;   //!
   TBranch        *b_trkMETPhiPV;   //!
   TBranch        *b_trkMETPV;   //!
   TBranch        *b_trkMETx;   //!
   TBranch        *b_trkMETy;   //!
   TBranch        *b_trkMETPhi;   //!
   TBranch        *b_trkMET;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleIsEcalDriven;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleEtVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx_x;   //!
   TBranch        *b_eleVtx_y;   //!
   TBranch        *b_eleVtx_z;   //!
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
   TBranch        *b_eleE2ndMax;   //!
   TBranch        *b_eleETop;   //!
   TBranch        *b_eleEBottom;   //!
   TBranch        *b_eleELeft;   //!
   TBranch        *b_eleERight;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Max;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleSeedEta;   //!
   TBranch        *b_eleSeedE;   //!
   TBranch        *b_eleSeedPhi;   //!
   TBranch        *b_eleCrysEta;   //!
   TBranch        *b_eleCrysPhi;   //!
   TBranch        *b_eleCrysIEta;   //!
   TBranch        *b_eleCrysIPhi;   //!
   TBranch        *b_eleRegrE;   //!
   TBranch        *b_eleRegrEerr;   //!
   TBranch        *b_elePhoRegrE;   //!
   TBranch        *b_elePhoRegrEerr;   //!
   TBranch        *b_eleSeedTime;   //!
  //TBranch        *b_eleGSFPt;   //!
  //TBranch        *b_eleGSFEta;   //!
  //TBranch        *b_eleGSFPhi;   //!
  //TBranch        *b_eleGSFCharge;   //!
  //TBranch        *b_eleGSFChi2NDF;   //!
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
   TBranch        *b_eleModIsoTrk;   //!
   TBranch        *b_eleModIsoEcal;   //!
   TBranch        *b_eleModIsoHcal;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleConvVtxFit;   //!
   TBranch        *b_eleIP3D;   //!
   TBranch        *b_eleIP3DErr;   //!
   TBranch        *b_eleIDMVANonTrig;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_elePFChIso04;   //!
   TBranch        *b_elePFPhoIso04;   //!
   TBranch        *b_elePFNeuIso04;   //!
   TBranch        *b_eleESEffSigmaRR_x;   //!
   TBranch        *b_eleESEffSigmaRR_y;   //!
   TBranch        *b_eleESEffSigmaRR_z;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoTrgFilter;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoSCPos_x;   //!
   TBranch        *b_phoSCPos_y;   //!
   TBranch        *b_phoSCPos_z;   //!
   TBranch        *b_phoCaloPos_x;   //!
   TBranch        *b_phoCaloPos_y;   //!
   TBranch        *b_phoCaloPos_z;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx_x;   //!
   TBranch        *b_phoVtx_y;   //!
   TBranch        *b_phoVtx_z;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoNClus;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR0312;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
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
   TBranch        *b_phoETop;   //!
   TBranch        *b_phoEBottom;   //!
   TBranch        *b_phoELeft;   //!
   TBranch        *b_phoERight;   //!
   TBranch        *b_phoE2ndMax;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE3x1;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoSeedE;   //!
   TBranch        *b_phoSeedEta;   //!
   TBranch        *b_phoSeedPhi;   //!
   TBranch        *b_phoCrysEta;   //!
   TBranch        *b_phoCrysPhi;   //!
   TBranch        *b_phoCrysIEta;   //!
   TBranch        *b_phoCrysIPhi;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoSCRChIso;   //!
   TBranch        *b_phoSCRPhoIso;   //!
   TBranch        *b_phoSCRNeuIso;   //!
   TBranch        *b_phoSCRChIso04;   //!
   TBranch        *b_phoSCRPhoIso04;   //!
   TBranch        *b_phoSCRNeuIso04;   //!
   TBranch        *b_phoRandConeChIso;   //!
   TBranch        *b_phoRandConePhoIso;   //!
   TBranch        *b_phoRandConeNeuIso;   //!
   TBranch        *b_phoRandConeChIso04;   //!
   TBranch        *b_phoRandConePhoIso04;   //!
   TBranch        *b_phoRandConeNeuIso04;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrEerr;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoLICTD;   //!
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
   TBranch        *b_pho_pfSLConvPos_x;   //!
   TBranch        *b_pho_pfSLConvPos_y;   //!
   TBranch        *b_pho_pfSLConvPos_z;   //!
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
   TBranch        *b_phoConvTrkd0_x;   //!
   TBranch        *b_phoConvTrkd0_y;   //!
   TBranch        *b_phoConvTrkPin_x;   //!
   TBranch        *b_phoConvTrkPin_y;   //!
   TBranch        *b_phoConvTrkPout_x;   //!
   TBranch        *b_phoConvTrkPout_y;   //!
   TBranch        *b_phoConvTrkdz_x;   //!
   TBranch        *b_phoConvTrkdz_y;   //!
   TBranch        *b_phoConvTrkdzErr_x;   //!
   TBranch        *b_phoConvTrkdzErr_y;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvNTrks;   //!
   TBranch        *b_phoConvCharge1;   //!
   TBranch        *b_phoConvCharge2;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4_0;   //!
   TBranch        *b_phoConvP4_1;   //!
   TBranch        *b_phoConvP4_2;   //!
   TBranch        *b_phoConvP4_3;   //!
   TBranch        *b_phoConvVtx_x;   //!
   TBranch        *b_phoConvVtx_y;   //!
   TBranch        *b_phoConvVtx_z;   //!
   TBranch        *b_phoConvVtxErr_x;   //!
   TBranch        *b_phoConvVtxErr_y;   //!
   TBranch        *b_phoConvVtxErr_z;   //!
   TBranch        *b_phoConvPairMomentum_x;   //!
   TBranch        *b_phoConvPairMomentum_y;   //!
   TBranch        *b_phoConvPairMomentum_z;   //!
   TBranch        *b_phoConvRefittedMomentum_x;   //!
   TBranch        *b_phoConvRefittedMomentum_y;   //!
   TBranch        *b_phoConvRefittedMomentum_z;   //!
   TBranch        *b_SingleLegConv;   //!
   TBranch        *b_phoPFConvVtx_x;   //!
   TBranch        *b_phoPFConvVtx_y;   //!
   TBranch        *b_phoPFConvVtx_z;   //!
   TBranch        *b_phoPFConvMom_x;   //!
   TBranch        *b_phoPFConvMom_y;   //!
   TBranch        *b_phoPFConvMom_z;   //!
   TBranch        *b_phoESEffSigmaRR_x;   //!
   TBranch        *b_phoESEffSigmaRR_y;   //!
   TBranch        *b_phoESEffSigmaRR_z;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muVtx_x;   //!
   TBranch        *b_muVtx_y;   //!
   TBranch        *b_muVtx_z;   //!
   TBranch        *b_muVtxGlb_x;   //!
   TBranch        *b_muVtxGlb_y;   //!
   TBranch        *b_muVtxGlb_z;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_mucktPt;   //!
   TBranch        *b_mucktPtErr;   //!
   TBranch        *b_mucktEta;   //!
   TBranch        *b_mucktPhi;   //!
   TBranch        *b_mucktdxy;   //!
   TBranch        *b_mucktdz;   //!
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
   TBranch        *b_muInnerPt;   //!
   TBranch        *b_muInnerPtErr;   //!
   TBranch        *b_muNumberOfValidTrkLayers;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelLayers;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_muIP3D;   //!
   TBranch        *b_muIP3DErr;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_tauDecayModeFinding;   //!
   TBranch        *b_tauAgainstElectronLooseMVA3;   //!
   TBranch        *b_tauAgainstElectronMediumMVA3;   //!
   TBranch        *b_tauAgainstElectronTightMVA3;   //!
   TBranch        *b_tauAgainstElectronVTightMVA3;   //!
   TBranch        *b_tauAgainstElectronDeadECAL;   //!
   TBranch        *b_tauAgainstMuonLoose2;   //!
   TBranch        *b_tauAgainstMuonMedium2;   //!
   TBranch        *b_tauAgainstMuonTight2;   //!
   TBranch        *b_tauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tauLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEt;   //!
   TBranch        *b_tauCharge;   //!
   TBranch        *b_tauDecayMode;   //!
   TBranch        *b_tauEMFraction;   //!
   TBranch        *b_tauHCAL3x3OverPLead;   //!
   TBranch        *b_tauHCALMaxOverPLead;   //!
   TBranch        *b_tauHCALTotOverPLead;   //!
   TBranch        *b_tauIsolationPFChargedHadrCandsPtSum;   //!
   TBranch        *b_tauIsolationPFGammaCandsEtSum;   //!
   TBranch        *b_tauLeadPFChargedHadrCandsignedSipt;   //!
   TBranch        *b_tauLeadChargedHadronExists;   //!
   TBranch        *b_tauLeadChargedHadronEta;   //!
   TBranch        *b_tauLeadChargedHadronPhi;   //!
   TBranch        *b_tauLeadChargedHadronPt;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho25_neu;   //!
   TBranch        *b_rho25_muPFiso;   //!
   TBranch        *b_rho25_elePFiso;   //!
   TBranch        *b_rho2011;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_QGTag_MLP;   //!
   TBranch        *b_QGTag_likelihood;   //!
   TBranch        *b_nCA8Jet;   //!
   TBranch        *b_CA8JetPt;   //!
   TBranch        *b_CA8JetEta;   //!
   TBranch        *b_CA8JetPhi;   //!
   TBranch        *b_CA8JetMass;   //!
   TBranch        *b_CA8JetArea;   //!
   TBranch        *b_CA8Jet_tau1;   //!
   TBranch        *b_CA8Jet_tau2;   //!
   TBranch        *b_CA8Jet_tau3;   //!
   TBranch        *b_CA8JetCHF;   //!
   TBranch        *b_CA8JetNHF;   //!
   TBranch        *b_CA8JetCEF;   //!
   TBranch        *b_CA8JetNEF;   //!
   TBranch        *b_CA8JetNCH;   //!
   TBranch        *b_CA8Jetnconstituents;   //!
   TBranch        *b_CA8prunedJetMass;   //!
   TBranch        *b_CA8prunedJet_nSubJets;   //!
   TBranch        *b_CA8prunedJet_SubjetPt;   //!
   TBranch        *b_CA8prunedJet_SubjetEta;   //!
   TBranch        *b_CA8prunedJet_SubjetPhi;   //!
   TBranch        *b_CA8prunedJet_SubjetMass;   //!
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
   TBranch        *b_jetMVAsExt_simple;   //!
   TBranch        *b_jetWPLevelsExt_simple;   //!
   TBranch        *b_jetMVAsExt_full;   //!
   TBranch        *b_jetWPLevelsExt_full;   //!
   TBranch        *b_jetMVAsExt_cutBased;   //!
   TBranch        *b_jetWPLevelsExt_cutBased;   //!
   TBranch        *b_jetMVAsExt_philv1;   //!
   TBranch        *b_jetWPLevelsExt_philv1;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtx3dL;   //!
   TBranch        *b_jetVtx3deL;   //!
   TBranch        *b_jetSoftLeptPt;   //!
   TBranch        *b_jetSoftLeptPtRel;   //!
   TBranch        *b_jetSoftLeptdR;   //!
   TBranch        *b_jetSoftLeptIdlooseMu;   //!
   TBranch        *b_jetSoftLeptIdEle95;   //!
   TBranch        *b_jetDPhiMETJet;   //!
   TBranch        *b_jetPuJetIdL;   //!
   TBranch        *b_jetPuJetIdM;   //!
   TBranch        *b_jetPuJetIdT;   //!
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
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4_x;   //!
   TBranch        *b_convP4_y;   //!
   TBranch        *b_convP4_z;   //!
   TBranch        *b_convP4_E;   //!
   TBranch        *b_convVtx_x;   //!
   TBranch        *b_convVtx_y;   //!
   TBranch        *b_convVtx_z;   //!
   TBranch        *b_convVtxErr_x;   //!
   TBranch        *b_convVtxErr_y;   //!
   TBranch        *b_convVtxErr_z;   //!
   TBranch        *b_convPairMomentum_x;   //!
   TBranch        *b_convPairMomentum_y;   //!
   TBranch        *b_convPairMomentum_z;   //!
   TBranch        *b_convRefittedMomentum_x;   //!
   TBranch        *b_convRefittedMomentum_y;   //!
   TBranch        *b_convRefittedMomentum_z;   //!
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
   TBranch        *b_convNHitsBeforeVtx_0;   //!
   TBranch        *b_convNHitsBeforeVtx_1;   //!
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

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree) : fChain(0), fChainSkim(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ggtree_mc.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ggtree_mc.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ggtree_mc.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);   
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtxbs_x = 0;
   vtxbs_y = 0;
   vtxbs_z = 0;
   mcPID = 0;
   mcVtx_x = 0;
   mcVtx_y = 0;
   mcVtx_z = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcDecayType = 0;
   mcParentage = 0;
   mcStatus = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   trkMETx = 0;
   trkMETy = 0;
   trkMETPhi = 0;
   trkMET = 0;
   eleTrg = 0;
   eleClass = 0;
   eleIsEcalDriven = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleEcalEn = 0;
   eleSCRawEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleEtaVtx = 0;
   elePhiVtx = 0;
   eleEtVtx = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleVtx_x = 0;
   eleVtx_y = 0;
   eleVtx_z = 0;
   eleD0 = 0;
   eleDz = 0;
   eleD0GV = 0;
   eleDzGV = 0;
   eleD0Vtx = 0;
   eleDzVtx = 0;
   eleHoverE = 0;
   eleHoverE12 = 0;
   eleEoverP = 0;
   elePin = 0;
   elePout = 0;
   eleTrkMomErr = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIPhi = 0;
   eleSigmaIPhiIPhi = 0;
   eleEmax = 0;
   eleE2ndMax = 0;
   eleETop = 0;
   eleEBottom = 0;
   eleELeft = 0;
   eleERight = 0;
   eleE1x5 = 0;
   eleE3x3 = 0;
   eleE5x5 = 0;
   eleE2x5Max = 0;
   eleE2x5Top = 0;
   eleE2x5Bottom = 0;
   eleE2x5Left = 0;
   eleE2x5Right = 0;
   eleSeedEta = 0;
   eleSeedE = 0;
   eleSeedPhi = 0;
   eleCrysEta = 0;
   eleCrysPhi = 0;
   eleCrysIEta = 0;
   eleCrysIPhi = 0;
   eleRegrE = 0;
   eleRegrEerr = 0;
   elePhoRegrE = 0;
   elePhoRegrEerr = 0;
   eleSeedTime = 0;
   //eleGSFPt = 0;
   //eleGSFEta = 0;
   //eleGSFPhi = 0;
   //eleGSFCharge = 0;
   //eleGSFChi2NDF = 0;
   eleRecoFlag = 0;
   elePos = 0;
   eleGenIndex = 0;
   eleGenGMomPID = 0;
   eleGenMomPID = 0;
   eleGenMomPt = 0;
   eleIsoTrkDR03 = 0;
   eleIsoEcalDR03 = 0;
   eleIsoHcalDR03 = 0;
   eleIsoHcalDR0312 = 0;
   eleIsoTrkDR04 = 0;
   eleIsoEcalDR04 = 0;
   eleIsoHcalDR04 = 0;
   eleIsoHcalDR0412 = 0;
   eleModIsoTrk = 0;
   eleModIsoEcal = 0;
   eleModIsoHcal = 0;
   eleMissHits = 0;
   eleConvDist = 0;
   eleConvDcot = 0;
   eleConvVtxFit = 0;
   eleIP3D = 0;
   eleIP3DErr = 0;
   eleIDMVANonTrig = 0;
   eleIDMVATrig = 0;
   elePFChIso03 = 0;
   elePFPhoIso03 = 0;
   elePFNeuIso03 = 0;
   elePFChIso04 = 0;
   elePFPhoIso04 = 0;
   elePFNeuIso04 = 0;
   eleESEffSigmaRR_x = 0;
   eleESEffSigmaRR_y = 0;
   eleESEffSigmaRR_z = 0;
   phoTrg = 0;
   phoTrgFilter = 0;
   phoIsPhoton = 0;
   phoSCPos_x = 0;
   phoSCPos_y = 0;
   phoSCPos_z = 0;
   phoCaloPos_x = 0;
   phoCaloPos_y = 0;
   phoCaloPos_z = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoVtx_x = 0;
   phoVtx_y = 0;
   phoVtx_z = 0;
   phoPhi = 0;
   phoEtVtx = 0;
   phoEtaVtx = 0;
   phoPhiVtx = 0;
   phoR9 = 0;
   phoNClus = 0;
   phoTrkIsoHollowDR03 = 0;
   phoEcalIsoDR03 = 0;
   phoHcalIsoDR03 = 0;
   phoHcalIsoDR0312 = 0;
   phoTrkIsoHollowDR04 = 0;
   phoCiCdRtoTrk = 0;
   phoEcalIsoDR04 = 0;
   phoHcalIsoDR04 = 0;
   phoHcalIsoDR0412 = 0;
   phoHoverE = 0;
   phoHoverE12 = 0;
   phoEleVeto = 0;
   phoSigmaIEtaIEta = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoCiCPF4phopfIso03 = 0;
   phoCiCPF4phopfIso04 = 0;
   phoCiCPF4chgpfIso02 = 0;
   phoCiCPF4chgpfIso03 = 0;
   phoCiCPF4chgpfIso04 = 0;
   phoEmax = 0;
   phoETop = 0;
   phoEBottom = 0;
   phoELeft = 0;
   phoERight = 0;
   phoE2ndMax = 0;
   phoE3x3 = 0;
   phoE3x1 = 0;
   phoE1x3 = 0;
   phoE5x5 = 0;
   phoE1x5 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE2x5Top = 0;
   phoE2x5Bottom = 0;
   phoE2x5Left = 0;
   phoE2x5Right = 0;
   phoSeedE = 0;
   phoSeedEta = 0;
   phoSeedPhi = 0;
   phoCrysEta = 0;
   phoCrysPhi = 0;
   phoCrysIEta = 0;
   phoCrysIPhi = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoSCRChIso = 0;
   phoSCRPhoIso = 0;
   phoSCRNeuIso = 0;
   phoSCRChIso04 = 0;
   phoSCRPhoIso04 = 0;
   phoSCRNeuIso04 = 0;
   phoRandConeChIso = 0;
   phoRandConePhoIso = 0;
   phoRandConeNeuIso = 0;
   phoRandConeChIso04 = 0;
   phoRandConePhoIso04 = 0;
   phoRandConeNeuIso04 = 0;
   phoRegrE = 0;
   phoRegrEerr = 0;
   phoSeedTime = 0;
   phoSeedDetId1 = 0;
   phoSeedDetId2 = 0;
   phoLICTD = 0;
   phoRecoFlag = 0;
   phoPos = 0;
   phoGenIndex = 0;
   phoGenGMomPID = 0;
   phoGenMomPID = 0;
   phoGenMomPt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoSCEt = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phoOverlap = 0;
   phohasPixelSeed = 0;
   pho_hasConvPf = 0;
   pho_hasSLConvPf = 0;
   pho_pfconvVtxZ = 0;
   pho_pfconvVtxZErr = 0;
   pho_nSLConv = 0;
   pho_pfSLConvPos_x = 0;
   pho_pfSLConvPos_y = 0;
   pho_pfSLConvPos_z = 0;
   pho_pfSLConvVtxZ = 0;
   phoIsConv = 0;
   phoNConv = 0;
   phoConvInvMass = 0;
   phoConvCotTheta = 0;
   phoConvEoverP = 0;
   phoConvZofPVfromTrks = 0;
   phoConvMinDist = 0;
   phoConvdPhiAtVtx = 0;
   phoConvdPhiAtCalo = 0;
   phoConvdEtaAtCalo = 0;
   phoConvTrkd0_x = 0;
   phoConvTrkd0_y = 0;
   phoConvTrkPin_x = 0;
   phoConvTrkPin_y = 0;
   phoConvTrkPout_x = 0;
   phoConvTrkPout_y = 0;
   phoConvTrkdz_x = 0;
   phoConvTrkdz_y = 0;
   phoConvTrkdzErr_x = 0;
   phoConvTrkdzErr_y = 0;
   phoConvChi2 = 0;
   phoConvChi2Prob = 0;
   phoConvNTrks = 0;
   phoConvCharge1 = 0;
   phoConvCharge2 = 0;
   phoConvValidVtx = 0;
   phoConvLikeLihood = 0;
   phoConvP4_0 = 0;
   phoConvP4_1 = 0;
   phoConvP4_2 = 0;
   phoConvP4_3 = 0;
   phoConvVtx_x = 0;
   phoConvVtx_y = 0;
   phoConvVtx_z = 0;
   phoConvVtxErr_x = 0;
   phoConvVtxErr_y = 0;
   phoConvVtxErr_z = 0;
   phoConvPairMomentum_x = 0;
   phoConvPairMomentum_y = 0;
   phoConvPairMomentum_z = 0;
   phoConvRefittedMomentum_x = 0;
   phoConvRefittedMomentum_y = 0;
   phoConvRefittedMomentum_z = 0;
   SingleLegConv = 0;
   phoPFConvVtx_x = 0;
   phoPFConvVtx_y = 0;
   phoPFConvVtx_z = 0;
   phoPFConvMom_x = 0;
   phoPFConvMom_y = 0;
   phoPFConvMom_z = 0;
   phoESEffSigmaRR_x = 0;
   phoESEffSigmaRR_y = 0;
   phoESEffSigmaRR_z = 0;
   muTrg = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muPt = 0;
   muPz = 0;
   muVtx_x = 0;
   muVtx_y = 0;
   muVtx_z = 0;
   muVtxGlb_x = 0;
   muVtxGlb_y = 0;
   muVtxGlb_z = 0;
   muGenIndex = 0;
   mucktPt = 0;
   mucktPtErr = 0;
   mucktEta = 0;
   mucktPhi = 0;
   mucktdxy = 0;
   mucktdz = 0;
   muIsoTrk = 0;
   muIsoCalo = 0;
   muIsoEcal = 0;
   muIsoHcal = 0;
   muChi2NDF = 0;
   muInnerChi2NDF = 0;
   muPFIsoR04_CH = 0;
   muPFIsoR04_NH = 0;
   muPFIsoR04_Pho = 0;
   muPFIsoR04_PU = 0;
   muPFIsoR04_CPart = 0;
   muPFIsoR04_NHHT = 0;
   muPFIsoR04_PhoHT = 0;
   muPFIsoR03_CH = 0;
   muPFIsoR03_NH = 0;
   muPFIsoR03_Pho = 0;
   muPFIsoR03_PU = 0;
   muPFIsoR03_CPart = 0;
   muPFIsoR03_NHHT = 0;
   muPFIsoR03_PhoHT = 0;
   muType = 0;
   muD0 = 0;
   muDz = 0;
   muD0GV = 0;
   muDzGV = 0;
   muD0Vtx = 0;
   muDzVtx = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muInnerD0GV = 0;
   muInnerDzGV = 0;
   muInnerPt = 0;
   muInnerPtErr = 0;
   muNumberOfValidTrkLayers = 0;
   muNumberOfValidTrkHits = 0;
   muNumberOfValidPixelLayers = 0;
   muNumberOfValidPixelHits = 0;
   muNumberOfValidMuonHits = 0;
   muStations = 0;
   muChambers = 0;
   muIP3D = 0;
   muIP3DErr = 0;
   tauDecayModeFinding = 0;
   tauAgainstElectronLooseMVA3 = 0;
   tauAgainstElectronMediumMVA3 = 0;
   tauAgainstElectronTightMVA3 = 0;
   tauAgainstElectronVTightMVA3 = 0;
   tauAgainstElectronDeadECAL = 0;
   tauAgainstMuonLoose2 = 0;
   tauAgainstMuonMedium2 = 0;
   tauAgainstMuonTight2 = 0;
   tauCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   tauLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauEta = 0;
   tauPhi = 0;
   tauPt = 0;
   tauEt = 0;
   tauCharge = 0;
   tauDecayMode = 0;
   tauEMFraction = 0;
   tauHCAL3x3OverPLead = 0;
   tauHCALMaxOverPLead = 0;
   tauHCALTotOverPLead = 0;
   tauIsolationPFChargedHadrCandsPtSum = 0;
   tauIsolationPFGammaCandsEtSum = 0;
   tauLeadPFChargedHadrCandsignedSipt = 0;
   tauLeadChargedHadronExists = 0;
   tauLeadChargedHadronEta = 0;
   tauLeadChargedHadronPhi = 0;
   tauLeadChargedHadronPt = 0;
   CA8JetPt = 0;
   CA8JetEta = 0;
   CA8JetPhi = 0;
   CA8JetMass = 0;
   CA8JetArea = 0;
   CA8Jet_tau1 = 0;
   CA8Jet_tau2 = 0;
   CA8Jet_tau3 = 0;
   CA8JetCHF = 0;
   CA8JetNHF = 0;
   CA8JetCEF = 0;
   CA8JetNEF = 0;
   CA8JetNCH = 0;
   CA8Jetnconstituents = 0;
   CA8prunedJetMass = 0;
   CA8prunedJet_nSubJets = 0;
   CA8prunedJet_SubjetPt = 0;
   CA8prunedJet_SubjetEta = 0;
   CA8prunedJet_SubjetPhi = 0;
   CA8prunedJet_SubjetMass = 0;
   jetTrg = 0;
   jetEn = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetCharge = 0;
   jetEt = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetArea = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   jetCombinedSecondaryVtxBJetTags = 0;
   jetCombinedSecondaryVtxMVABJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetJetBProbabilityBJetTags = 0;
   jetBetaStar = 0;
   jetPFLooseId = 0;
   jetDRMean = 0;
   jetDR2Mean = 0;
   jetDZ = 0;
   jetFrac01 = 0;
   jetFrac02 = 0;
   jetFrac03 = 0;
   jetFrac04 = 0;
   jetFrac05 = 0;
   jetFrac06 = 0;
   jetFrac07 = 0;
   jetBeta = 0;
   jetBetaStarCMG = 0;
   jetBetaStarClassic = 0;
   jetBetaExt = 0;
   jetBetaStarCMGExt = 0;
   jetBetaStarClassicExt = 0;
   jetNNeutrals = 0;
   jetNCharged = 0;
   jetMVAs = 0;
   jetWPLevels = 0;
   jetMVAsExt_simple = 0;
   jetWPLevelsExt_simple = 0;
   jetMVAsExt_full = 0;
   jetWPLevelsExt_full = 0;
   jetMVAsExt_cutBased = 0;
   jetWPLevelsExt_cutBased = 0;
   jetMVAsExt_philv1 = 0;
   jetWPLevelsExt_philv1 = 0;
   jetMt = 0;
   jetJECUnc = 0;
   jetLeadTrackPt = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtx3dL = 0;
   jetVtx3deL = 0;
   jetSoftLeptPt = 0;
   jetSoftLeptPtRel = 0;
   jetSoftLeptdR = 0;
   jetSoftLeptIdlooseMu = 0;
   jetSoftLeptIdEle95 = 0;
   jetDPhiMETJet = 0;
   jetPuJetIdL = 0;
   jetPuJetIdM = 0;
   jetPuJetIdT = 0;
   jetPartonID = 0;
   jetGenJetIndex = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   convP4_x = 0;
   convP4_y = 0;
   convP4_z = 0;
   convP4_E = 0;
   convVtx_x = 0;
   convVtx_y = 0;
   convVtx_z = 0;
   convVtxErr_x = 0;
   convVtxErr_y = 0;
   convVtxErr_z = 0;
   convPairMomentum_x = 0;
   convPairMomentum_y = 0;
   convPairMomentum_z = 0;
   convRefittedMomentum_x = 0;
   convRefittedMomentum_y = 0;
   convRefittedMomentum_z = 0;
   convNTracks = 0;
   convPairInvMass = 0;
   convPairCotThetaSep = 0;
   convEoverP = 0;
   convDistOfMinApproach = 0;
   convDPhiTrksAtVtx = 0;
   convDPhiTrksAtEcal = 0;
   convDEtaTrksAtEcal = 0;
   convDxy = 0;
   convDz = 0;
   convLxy = 0;
   convLz = 0;
   convZofPrimVtxFromTrks = 0;
   convNHitsBeforeVtx_0 = 0;
   convNHitsBeforeVtx_1 = 0;
   convNSharedHits = 0;
   convValidVtx = 0;
   convMVALikelihood = 0;
   convChi2 = 0;
   convChi2Probability = 0;
   convTk1Dz = 0;
   convTk2Dz = 0;
   convTk1DzErr = 0;
   convTk2DzErr = 0;
   convCh1Ch2 = 0;
   convTk1D0 = 0;
   convTk1Pout = 0;
   convTk1Pin = 0;
   convTk2D0 = 0;
   convTk2Pout = 0;
   convTk2Pin = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
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
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs_x", &vtxbs_x, &b_vtxbs_x);
   fChain->SetBranchAddress("vtxbs_y", &vtxbs_y, &b_vtxbs_y);
   fChain->SetBranchAddress("vtxbs_z", &vtxbs_z, &b_vtxbs_z);
   fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx_x", &mcVtx_x, &b_mcVtx_x);
   fChain->SetBranchAddress("mcVtx_y", &mcVtx_y, &b_mcVtx_y);
   fChain->SetBranchAddress("mcVtx_z", &mcVtx_z, &b_mcVtx_z);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcDecayType", &mcDecayType, &b_mcDecayType);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfType01MET", &pfType01MET, &b_pfType01MET);
   fChain->SetBranchAddress("pfType01METPhi", &pfType01METPhi, &b_pfType01METPhi);
   fChain->SetBranchAddress("pfType01METsumEt", &pfType01METsumEt, &b_pfType01METsumEt);
   fChain->SetBranchAddress("pfType01METmEtSig", &pfType01METmEtSig, &b_pfType01METmEtSig);
   fChain->SetBranchAddress("pfType01METSig", &pfType01METSig, &b_pfType01METSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("trkMETxPV", &trkMETxPV, &b_trkMETxPV);
   fChain->SetBranchAddress("trkMETyPV", &trkMETyPV, &b_trkMETyPV);
   fChain->SetBranchAddress("trkMETPhiPV", &trkMETPhiPV, &b_trkMETPhiPV);
   fChain->SetBranchAddress("trkMETPV", &trkMETPV, &b_trkMETPV);
   fChain->SetBranchAddress("trkMETx", &trkMETx, &b_trkMETx);
   fChain->SetBranchAddress("trkMETy", &trkMETy, &b_trkMETy);
   fChain->SetBranchAddress("trkMETPhi", &trkMETPhi, &b_trkMETPhi);
   fChain->SetBranchAddress("trkMET", &trkMET, &b_trkMET);
   fChain->SetBranchAddress("metFilters", metFilters, &b_metFilters);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", &eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleClass", &eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleIsEcalDriven", &eleIsEcalDriven, &b_eleIsEcalDriven);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleEtaVtx", &eleEtaVtx, &b_eleEtaVtx);
   fChain->SetBranchAddress("elePhiVtx", &elePhiVtx, &b_elePhiVtx);
   fChain->SetBranchAddress("eleEtVtx", &eleEtVtx, &b_eleEtVtx);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx_x", &eleVtx_x, &b_eleVtx_x);
   fChain->SetBranchAddress("eleVtx_y", &eleVtx_y, &b_eleVtx_y);
   fChain->SetBranchAddress("eleVtx_z", &eleVtx_z, &b_eleVtx_z);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleD0GV", &eleD0GV, &b_eleD0GV);
   fChain->SetBranchAddress("eleDzGV", &eleDzGV, &b_eleDzGV);
   fChain->SetBranchAddress("eleD0Vtx", &eleD0Vtx, &b_eleD0Vtx);
   fChain->SetBranchAddress("eleDzVtx", &eleDzVtx, &b_eleDzVtx);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleHoverE12", &eleHoverE12, &b_eleHoverE12);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", &elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", &elePout, &b_elePout);
   fChain->SetBranchAddress("eleTrkMomErr", &eleTrkMomErr, &b_eleTrkMomErr);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleEmax", &eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE2ndMax", &eleE2ndMax, &b_eleE2ndMax);
   fChain->SetBranchAddress("eleETop", &eleETop, &b_eleETop);
   fChain->SetBranchAddress("eleEBottom", &eleEBottom, &b_eleEBottom);
   fChain->SetBranchAddress("eleELeft", &eleELeft, &b_eleELeft);
   fChain->SetBranchAddress("eleERight", &eleERight, &b_eleERight);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE3x3", &eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Max", &eleE2x5Max, &b_eleE2x5Max);
   fChain->SetBranchAddress("eleE2x5Top", &eleE2x5Top, &b_eleE2x5Top);
   fChain->SetBranchAddress("eleE2x5Bottom", &eleE2x5Bottom, &b_eleE2x5Bottom);
   fChain->SetBranchAddress("eleE2x5Left", &eleE2x5Left, &b_eleE2x5Left);
   fChain->SetBranchAddress("eleE2x5Right", &eleE2x5Right, &b_eleE2x5Right);
   fChain->SetBranchAddress("eleSeedEta", &eleSeedEta, &b_eleSeedEta);
   fChain->SetBranchAddress("eleSeedE", &eleSeedE, &b_eleSeedE);
   fChain->SetBranchAddress("eleSeedPhi", &eleSeedPhi, &b_eleSeedPhi);
   fChain->SetBranchAddress("eleCrysEta", &eleCrysEta, &b_eleCrysEta);
   fChain->SetBranchAddress("eleCrysPhi", &eleCrysPhi, &b_eleCrysPhi);
   fChain->SetBranchAddress("eleCrysIEta", &eleCrysIEta, &b_eleCrysIEta);
   fChain->SetBranchAddress("eleCrysIPhi", &eleCrysIPhi, &b_eleCrysIPhi);
   fChain->SetBranchAddress("eleRegrE", &eleRegrE, &b_eleRegrE);
   fChain->SetBranchAddress("eleRegrEerr", &eleRegrEerr, &b_eleRegrEerr);
   fChain->SetBranchAddress("elePhoRegrE", &elePhoRegrE, &b_elePhoRegrE);
   fChain->SetBranchAddress("elePhoRegrEerr", &elePhoRegrEerr, &b_elePhoRegrEerr);
   fChain->SetBranchAddress("eleSeedTime", &eleSeedTime, &b_eleSeedTime);
   //fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   //fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   //fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   //fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   //fChain->SetBranchAddress("eleGSFChi2NDF", &eleGSFChi2NDF, &b_eleGSFChi2NDF);
   fChain->SetBranchAddress("eleRecoFlag", &eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("elePos", &elePos, &b_elePos);
   fChain->SetBranchAddress("eleGenIndex", &eleGenIndex, &b_eleGenIndex);
   fChain->SetBranchAddress("eleGenGMomPID", &eleGenGMomPID, &b_eleGenGMomPID);
   fChain->SetBranchAddress("eleGenMomPID", &eleGenMomPID, &b_eleGenMomPID);
   fChain->SetBranchAddress("eleGenMomPt", &eleGenMomPt, &b_eleGenMomPt);
   fChain->SetBranchAddress("eleIsoTrkDR03", &eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", &eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", &eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR0312", &eleIsoHcalDR0312, &b_eleIsoHcalDR0312);
   fChain->SetBranchAddress("eleIsoTrkDR04", &eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", &eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", &eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR0412", &eleIsoHcalDR0412, &b_eleIsoHcalDR0412);
   fChain->SetBranchAddress("eleModIsoTrk", &eleModIsoTrk, &b_eleModIsoTrk);
   fChain->SetBranchAddress("eleModIsoEcal", &eleModIsoEcal, &b_eleModIsoEcal);
   fChain->SetBranchAddress("eleModIsoHcal", &eleModIsoHcal, &b_eleModIsoHcal);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", &eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", &eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleConvVtxFit", &eleConvVtxFit, &b_eleConvVtxFit);
   fChain->SetBranchAddress("eleIP3D", &eleIP3D, &b_eleIP3D);
   fChain->SetBranchAddress("eleIP3DErr", &eleIP3DErr, &b_eleIP3DErr);
   fChain->SetBranchAddress("eleIDMVANonTrig", &eleIDMVANonTrig, &b_eleIDMVANonTrig);
   fChain->SetBranchAddress("eleIDMVATrig", &eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("elePFChIso03", &elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("elePFChIso04", &elePFChIso04, &b_elePFChIso04);
   fChain->SetBranchAddress("elePFPhoIso04", &elePFPhoIso04, &b_elePFPhoIso04);
   fChain->SetBranchAddress("elePFNeuIso04", &elePFNeuIso04, &b_elePFNeuIso04);
   fChain->SetBranchAddress("eleESEffSigmaRR_x", &eleESEffSigmaRR_x, &b_eleESEffSigmaRR_x);
   fChain->SetBranchAddress("eleESEffSigmaRR_y", &eleESEffSigmaRR_y, &b_eleESEffSigmaRR_y);
   fChain->SetBranchAddress("eleESEffSigmaRR_z", &eleESEffSigmaRR_z, &b_eleESEffSigmaRR_z);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", &phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoTrgFilter", &phoTrgFilter, &b_phoTrgFilter);
   fChain->SetBranchAddress("phoIsPhoton", &phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoSCPos_x", &phoSCPos_x, &b_phoSCPos_x);
   fChain->SetBranchAddress("phoSCPos_y", &phoSCPos_y, &b_phoSCPos_y);
   fChain->SetBranchAddress("phoSCPos_z", &phoSCPos_z, &b_phoSCPos_z);
   fChain->SetBranchAddress("phoCaloPos_x", &phoCaloPos_x, &b_phoCaloPos_x);
   fChain->SetBranchAddress("phoCaloPos_y", &phoCaloPos_y, &b_phoCaloPos_y);
   fChain->SetBranchAddress("phoCaloPos_z", &phoCaloPos_z, &b_phoCaloPos_z);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx_x", &phoVtx_x, &b_phoVtx_x);
   fChain->SetBranchAddress("phoVtx_y", &phoVtx_y, &b_phoVtx_y);
   fChain->SetBranchAddress("phoVtx_z", &phoVtx_z, &b_phoVtx_z);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", &phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", &phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", &phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoNClus", &phoNClus, &b_phoNClus);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", &phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", &phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", &phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR0312", &phoHcalIsoDR0312, &b_phoHcalIsoDR0312);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", &phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoCiCdRtoTrk", &phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", &phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", &phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR0412", &phoHcalIsoDR0412, &b_phoHcalIsoDR0412);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoHoverE12", &phoHoverE12, &b_phoHoverE12);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoCiCPF4phopfIso03", &phoCiCPF4phopfIso03, &b_phoCiCPF4phopfIso03);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04", &phoCiCPF4phopfIso04, &b_phoCiCPF4phopfIso04);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso02", &phoCiCPF4chgpfIso02, &b_phoCiCPF4chgpfIso02);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso03", &phoCiCPF4chgpfIso03, &b_phoCiCPF4chgpfIso03);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04", &phoCiCPF4chgpfIso04, &b_phoCiCPF4chgpfIso04);
   fChain->SetBranchAddress("phoEmax", &phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoETop", &phoETop, &b_phoETop);
   fChain->SetBranchAddress("phoEBottom", &phoEBottom, &b_phoEBottom);
   fChain->SetBranchAddress("phoELeft", &phoELeft, &b_phoELeft);
   fChain->SetBranchAddress("phoERight", &phoERight, &b_phoERight);
   fChain->SetBranchAddress("phoE2ndMax", &phoE2ndMax, &b_phoE2ndMax);
   fChain->SetBranchAddress("phoE3x3", &phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE3x1", &phoE3x1, &b_phoE3x1);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE2x5Top", &phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", &phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoE2x5Left", &phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Right", &phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoSeedE", &phoSeedE, &b_phoSeedE);
   fChain->SetBranchAddress("phoSeedEta", &phoSeedEta, &b_phoSeedEta);
   fChain->SetBranchAddress("phoSeedPhi", &phoSeedPhi, &b_phoSeedPhi);
   fChain->SetBranchAddress("phoCrysEta", &phoCrysEta, &b_phoCrysEta);
   fChain->SetBranchAddress("phoCrysPhi", &phoCrysPhi, &b_phoCrysPhi);
   fChain->SetBranchAddress("phoCrysIEta", &phoCrysIEta, &b_phoCrysIEta);
   fChain->SetBranchAddress("phoCrysIPhi", &phoCrysIPhi, &b_phoCrysIPhi);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoSCRChIso", &phoSCRChIso, &b_phoSCRChIso);
   fChain->SetBranchAddress("phoSCRPhoIso", &phoSCRPhoIso, &b_phoSCRPhoIso);
   fChain->SetBranchAddress("phoSCRNeuIso", &phoSCRNeuIso, &b_phoSCRNeuIso);
   fChain->SetBranchAddress("phoSCRChIso04", &phoSCRChIso04, &b_phoSCRChIso04);
   fChain->SetBranchAddress("phoSCRPhoIso04", &phoSCRPhoIso04, &b_phoSCRPhoIso04);
   fChain->SetBranchAddress("phoSCRNeuIso04", &phoSCRNeuIso04, &b_phoSCRNeuIso04);
   fChain->SetBranchAddress("phoRandConeChIso", &phoRandConeChIso, &b_phoRandConeChIso);
   fChain->SetBranchAddress("phoRandConePhoIso", &phoRandConePhoIso, &b_phoRandConePhoIso);
   fChain->SetBranchAddress("phoRandConeNeuIso", &phoRandConeNeuIso, &b_phoRandConeNeuIso);
   fChain->SetBranchAddress("phoRandConeChIso04", &phoRandConeChIso04, &b_phoRandConeChIso04);
   fChain->SetBranchAddress("phoRandConePhoIso04", &phoRandConePhoIso04, &b_phoRandConePhoIso04);
   fChain->SetBranchAddress("phoRandConeNeuIso04", &phoRandConeNeuIso04, &b_phoRandConeNeuIso04);
   fChain->SetBranchAddress("phoRegrE", &phoRegrE, &b_phoRegrE);
   fChain->SetBranchAddress("phoRegrEerr", &phoRegrEerr, &b_phoRegrEerr);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", &phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", &phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoLICTD", &phoLICTD, &b_phoLICTD);
   fChain->SetBranchAddress("phoRecoFlag", &phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", &phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoGenIndex", &phoGenIndex, &b_phoGenIndex);
   fChain->SetBranchAddress("phoGenGMomPID", &phoGenGMomPID, &b_phoGenGMomPID);
   fChain->SetBranchAddress("phoGenMomPID", &phoGenMomPID, &b_phoGenMomPID);
   fChain->SetBranchAddress("phoGenMomPt", &phoGenMomPt, &b_phoGenMomPt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoSCEt", &phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", &phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("pho_hasConvPf", &pho_hasConvPf, &b_pho_hasConvPf);
   fChain->SetBranchAddress("pho_hasSLConvPf", &pho_hasSLConvPf, &b_pho_hasSLConvPf);
   fChain->SetBranchAddress("pho_pfconvVtxZ", &pho_pfconvVtxZ, &b_pho_pfconvVtxZ);
   fChain->SetBranchAddress("pho_pfconvVtxZErr", &pho_pfconvVtxZErr, &b_pho_pfconvVtxZErr);
   fChain->SetBranchAddress("pho_nSLConv", &pho_nSLConv, &b_pho_nSLConv);
   fChain->SetBranchAddress("pho_pfSLConvPos_x", &pho_pfSLConvPos_x, &b_pho_pfSLConvPos_x);
   fChain->SetBranchAddress("pho_pfSLConvPos_y", &pho_pfSLConvPos_y, &b_pho_pfSLConvPos_y);
   fChain->SetBranchAddress("pho_pfSLConvPos_z", &pho_pfSLConvPos_z, &b_pho_pfSLConvPos_z);
   fChain->SetBranchAddress("pho_pfSLConvVtxZ", &pho_pfSLConvVtxZ, &b_pho_pfSLConvVtxZ);
   fChain->SetBranchAddress("phoIsConv", &phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", &phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", &phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", &phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", &phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", &phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", &phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", &phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", &phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", &phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0_x", &phoConvTrkd0_x, &b_phoConvTrkd0_x);
   fChain->SetBranchAddress("phoConvTrkd0_y", &phoConvTrkd0_y, &b_phoConvTrkd0_y);
   fChain->SetBranchAddress("phoConvTrkPin_x", &phoConvTrkPin_x, &b_phoConvTrkPin_x);
   fChain->SetBranchAddress("phoConvTrkPin_y", &phoConvTrkPin_y, &b_phoConvTrkPin_y);
   fChain->SetBranchAddress("phoConvTrkPout_x", &phoConvTrkPout_x, &b_phoConvTrkPout_x);
   fChain->SetBranchAddress("phoConvTrkPout_y", &phoConvTrkPout_y, &b_phoConvTrkPout_y);
   fChain->SetBranchAddress("phoConvTrkdz_x", &phoConvTrkdz_x, &b_phoConvTrkdz_x);
   fChain->SetBranchAddress("phoConvTrkdz_y", &phoConvTrkdz_y, &b_phoConvTrkdz_y);
   fChain->SetBranchAddress("phoConvTrkdzErr_x", &phoConvTrkdzErr_x, &b_phoConvTrkdzErr_x);
   fChain->SetBranchAddress("phoConvTrkdzErr_y", &phoConvTrkdzErr_y, &b_phoConvTrkdzErr_y);
   fChain->SetBranchAddress("phoConvChi2", &phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", &phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvNTrks", &phoConvNTrks, &b_phoConvNTrks);
   fChain->SetBranchAddress("phoConvCharge1", &phoConvCharge1, &b_phoConvCharge1);
   fChain->SetBranchAddress("phoConvCharge2", &phoConvCharge2, &b_phoConvCharge2);
   fChain->SetBranchAddress("phoConvValidVtx", &phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", &phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4_0", &phoConvP4_0, &b_phoConvP4_0);
   fChain->SetBranchAddress("phoConvP4_1", &phoConvP4_1, &b_phoConvP4_1);
   fChain->SetBranchAddress("phoConvP4_2", &phoConvP4_2, &b_phoConvP4_2);
   fChain->SetBranchAddress("phoConvP4_3", &phoConvP4_3, &b_phoConvP4_3);
   fChain->SetBranchAddress("phoConvVtx_x", &phoConvVtx_x, &b_phoConvVtx_x);
   fChain->SetBranchAddress("phoConvVtx_y", &phoConvVtx_y, &b_phoConvVtx_y);
   fChain->SetBranchAddress("phoConvVtx_z", &phoConvVtx_z, &b_phoConvVtx_z);
   fChain->SetBranchAddress("phoConvVtxErr_x", &phoConvVtxErr_x, &b_phoConvVtxErr_x);
   fChain->SetBranchAddress("phoConvVtxErr_y", &phoConvVtxErr_y, &b_phoConvVtxErr_y);
   fChain->SetBranchAddress("phoConvVtxErr_z", &phoConvVtxErr_z, &b_phoConvVtxErr_z);
   fChain->SetBranchAddress("phoConvPairMomentum_x", &phoConvPairMomentum_x, &b_phoConvPairMomentum_x);
   fChain->SetBranchAddress("phoConvPairMomentum_y", &phoConvPairMomentum_y, &b_phoConvPairMomentum_y);
   fChain->SetBranchAddress("phoConvPairMomentum_z", &phoConvPairMomentum_z, &b_phoConvPairMomentum_z);
   fChain->SetBranchAddress("phoConvRefittedMomentum_x", &phoConvRefittedMomentum_x, &b_phoConvRefittedMomentum_x);
   fChain->SetBranchAddress("phoConvRefittedMomentum_y", &phoConvRefittedMomentum_y, &b_phoConvRefittedMomentum_y);
   fChain->SetBranchAddress("phoConvRefittedMomentum_z", &phoConvRefittedMomentum_z, &b_phoConvRefittedMomentum_z);
   fChain->SetBranchAddress("SingleLegConv", &SingleLegConv, &b_SingleLegConv);
   fChain->SetBranchAddress("phoPFConvVtx_x", &phoPFConvVtx_x, &b_phoPFConvVtx_x);
   fChain->SetBranchAddress("phoPFConvVtx_y", &phoPFConvVtx_y, &b_phoPFConvVtx_y);
   fChain->SetBranchAddress("phoPFConvVtx_z", &phoPFConvVtx_z, &b_phoPFConvVtx_z);
   fChain->SetBranchAddress("phoPFConvMom_x", &phoPFConvMom_x, &b_phoPFConvMom_x);
   fChain->SetBranchAddress("phoPFConvMom_y", &phoPFConvMom_y, &b_phoPFConvMom_y);
   fChain->SetBranchAddress("phoPFConvMom_z", &phoPFConvMom_z, &b_phoPFConvMom_z);
   fChain->SetBranchAddress("phoESEffSigmaRR_x", &phoESEffSigmaRR_x, &b_phoESEffSigmaRR_x);
   fChain->SetBranchAddress("phoESEffSigmaRR_y", &phoESEffSigmaRR_y, &b_phoESEffSigmaRR_y);
   fChain->SetBranchAddress("phoESEffSigmaRR_z", &phoESEffSigmaRR_z, &b_phoESEffSigmaRR_z);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", &muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", &muPz, &b_muPz);
   fChain->SetBranchAddress("muVtx_x", &muVtx_x, &b_muVtx_x);
   fChain->SetBranchAddress("muVtx_y", &muVtx_y, &b_muVtx_y);
   fChain->SetBranchAddress("muVtx_z", &muVtx_z, &b_muVtx_z);
   fChain->SetBranchAddress("muVtxGlb_x", &muVtxGlb_x, &b_muVtxGlb_x);
   fChain->SetBranchAddress("muVtxGlb_y", &muVtxGlb_y, &b_muVtxGlb_y);
   fChain->SetBranchAddress("muVtxGlb_z", &muVtxGlb_z, &b_muVtxGlb_z);
   fChain->SetBranchAddress("muGenIndex", &muGenIndex, &b_muGenIndex);
   fChain->SetBranchAddress("mucktPt", &mucktPt, &b_mucktPt);
   fChain->SetBranchAddress("mucktPtErr", &mucktPtErr, &b_mucktPtErr);
   fChain->SetBranchAddress("mucktEta", &mucktEta, &b_mucktEta);
   fChain->SetBranchAddress("mucktPhi", &mucktPhi, &b_mucktPhi);
   fChain->SetBranchAddress("mucktdxy", &mucktdxy, &b_mucktdxy);
   fChain->SetBranchAddress("mucktdz", &mucktdz, &b_mucktdz);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", &muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", &muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", &muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerChi2NDF", &muInnerChi2NDF, &b_muInnerChi2NDF);
   fChain->SetBranchAddress("muPFIsoR04_CH", &muPFIsoR04_CH, &b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH", &muPFIsoR04_NH, &b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho", &muPFIsoR04_Pho, &b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("muPFIsoR04_PU", &muPFIsoR04_PU, &b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR04_CPart", &muPFIsoR04_CPart, &b_muPFIsoR04_CPart);
   fChain->SetBranchAddress("muPFIsoR04_NHHT", &muPFIsoR04_NHHT, &b_muPFIsoR04_NHHT);
   fChain->SetBranchAddress("muPFIsoR04_PhoHT", &muPFIsoR04_PhoHT, &b_muPFIsoR04_PhoHT);
   fChain->SetBranchAddress("muPFIsoR03_CH", &muPFIsoR03_CH, &b_muPFIsoR03_CH);
   fChain->SetBranchAddress("muPFIsoR03_NH", &muPFIsoR03_NH, &b_muPFIsoR03_NH);
   fChain->SetBranchAddress("muPFIsoR03_Pho", &muPFIsoR03_Pho, &b_muPFIsoR03_Pho);
   fChain->SetBranchAddress("muPFIsoR03_PU", &muPFIsoR03_PU, &b_muPFIsoR03_PU);
   fChain->SetBranchAddress("muPFIsoR03_CPart", &muPFIsoR03_CPart, &b_muPFIsoR03_CPart);
   fChain->SetBranchAddress("muPFIsoR03_NHHT", &muPFIsoR03_NHHT, &b_muPFIsoR03_NHHT);
   fChain->SetBranchAddress("muPFIsoR03_PhoHT", &muPFIsoR03_PhoHT, &b_muPFIsoR03_PhoHT);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muD0GV", &muD0GV, &b_muD0GV);
   fChain->SetBranchAddress("muDzGV", &muDzGV, &b_muDzGV);
   fChain->SetBranchAddress("muD0Vtx", &muD0Vtx, &b_muD0Vtx);
   fChain->SetBranchAddress("muDzVtx", &muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muInnerD0GV", &muInnerD0GV, &b_muInnerD0GV);
   fChain->SetBranchAddress("muInnerDzGV", &muInnerDzGV, &b_muInnerDzGV);
   fChain->SetBranchAddress("muInnerPt", &muInnerPt, &b_muInnerPt);
   fChain->SetBranchAddress("muInnerPtErr", &muInnerPtErr, &b_muInnerPtErr);
   fChain->SetBranchAddress("muNumberOfValidTrkLayers", &muNumberOfValidTrkLayers, &b_muNumberOfValidTrkLayers);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", &muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelLayers", &muNumberOfValidPixelLayers, &b_muNumberOfValidPixelLayers);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", &muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", &muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", &muChambers, &b_muChambers);
   fChain->SetBranchAddress("muIP3D", &muIP3D, &b_muIP3D);
   fChain->SetBranchAddress("muIP3DErr", &muIP3DErr, &b_muIP3DErr);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("tauDecayModeFinding", &tauDecayModeFinding, &b_tauDecayModeFinding);
   fChain->SetBranchAddress("tauAgainstElectronLooseMVA3", &tauAgainstElectronLooseMVA3, &b_tauAgainstElectronLooseMVA3);
   fChain->SetBranchAddress("tauAgainstElectronMediumMVA3", &tauAgainstElectronMediumMVA3, &b_tauAgainstElectronMediumMVA3);
   fChain->SetBranchAddress("tauAgainstElectronTightMVA3", &tauAgainstElectronTightMVA3, &b_tauAgainstElectronTightMVA3);
   fChain->SetBranchAddress("tauAgainstElectronVTightMVA3", &tauAgainstElectronVTightMVA3, &b_tauAgainstElectronVTightMVA3);
   fChain->SetBranchAddress("tauAgainstElectronDeadECAL", &tauAgainstElectronDeadECAL, &b_tauAgainstElectronDeadECAL);
   fChain->SetBranchAddress("tauAgainstMuonLoose2", &tauAgainstMuonLoose2, &b_tauAgainstMuonLoose2);
   fChain->SetBranchAddress("tauAgainstMuonMedium2", &tauAgainstMuonMedium2, &b_tauAgainstMuonMedium2);
   fChain->SetBranchAddress("tauAgainstMuonTight2", &tauAgainstMuonTight2, &b_tauAgainstMuonTight2);
   fChain->SetBranchAddress("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tauCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tauLooseCombinedIsolationDeltaBetaCorr3Hits", &tauLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tauLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauMediumCombinedIsolationDeltaBetaCorr3Hits", &tauMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tauMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauTightCombinedIsolationDeltaBetaCorr3Hits", &tauTightCombinedIsolationDeltaBetaCorr3Hits, &b_tauTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEt", &tauEt, &b_tauEt);
   fChain->SetBranchAddress("tauCharge", &tauCharge, &b_tauCharge);
   fChain->SetBranchAddress("tauDecayMode", &tauDecayMode, &b_tauDecayMode);
   fChain->SetBranchAddress("tauEMFraction", &tauEMFraction, &b_tauEMFraction);
   fChain->SetBranchAddress("tauHCAL3x3OverPLead", &tauHCAL3x3OverPLead, &b_tauHCAL3x3OverPLead);
   fChain->SetBranchAddress("tauHCALMaxOverPLead", &tauHCALMaxOverPLead, &b_tauHCALMaxOverPLead);
   fChain->SetBranchAddress("tauHCALTotOverPLead", &tauHCALTotOverPLead, &b_tauHCALTotOverPLead);
   fChain->SetBranchAddress("tauIsolationPFChargedHadrCandsPtSum", &tauIsolationPFChargedHadrCandsPtSum, &b_tauIsolationPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("tauIsolationPFGammaCandsEtSum", &tauIsolationPFGammaCandsEtSum, &b_tauIsolationPFGammaCandsEtSum);
   fChain->SetBranchAddress("tauLeadPFChargedHadrCandsignedSipt", &tauLeadPFChargedHadrCandsignedSipt, &b_tauLeadPFChargedHadrCandsignedSipt);
   fChain->SetBranchAddress("tauLeadChargedHadronExists", &tauLeadChargedHadronExists, &b_tauLeadChargedHadronExists);
   fChain->SetBranchAddress("tauLeadChargedHadronEta", &tauLeadChargedHadronEta, &b_tauLeadChargedHadronEta);
   fChain->SetBranchAddress("tauLeadChargedHadronPhi", &tauLeadChargedHadronPhi, &b_tauLeadChargedHadronPhi);
   fChain->SetBranchAddress("tauLeadChargedHadronPt", &tauLeadChargedHadronPt, &b_tauLeadChargedHadronPt);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho25_neu", &rho25_neu, &b_rho25_neu);
   fChain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso, &b_rho25_muPFiso);
   fChain->SetBranchAddress("rho25_elePFiso", &rho25_elePFiso, &b_rho25_elePFiso);
   fChain->SetBranchAddress("rho2011", &rho2011, &b_rho2011);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("QGTag_MLP", &QGTag_MLP, &b_QGTag_MLP);
   fChain->SetBranchAddress("QGTag_likelihood", &QGTag_likelihood, &b_QGTag_likelihood);
   fChain->SetBranchAddress("nCA8Jet", &nCA8Jet, &b_nCA8Jet);
   fChain->SetBranchAddress("CA8JetPt", &CA8JetPt, &b_CA8JetPt);
   fChain->SetBranchAddress("CA8JetEta", &CA8JetEta, &b_CA8JetEta);
   fChain->SetBranchAddress("CA8JetPhi", &CA8JetPhi, &b_CA8JetPhi);
   fChain->SetBranchAddress("CA8JetMass", &CA8JetMass, &b_CA8JetMass);
   fChain->SetBranchAddress("CA8JetArea", &CA8JetArea, &b_CA8JetArea);
   fChain->SetBranchAddress("CA8Jet_tau1", &CA8Jet_tau1, &b_CA8Jet_tau1);
   fChain->SetBranchAddress("CA8Jet_tau2", &CA8Jet_tau2, &b_CA8Jet_tau2);
   fChain->SetBranchAddress("CA8Jet_tau3", &CA8Jet_tau3, &b_CA8Jet_tau3);
   fChain->SetBranchAddress("CA8JetCHF", &CA8JetCHF, &b_CA8JetCHF);
   fChain->SetBranchAddress("CA8JetNHF", &CA8JetNHF, &b_CA8JetNHF);
   fChain->SetBranchAddress("CA8JetCEF", &CA8JetCEF, &b_CA8JetCEF);
   fChain->SetBranchAddress("CA8JetNEF", &CA8JetNEF, &b_CA8JetNEF);
   fChain->SetBranchAddress("CA8JetNCH", &CA8JetNCH, &b_CA8JetNCH);
   fChain->SetBranchAddress("CA8Jetnconstituents", &CA8Jetnconstituents, &b_CA8Jetnconstituents);
   fChain->SetBranchAddress("CA8prunedJetMass", &CA8prunedJetMass, &b_CA8prunedJetMass);
   fChain->SetBranchAddress("CA8prunedJet_nSubJets", &CA8prunedJet_nSubJets, &b_CA8prunedJet_nSubJets);
   fChain->SetBranchAddress("CA8prunedJet_SubjetPt", &CA8prunedJet_SubjetPt, &b_CA8prunedJet_SubjetPt);
   fChain->SetBranchAddress("CA8prunedJet_SubjetEta", &CA8prunedJet_SubjetEta, &b_CA8prunedJet_SubjetEta);
   fChain->SetBranchAddress("CA8prunedJet_SubjetPhi", &CA8prunedJet_SubjetPhi, &b_CA8prunedJet_SubjetPhi);
   fChain->SetBranchAddress("CA8prunedJet_SubjetMass", &CA8prunedJet_SubjetMass, &b_CA8prunedJet_SubjetMass);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", &jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCharge", &jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetEt", &jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", &jetCombinedSecondaryVtxMVABJetTags, &b_jetCombinedSecondaryVtxMVABJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags, &b_jetJetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetBetaStar", &jetBetaStar, &b_jetBetaStar);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetDRMean", &jetDRMean, &b_jetDRMean);
   fChain->SetBranchAddress("jetDR2Mean", &jetDR2Mean, &b_jetDR2Mean);
   fChain->SetBranchAddress("jetDZ", &jetDZ, &b_jetDZ);
   fChain->SetBranchAddress("jetFrac01", &jetFrac01, &b_jetFrac01);
   fChain->SetBranchAddress("jetFrac02", &jetFrac02, &b_jetFrac02);
   fChain->SetBranchAddress("jetFrac03", &jetFrac03, &b_jetFrac03);
   fChain->SetBranchAddress("jetFrac04", &jetFrac04, &b_jetFrac04);
   fChain->SetBranchAddress("jetFrac05", &jetFrac05, &b_jetFrac05);
   fChain->SetBranchAddress("jetFrac06", &jetFrac06, &b_jetFrac06);
   fChain->SetBranchAddress("jetFrac07", &jetFrac07, &b_jetFrac07);
   fChain->SetBranchAddress("jetBeta", &jetBeta, &b_jetBeta);
   fChain->SetBranchAddress("jetBetaStarCMG", &jetBetaStarCMG, &b_jetBetaStarCMG);
   fChain->SetBranchAddress("jetBetaStarClassic", &jetBetaStarClassic, &b_jetBetaStarClassic);
   fChain->SetBranchAddress("jetBetaExt", &jetBetaExt, &b_jetBetaExt);
   fChain->SetBranchAddress("jetBetaStarCMGExt", &jetBetaStarCMGExt, &b_jetBetaStarCMGExt);
   fChain->SetBranchAddress("jetBetaStarClassicExt", &jetBetaStarClassicExt, &b_jetBetaStarClassicExt);
   fChain->SetBranchAddress("jetNNeutrals", &jetNNeutrals, &b_jetNNeutrals);
   fChain->SetBranchAddress("jetNCharged", &jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetMVAs", &jetMVAs, &b_jetMVAs);
   fChain->SetBranchAddress("jetWPLevels", &jetWPLevels, &b_jetWPLevels);
   fChain->SetBranchAddress("jetMVAsExt_simple", &jetMVAsExt_simple, &b_jetMVAsExt_simple);
   fChain->SetBranchAddress("jetWPLevelsExt_simple", &jetWPLevelsExt_simple, &b_jetWPLevelsExt_simple);
   fChain->SetBranchAddress("jetMVAsExt_full", &jetMVAsExt_full, &b_jetMVAsExt_full);
   fChain->SetBranchAddress("jetWPLevelsExt_full", &jetWPLevelsExt_full, &b_jetWPLevelsExt_full);
   fChain->SetBranchAddress("jetMVAsExt_cutBased", &jetMVAsExt_cutBased, &b_jetMVAsExt_cutBased);
   fChain->SetBranchAddress("jetWPLevelsExt_cutBased", &jetWPLevelsExt_cutBased, &b_jetWPLevelsExt_cutBased);
   fChain->SetBranchAddress("jetMVAsExt_philv1", &jetMVAsExt_philv1, &b_jetMVAsExt_philv1);
   fChain->SetBranchAddress("jetWPLevelsExt_philv1", &jetWPLevelsExt_philv1, &b_jetWPLevelsExt_philv1);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtx3dL", &jetVtx3dL, &b_jetVtx3dL);
   fChain->SetBranchAddress("jetVtx3deL", &jetVtx3deL, &b_jetVtx3deL);
   fChain->SetBranchAddress("jetSoftLeptPt", &jetSoftLeptPt, &b_jetSoftLeptPt);
   fChain->SetBranchAddress("jetSoftLeptPtRel", &jetSoftLeptPtRel, &b_jetSoftLeptPtRel);
   fChain->SetBranchAddress("jetSoftLeptdR", &jetSoftLeptdR, &b_jetSoftLeptdR);
   fChain->SetBranchAddress("jetSoftLeptIdlooseMu", &jetSoftLeptIdlooseMu, &b_jetSoftLeptIdlooseMu);
   fChain->SetBranchAddress("jetSoftLeptIdEle95", &jetSoftLeptIdEle95, &b_jetSoftLeptIdEle95);
   fChain->SetBranchAddress("jetDPhiMETJet", &jetDPhiMETJet, &b_jetDPhiMETJet);
   fChain->SetBranchAddress("jetPuJetIdL", &jetPuJetIdL, &b_jetPuJetIdL);
   fChain->SetBranchAddress("jetPuJetIdM", &jetPuJetIdM, &b_jetPuJetIdM);
   fChain->SetBranchAddress("jetPuJetIdT", &jetPuJetIdT, &b_jetPuJetIdT);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("convP4_x", &convP4_x, &b_convP4_x);
   fChain->SetBranchAddress("convP4_y", &convP4_y, &b_convP4_y);
   fChain->SetBranchAddress("convP4_z", &convP4_z, &b_convP4_z);
   fChain->SetBranchAddress("convP4_E", &convP4_E, &b_convP4_E);
   fChain->SetBranchAddress("convVtx_x", &convVtx_x, &b_convVtx_x);
   fChain->SetBranchAddress("convVtx_y", &convVtx_y, &b_convVtx_y);
   fChain->SetBranchAddress("convVtx_z", &convVtx_z, &b_convVtx_z);
   fChain->SetBranchAddress("convVtxErr_x", &convVtxErr_x, &b_convVtxErr_x);
   fChain->SetBranchAddress("convVtxErr_y", &convVtxErr_y, &b_convVtxErr_y);
   fChain->SetBranchAddress("convVtxErr_z", &convVtxErr_z, &b_convVtxErr_z);
   fChain->SetBranchAddress("convPairMomentum_x", &convPairMomentum_x, &b_convPairMomentum_x);
   fChain->SetBranchAddress("convPairMomentum_y", &convPairMomentum_y, &b_convPairMomentum_y);
   fChain->SetBranchAddress("convPairMomentum_z", &convPairMomentum_z, &b_convPairMomentum_z);
   fChain->SetBranchAddress("convRefittedMomentum_x", &convRefittedMomentum_x, &b_convRefittedMomentum_x);
   fChain->SetBranchAddress("convRefittedMomentum_y", &convRefittedMomentum_y, &b_convRefittedMomentum_y);
   fChain->SetBranchAddress("convRefittedMomentum_z", &convRefittedMomentum_z, &b_convRefittedMomentum_z);
   fChain->SetBranchAddress("convNTracks", &convNTracks, &b_convNTracks);
   fChain->SetBranchAddress("convPairInvMass", &convPairInvMass, &b_convPairInvMass);
   fChain->SetBranchAddress("convPairCotThetaSep", &convPairCotThetaSep, &b_convPairCotThetaSep);
   fChain->SetBranchAddress("convEoverP", &convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convDistOfMinApproach", &convDistOfMinApproach, &b_convDistOfMinApproach);
   fChain->SetBranchAddress("convDPhiTrksAtVtx", &convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
   fChain->SetBranchAddress("convDPhiTrksAtEcal", &convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
   fChain->SetBranchAddress("convDEtaTrksAtEcal", &convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
   fChain->SetBranchAddress("convDxy", &convDxy, &b_convDxy);
   fChain->SetBranchAddress("convDz", &convDz, &b_convDz);
   fChain->SetBranchAddress("convLxy", &convLxy, &b_convLxy);
   fChain->SetBranchAddress("convLz", &convLz, &b_convLz);
   fChain->SetBranchAddress("convZofPrimVtxFromTrks", &convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
   fChain->SetBranchAddress("convNHitsBeforeVtx_0", &convNHitsBeforeVtx_0, &b_convNHitsBeforeVtx_0);
   fChain->SetBranchAddress("convNHitsBeforeVtx_1", &convNHitsBeforeVtx_1, &b_convNHitsBeforeVtx_1);
   fChain->SetBranchAddress("convNSharedHits", &convNSharedHits, &b_convNSharedHits);
   fChain->SetBranchAddress("convValidVtx", &convValidVtx, &b_convValidVtx);
   fChain->SetBranchAddress("convMVALikelihood", &convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convChi2", &convChi2, &b_convChi2);
   fChain->SetBranchAddress("convChi2Probability", &convChi2Probability, &b_convChi2Probability);
   fChain->SetBranchAddress("convTk1Dz", &convTk1Dz, &b_convTk1Dz);
   fChain->SetBranchAddress("convTk2Dz", &convTk2Dz, &b_convTk2Dz);
   fChain->SetBranchAddress("convTk1DzErr", &convTk1DzErr, &b_convTk1DzErr);
   fChain->SetBranchAddress("convTk2DzErr", &convTk2DzErr, &b_convTk2DzErr);
   fChain->SetBranchAddress("convCh1Ch2", &convCh1Ch2, &b_convCh1Ch2);
   fChain->SetBranchAddress("convTk1D0", &convTk1D0, &b_convTk1D0);
   fChain->SetBranchAddress("convTk1Pout", &convTk1Pout, &b_convTk1Pout);
   fChain->SetBranchAddress("convTk1Pin", &convTk1Pin, &b_convTk1Pin);
   fChain->SetBranchAddress("convTk2D0", &convTk2D0, &b_convTk2D0);
   fChain->SetBranchAddress("convTk2Pout", &convTk2Pout, &b_convTk2Pout);
   fChain->SetBranchAddress("convTk2Pin", &convTk2Pin, &b_convTk2Pin);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void EventTree::InitSkim()
{
  fChainSkim = fChain->CloneTree(0);
}
void EventTree::SetDirectory( TFile *fSkimFile )
{
  fSkimFile->mkdir("ggNtuplizer");
  fSkimFile->cd("ggNtuplizer");
  fChainSkim->SetDirectory( fSkimFile );
}
void EventTree::FillSkim( TFile *fSkimFile)
{
  //fSkimFile->cd("ggNtuplizer");
  fChainSkim->Fill();
}
void EventTree::AutoSave()
{
  fChainSkim->AutoSave();
}
void EventTree::Write()
{
  fChainSkim->Write();
}
#endif // #ifdef EventTree_cxx
