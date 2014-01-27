#ifndef MUONSELECTOR_H
#define MUONSELECTOR_H

#include "BaseSelector.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class MuonSelector: public BaseSelector {

public:
  MuonSelector() {}
  ~MuonSelector() {}

  bool PassLoose(int imu)
  {
    // check muon type
    static const unsigned int GlobalMuon     =  1<<1;
    static const unsigned int TrackerMuon    =  1<<2;
    static const unsigned int PFMuon =  1<<5;
    bool isGlobalMuon  = fReader->muType->at(imu) & GlobalMuon;
    bool isTrackerMuon = fReader->muType->at(imu) & TrackerMuon;
    bool isPFMuon      = fReader->muType->at(imu) & PFMuon;

    frelIso = (fReader->muPFIsoR04_CH->at(imu) + fReader->muPFIsoR04_NH->at(imu))/ fReader->muPt->at(imu);
    // delta beta corrections
    //  I = [sumChargedHadronPt+ max(0.,sumNeutralHadronPt+sumPhotonPt-0.5sumPUPt]/pt
    frelIsocorr = ( fReader->muPFIsoR04_CH->at(imu) + fmax(0.0, fReader->muPFIsoR04_NH->at(imu) +fReader->muPFIsoR04_Pho->at(imu) -0.5*fReader->muPFIsoR04_PU->at(imu)) )/ fReader->muPt->at(imu);
    // loose muon selection
    if ( fReader->muPt->at(imu) > 10.0 &&
         fabs( fReader->muEta->at(imu) ) < 2.5 &&
         frelIsocorr < 0.2 &&
         isPFMuon && ( isGlobalMuon || isTrackerMuon) )
      return true;
    else
      return false;
  }

  bool PassTight(int imu)
  {
    // check muon type
    static const unsigned int GlobalMuon     =  1<<1;
    static const unsigned int TrackerMuon    =  1<<2;
    static const unsigned int PFMuon =  1<<5;
    bool isGlobalMuon  = fReader->muType->at(imu) & GlobalMuon;
    bool isTrackerMuon = fReader->muType->at(imu) & TrackerMuon;
    bool isPFMuon      = fReader->muType->at(imu) & PFMuon;

    if ( PassLoose( imu ) &&
         fReader->muPt->at(imu) > 26.0 &&
         fabs(fReader->muEta->at(imu)) < 2.1 &&
         fReader->muChi2NDF->at(imu) < 10 &&
         fReader->muNumberOfValidTrkLayers->at(imu) > 5 &&
         fReader->muNumberOfValidMuonHits->at(imu) > 0 &&
         fReader->muD0->at(imu) < 0.2 &&
         fabs( fReader->muDz->at(imu) ) < 0.5 && //check this
         fReader->muNumberOfValidPixelHits->at(imu) > 0 &&
         fReader->muStations->at(imu) > 1 &&
         frelIsocorr < 0.12 &&
         isPFMuon && isGlobalMuon && isTrackerMuon )

      return true;
    else
      return false;
  }

  float get_relIso() { return frelIso; }
  float get_relIsocorr() { return frelIsocorr; }

  private:
    float frelIso;
    float frelIsocorr;
};
#endif
