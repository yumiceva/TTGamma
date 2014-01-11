#ifndef EVTCLEANING_H
#define EVTCLEANING_H

#include "BaseSelector.h"
#include<iostream>

using namespace std;

class EvtCleaning: public BaseSelector {

public:
  EvtCleaning() {}
  ~EvtCleaning() {}

  bool Pass()
  {
    int HBHENoiseFilter = fReader->metFilters[1];
    int HcalLaserFilter = fReader->metFilters[2];
    int EcalDeadCellFilter=fReader->metFilters[3];
    int TrackingFailureFilter=fReader->metFilters[4];
    int EEBadScFilter=fReader->metFilters[5];
    int EcalLaserFilter=fReader->metFilters[6];
    int Manystripclus53X=fReader->metFilters[7];
    int Toomanystripclus53X=fReader->metFilters[8];
    int LogErrorTooManyClusters=fReader->metFilters[9];

    if (fVerbose) {
      cout << "Filters:" <<endl;
      cout << "HBHENoiseFilter= "<< HBHENoiseFilter << endl;
      cout << "HcalLaserFilter= "<< HcalLaserFilter << endl;
      cout << "EcalDeadCellFilter= "<< EcalDeadCellFilter << endl;
      cout << "TrackingFailureFilter= "<< TrackingFailureFilter << endl;
      cout << "EEBadScFilter= "<< EEBadScFilter << endl;
      cout << "EcalLaserFilter= "<< EcalLaserFilter<< endl;
      cout << "Manystripclus53X= "<< Manystripclus53X << endl;
      cout << "Toomanystripclus53X= "<< Toomanystripclus53X << endl;
      cout << "LogErrorTooManyClusters= "<< LogErrorTooManyClusters << endl;
    }
    if ( HBHENoiseFilter && HcalLaserFilter && EcalDeadCellFilter && TrackingFailureFilter && EEBadScFilter )
      return true;
    else
      return false;
  }

};
#endif
