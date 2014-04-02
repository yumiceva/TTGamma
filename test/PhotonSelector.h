#ifndef PHOTONSELECTOR_H
#define PHOTONSELECTOR_H

#include "BaseSelector.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class PhotonSelector: public BaseSelector {

public:
  PhotonSelector() {}
  ~PhotonSelector() {}

  bool PassLoose(int ie)
  {

    double eta = fReader->phoEta->at(ie);
    if ( fabs(eta) > 1.5 ) fInBarrel = false;
    else fInBarrel = true;

    float AEff03 = 0.00;

    //if (!fIsMC) {
    
    // loose photons
    if ( fReader->phoEt->at(ie) > 25.0 &&
         !( ( 1.4442 < fabs(eta) ) && fabs(eta) < 1.5660) &&
         fabs(fReader->eleSCEta->at(ie)) < 2.5 &&
         fReader->eleIDMVATrig->at(ie) > 0.0 && fReader->eleIDMVATrig->at(ie) < 1.0 &&
         frelIsocorr < 0.15 )
 
      return true;
    else
      return false;
  }

  bool PassTight(int ie)
  {
    if ( PassLoose( ie ) ) {

      if ( fReader->elePt->at(ie) > 30.0 &&
           fabs(fReader->eleSCEta->at(ie)) < 2.5 &&
           (fReader->phohasPixelSeed_->at(ie) == 0) &&
           (fReader->phoEleVeto_->at(phoInd) == 0) &&
           
           frelIsocorr < 0.1 )
      
        return true;
      else
        return false;
    }
  }

  float get_relIso() { return frelIso; }
  float get_relIsocorr() { return frelIsocorr; }
  bool InBarrel() { retrun fInBarrel; }

  private:
    float frelIso;
    float frelIsocorr;
    bool fInBarrel;
};
#endif
