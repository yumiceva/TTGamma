#ifndef ELECTRONSELECTOR_H
#define ELECTRONSELECTOR_H

#include "BaseSelector.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

using namespace std;

class ElectronSelector: public BaseSelector {

public:
  ElectronSelector() {}
  ~ElectronSelector() {}

  bool PassLoose(int ie)
  {
    float AEff03 = 0.00;

    if (!fIsMC) {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, fReader->eleSCEta->at(ie), ElectronEffectiveArea::kEleEAData2011);
    }
    else {
      AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, fReader->eleSCEta->at(ie), ElectronEffectiveArea::kEleEAFall11MC);
    }

    frelIso = (fReader->elePFChIso03->at(ie) + fReader->elePFNeuIso03->at(ie) + fReader->elePFPhoIso03->at(ie))/ fReader->elePt->at(ie);
    frelIsocorr = ( fReader->elePFChIso03->at(ie) + fmax(0.0, fReader->elePFNeuIso03->at(ie) +fReader->elePFPhoIso03->at(ie) - fReader->rho25_elePFiso*AEff03) )/ fReader->elePt->at(ie);

    // loose electrons
    if ( fReader->elePt->at(ie) > 20.0 &&
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
           !( (1.4442 < fabs( fReader->eleSCEta->at(ie) )) && fabs( fReader->eleSCEta->at(ie)) < 1.5660) &&
           fReader->eleD0->at(ie) < 0.02 &&
           fReader->eleIDMVATrig->at(ie) > 0.5 && fReader->eleIDMVATrig->at(ie) < 1.0 &&
           fReader->eleMissHits->at(ie) <= 0 &&
           frelIsocorr < 0.1 )
      
        return true;
      else
        return false;
    }
  }

  float get_relIso() { return frelIso; }
  float get_relIsocorr() { return frelIsocorr; }

  private:
    float frelIso;
    float frelIsocorr;
};
#endif
