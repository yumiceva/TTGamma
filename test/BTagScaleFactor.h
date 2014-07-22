#ifndef BTAGSCALEFACTOR_H
#define BTAGSCALEFACTOR_H

#include "BaseScaleFactor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class BTagScaleFactor: public BaseScaleFactor {

public:

  int btagvar012_g;

  BTagScaleFactor() 
  {
    btagvar012_g = 0;
  }
  ~BTagScaleFactor() {}

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
  // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
  // weight for >=1 btag :  1 - prod(1-SFi) over all b-tagged jets

  float lfJetCSVM(float x, float jeteta);
  float lfJetCSVMmin(float x, float jeteta);
  float lfJetCSVMmax(float x, float jeteta);
  float lfJetSF(float jetpt, float jeteta);
  float bSFerr(float jetpt);
  float bJetSF(float jetpt);
  float cJetSF(float jetpt);

    private:

};

#endif

