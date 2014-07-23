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
    btagvar012_g = 1; // central
  }
  ~BTagScaleFactor() {}

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
  // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
  // weight for >=1 btag :  1 - prod(1-SFi) over all b-tagged jets

  double lfJetCSVM(double x, double jeteta);
  double lfJetCSVMmin(double x, double jeteta);
  double lfJetCSVMmax(double x, double jeteta);
  double lfJetSF(double jetpt, double jeteta);
  double bSFerr(double jetpt);
  double bJetSF(double jetpt);
  double cJetSF(double jetpt);

    private:

};

#endif

