#ifndef MUONSCALEFACTOR_H
#define MUONSCALEFACTOR_H

#include "BaseScaleFactor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

class MuonScaleFactor: public BaseScaleFactor {

public:
  MuonScaleFactor() {}
  ~MuonScaleFactor() {}

  float GetSF( float eta )
  {

    float SF = 0.;
    // taken from https://indico.cern.ch/event/257000/material/slides/0?contribId=2
    if ( fabs(eta) >= 0   && fabs(eta) < 0.9 ) SF = 0.9837;
    if ( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) SF = 0.9656;
    if ( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) SF = 0.9962;

    fSF = SF;
    return SF;

  }

  float GetSigmaSF( float eta )
  {

    float sigma = 0.;
    
    if ( fabs(eta) >= 0   && fabs(eta) < 0.9 ) sigma = 0.9837*0.02;
    if ( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) sigma = 0.9656*0.02;
    if ( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) sigma = 0.9962*0.02;

    return sigma;

  }
  

  private:

};
#endif
