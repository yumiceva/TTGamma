#ifndef BASESCALEFACTOR_H
#define BASESCALEFACTOR_H

class BaseScaleFactor {

public:
  BaseScaleFactor () {}
  ~BaseScaleFactor () {}

  void Init()
  {
    fSF = 0.;
  }
  float GetSF() {}
  float GetSigmaSF() {}

protected:
  float fSF;
  float fsigma;
};
#endif
