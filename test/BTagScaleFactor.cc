#define BTAGSCALEFACTOR_CC

#include "BTagScaleFactor.h"

double BTagScaleFactor::lfJetCSVM(double x, double jeteta){
  if(jeteta < 0.8) return ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
  if(jeteta < 1.6) return ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
  return ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
}

double BTagScaleFactor::lfJetCSVMmin(double x, double jeteta){
  if(jeteta < 0.8) return ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
  if(jeteta < 1.6) return ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
  return ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
}

double BTagScaleFactor::lfJetCSVMmax(double x, double jeteta){
  if(jeteta < 0.8) return ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
  if(jeteta < 1.6) return ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
  return ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
}

double BTagScaleFactor::lfJetSF(double jetpt, double jeteta)
{
  bool pthigh = false;
  if(jeteta < 1.6 && jetpt > 1000.0) {jetpt = 1000.0; pthigh = true;}
  if(jeteta >= 1.6 && jetpt > 650.0) {jetpt = 650.0; pthigh = true;}

  if(!pthigh){
    if( btagvar012_g == 1 ) return lfJetCSVM(jetpt, jeteta);
    if( btagvar012_g == 0 ) return lfJetCSVMmin(jetpt, jeteta);
    if( btagvar012_g == 2 ) return lfJetCSVMmax(jetpt, jeteta);
  } else{ // in case jetpt is too high:
    if( btagvar012_g == 1 ) return lfJetCSVM(jetpt, jeteta);
    if( btagvar012_g == 0 ) return 2.0*lfJetCSVMmin(jetpt, jeteta) - lfJetCSVM(jetpt, jeteta);
    if( btagvar012_g == 2 ) return 2.0*lfJetCSVMmax(jetpt, jeteta) - lfJetCSVM(jetpt, jeteta);
  }
  std::cout << "should not be here!" << std::endl;
  return 1.0;

}

double BTagScaleFactor::bSFerr(double jetpt)
{
  //Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  static double ptmax[16] = {30, 40, 50, 60, 
                             70, 80, 100, 120, 
                             160, 210, 260, 320, 
                             400, 500, 600, 800};
  static double SFb_error[17] = {0.0415694, 0.023429, 0.0261074, 0.0239251, 
                                 0.0232416, 0.0197251, 0.0217319, 0.0198108, 
                                 0.0193, 0.0276144, 0.020583, 0.026915, 
                                 0.0312739, 0.0415054, 0.074056, 0.0598311, 
                                 0.0598311*2};

  int ptInd = 0;
  if( btagvar012_g == 1) return 0.0;

  for(int ipt=0; ipt<16; ipt++) 
    if(jetpt > ptmax[ipt]) ptInd++;

  if( btagvar012_g == 0 ) return -1.0*SFb_error[ptInd];
  if( btagvar012_g == 2 ) return  SFb_error[ptInd];

  return 0.0;

}

double BTagScaleFactor::bJetSF(double jetpt){
  return (0.939158+(0.000158694*jetpt))+(-2.53962e-07*(jetpt*jetpt)) + bSFerr(jetpt);
}

double BTagScaleFactor::cJetSF(double jetpt){
  return (0.939158+(0.000158694*jetpt))+(-2.53962e-07*(jetpt*jetpt)) + 2.0*bSFerr(jetpt);
}
