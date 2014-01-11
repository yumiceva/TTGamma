#ifndef BASESELECTOR_H
#define BASESELECTOR_H

#include"EventTree.h"

class BaseSelector {

public:
  BaseSelector () {}
  ~BaseSelector () {}

  void Init( EventTree *reader, bool verbose = false, bool IsMC = false)
  {
    fReader = reader;
    fVerbose = verbose;
    fIsMC = IsMC;
  }
  bool Pass() {}

protected:
  EventTree *fReader;
  bool fVerbose;
  bool fIsMC;
};
#endif
