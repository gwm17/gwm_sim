#ifndef BAD_DET_LIST_H
#define BAD_DET_LIST_H

#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

class BadDetectorList {
  public:
    BadDetectorList();
    ~BadDetectorList();

    Int_t GetTotal();
    Bool_t IsBad(Int_t);
    Bool_t IsBad(Int_t, Int_t);

  private:
    Int_t total;
    Int_t **DetList;
};

#endif
