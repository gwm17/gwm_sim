#ifndef DALITZ_H
#define DALITZ_H

#include <TROOT.h>
#include <TH2.h>
#include <string>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <TGraph.h>
#include <unordered_set>

using namespace std;

class MakeEff {

  public:
    MakeEff();
    ~MakeEff();
    TH2F* makeEffHisto(TH2F *free, TH2F *detect, int Ebin);
  
  private:
    Double_t interpBinContent(Int_t binx, Int_t biny, Int_t gbin);
    bool checkBinValid(Int_t gbin);
    TH2F* heff;
    TH2F* hfree;
    TH2F* hdetect;
    Int_t depth;
    unordered_set<Int_t> interp_set;
    
};

#endif
