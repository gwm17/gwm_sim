#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <TROOT.h>
#include <TLorentzVector.h>
#include <vector>
#include <string>
#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include "thetaphi.h"
#include "QQQ.h"
#include "SX3.h"
#include "qqq_hit.h"
#include "LookUp.h"
#include "BadDetectorList.h"
#include "nucleus.h"

using namespace std;

class calculator {

  public:
    calculator(Float_t mb, Float_t mt, Float_t mr, Float_t me, Float_t mbr1, Float_t mbr2,
               Float_t rxMin, Float_t rxMax);
    ~calculator();
    void ResetNucleus(nucleus &nuc);
    void InitializeBeam(Float_t minKE, Float_t maxKE);
    bool GetSX3Coords(nucleus &nuc);
    bool GetQ3Coords(nucleus &nuc);
    track GetIntPoint(nucleus nuc, LookUp *nuc_eloss);
    void Step1_CalcRxn();
    bool Step2_CalcBreakup();
    LookUp *beam_eloss, *eject_eloss, *break1_eloss, *break2_eloss;
    nucleus beam, target, eject, recoil, break1, break2;
    Float_t x_rxn, y_rxn, z_rxn;
    Float_t minBeamKE, maxBeamKE;
    Float_t maxEx, minEx;

  private:
    TRandom3 *random;

    void GetPCCoords(nucleus &nuc);

    TLorentzVector targ_LV, beam_LV, recoil_LV, eject_LV, break1_LV, break2_LV;
    QQQ Q3; SX3 S3;

    BadDetectorList BadDetList;
    const float ana_length = 55.0545; //length of ANASEN for checks

};

#endif
