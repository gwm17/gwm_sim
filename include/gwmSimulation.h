#ifndef GWM_SIMULATION_H
#define GWM_SIMULATION_H

#include <TROOT.h>
#include <vector>
#include <string>
#include <iostream>
#include <TRandom3.h>
#include <THashTable.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>
#include "calculator.h"
#include "Reconstruct.h"
#include "nucleus.h"


using namespace std;

class simulation {

  public:
    simulation();
    ~simulation();
    void run(char *outName);

  private:
    void MyFill(string name,int binsX,double lowX,double highX,double valueX); 
    void MyFill(string name,int binsX,double lowX,double highX,double valueX,
                            int binsY,double lowY,double highY,double valueY); 
    void GetReactions(string listName);
    void initEloss();
    void resetNucleus(nucleus &nuc);
    void resetNucleus(track &tr);
    void resetNucleus(rcNucleus &rc);
    void resetNuclei();
    Int_t Go2StepCalc(calculator &calc, Int_t nparticles);
    Int_t Go1StepCalc(calculator &calc, Int_t nparticles);
    void getDalitzEff();

    /********Reaction Parameters********/
    vector<Float_t> mt, mb, mr, me, mbr1, mbr2, rxMax, rxMin;
    LookUp *be7_eloss, *he4_eloss, *he3_eloss, *p_eloss;
    Double_t minBeamKE, maxBeamKE;

    /*******Simulation Parameters*******/
    int nBeamParticles, nDetected, nMissed;
    string be7eloss_name, he4eloss_name, he3eloss_name, peloss_name;
    const float m_7be = 6534.1836677282, m_alpha = 3727.37929745092,
                m_p = 938.27206671856, m_3he = 2808.3915032078;
    
    /*******ANASEN Parameters**********/
    const float totalLength = 55.45;  //cm

    /******ROOT storage***************/
    Int_t rxnTag, detectFlag;
    nucleus beam, target, recoil, ejectile, break1, break2;
    track trEjectile, trBreak1, trBreak2;
    rcNucleus rcRecoil, rcRecoil_qval;
    Float_t z_rxn, x_rxn, y_rxn;
    TTree *outTree;
    THashTable *rootObj;

};

#endif
