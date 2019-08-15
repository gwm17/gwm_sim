/*Reconstruct.h
 *Class for implementing kinematic reconstruction of reaction events in ANASEN detector
 *Should be customized for each unique experiment
 *For use in the analyzer class 
 *
 * Gordon M. -- April 2019
 * Based on previous version by M. Anastasiou
 */

#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H

#include <TROOT.h>
#include <vector> 
#include <TList.h>
#include "LookUp.h"
#include "nucleus.h"

using namespace std;

class Reconstruct {

  public:
    /*******Class wide variables********/
    nucleus beam, target, recoil, eject, break1, break2;
    LookUp* ELoss_eject;
    LookUp* ELoss_break1;
    LookUp* ELoss_break2;
    const float maxBeamKE = 17.19;

    /********Construction Area********/
    Reconstruct(nucleus bm, nucleus targ, nucleus ejt, nucleus rec, nucleus br1, nucleus br2);
    ~Reconstruct();

    /*********Functions*************/
    //Double_t CalcElastic(Double_t &EnergyProj,Double_t &EjectTheta,Track Tr,Int_t c);
    void CalcRecoil(track trEject, rcNucleus &rec);
    void CalcRecoil_MultiParticle(track trEject, track trBr1, track trBr2, 
                                  rcNucleus &rec_qval,rcNucleus &rec);
    void CalcRecoil_from_Qvalue(Double_t Qvalue, track trEject, rcNucleus &rec_qval);
};


#endif
