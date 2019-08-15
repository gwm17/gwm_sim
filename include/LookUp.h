/*LookUp.h
 *Class designed to take in srim files and calculate the energy loss of a particle 
 *moving through gas a certain distance. Intended for use in ANASEN analysis. 
 *
 *N. Rijal -- June 2016
 *Facelift by G.M. April 2019
 *
 *Modified to provide also the straggling of the particle through the gas for use in ANASEN
 *simulation with the beam particles -- G.M. June 2019 NOTE: This requires that all of the info
 *(columns) from SRIM be included in the eloss file!
 *
 */

#ifndef __LOOKUP_H__
#define __LOOKUP_H__
#include "TROOT.h"
#include "TGraph.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class LookUp{

  public:

    LookUp();
    ~LookUp();
    LookUp(string Eloss_file, Double_t InputMass);
    Double_t GetEnergyLoss(Double_t initial_energy, Double_t distance);
    Double_t GetInitialEnergy(Double_t FinalEnergy, Double_t PathLength, Double_t StepSize);
    Double_t GetFinalEnergy(Double_t InitialEnergy, Double_t PathLength, Double_t StepSize);
    Double_t GetDistance(Double_t InitialE, Double_t FinalE, Double_t StepSize);
    Double_t GetDistance(Double_t InitialE, Double_t FinalE, Double_t StepSize, Int_t Max);
    Double_t GetPathLength(Float_t InitialEnergy, Float_t FinalEnergy, Float_t DeltaT);
    Double_t LoadRange(Float_t energy1);
    Double_t GetTimeOfFlight(Double_t InitialEnergy, Double_t PathLength, Double_t StepSize);
    void SetIonMass(Double_t IonMass);  
    void InitializeLookupTables(Double_t MaximumEnergy, Double_t MaximumDistance, Double_t DeltaE, Double_t DeltaD);
    void PrintLookupTables();
    Double_t GetLookupEnergy(Double_t InitialEnergy, Double_t distance);
    bool GoodELossFile;
    TGraph* EvD;
    void GetStraggle(Double_t z_dist, Double_t &z, Double_t &y);

  private:

    Double_t c;
    Double_t IonMass;
    Double_t IonEnergy;
    Double_t dEdx_e;
    Double_t dEdx_n;
    Double_t range;
    Double_t lgstrg;
    Double_t ltstrg;
    Double_t EtoDtab;
    Double_t DtoEtab; 

    vector<Double_t> IonEnergy_v;
    vector<Double_t> dEdx_e_v;
    vector<Double_t> dEdx_n_v;

    Double_t MaximumEnergy;
    Double_t MaximumDistance;
    Double_t DeltaD;
    Double_t DeltaE;

    vector<Double_t> EtoDtab_v;
    vector<Double_t> DtoEtab_v;

    vector<Double_t> Range;
    vector<Double_t> longStrag;
    vector<Double_t> latStrag;

    int points;
    int last_point;
    int points1;
    int last_point1;
    bool Energy_in_range;
};


#endif 
