#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct nucleus {
  /***Sim variables***/
  Float_t mass = -1.0, px = -1e6, py = -1e6, pz = -1e6, p = -1e6;
  Float_t E = -1e6, KE = -1e6, Ex = -1.0;
  Float_t Ecm = -1e6;
  Float_t a2_mass_sq = -1.0, ap_mass_sq = -1.0;
  Float_t theta = -1e6, phi = -1e6, theta_cm = -1e6, phi_cm = -1e6;
  Float_t eject_theta_cm = -1e6;
  Int_t detID = -1;
  Int_t wireID = -1;
  vector<Float_t> SX3xyz, PCxyz, Q3xyz;
 
};

struct track {
  /***Tracking variables***/
  Float_t IntPoint = -10., SiR = -10., PCR = -10., SiZ = -10., PCZ = -10.;
  Float_t theta = -10., SiPhi = -10., path = -10., beamKE = -10., SiE = -10.;
};

struct rcNucleus {
  /***Reconstruction variables***/
  Float_t px = -10., py = -10., pz = -10., p = -10., E = -10., KE = -10.,
          theta = -10., phi = -10., Ex = -10.;
  Float_t beamKE = -10., beamPz = -10., IntPoint = -10., deltaPhi = -10, deltaTheta = -10.; 
  Float_t mass_sq1 = -10, mass_sq2 = -10;
  Float_t ejectKE = -10, break1KE = -10, break2KE = -10;
};

#endif
