#include "nucFuncs.h"


void resetNucleus(nucleus &nuc, int massFlag) {
  if(massFlag) {
    nuc.mass = -10;
  }
  nuc.px = -1e6; nuc.py = -1e6; nuc.pz = -1e6;
  nuc.E = -1e6; nuc.KE = -1e6; nuc.Ex = -1.0;
  nuc.Ecm = -1e6;
  nuc.theta = -1e6; nuc.phi = -1e6;
  nuc.eject_theta_cm = -1e6;
  nuc.a2_mass_sq = -1.0; nuc.ap_mass_sq = -1.0;
  nuc.detID = -1;
  nuc.wireID = -1;
  nuc.SX3xyz.resize(0); nuc.Q3xyz.resize(0);
  nuc.PCxyz.resize(0);
};

void resetRcNucleus(rcNucleus &nuc) {
  nuc.px = -10.; nuc.py = -10.; nuc.pz = -10.; nuc.p = -10.; nuc.E = -10.; nuc.KE = -10.;
  nuc.theta = -10.; nuc.phi = -10.; nuc.Ex = -10.;
  nuc.beamKE = -10.; nuc.beamPz = -10.; nuc.IntPoint = -10.;
  nuc.deltaPhi = -10; nuc.deltaTheta = -10.; 
  nuc.mass_sq1 = -10; nuc.mass_sq2 = -10;
  nuc.ejectKE = -10; nuc.break1KE = -10; nuc.break2KE = -10;
};

void resetTrack(track &tr) {
  tr.IntPoint = -10.; tr.SiR = -10.; tr.PCR = -10.; tr.SiZ = -10.; tr.PCZ = -10.;
  tr.theta = -10.; tr.SiPhi = -10.; tr.path = -10.; tr.beamKE = -10.; tr.SiE = -10.;
};
