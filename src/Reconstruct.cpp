/*Reconstruct.cpp
 *Class for implementing kinematic reconstruction of reaction events in ANASEN detector
 *Should be customized for each unique experiment
 *For use in the analyzer class 
 *
 * Gordon M. -- April 2019
 * Based on previous version by M. Anastasiou
 */

#include "Reconstruct.h"
#include "TLorentzVector.h"
#include <iostream>
#include "TMath.h"
#include "nucFuncs.h"

using namespace std;

//Constructor
Reconstruct::Reconstruct()
{
  ELoss_eject = NULL;
  ELoss_break2 = NULL;
  ELoss_break1 = NULL;
}

//Destructor
Reconstruct::~Reconstruct() {
  delete ELoss_eject;
  delete ELoss_break1;
  delete ELoss_break2;
}

void Reconstruct::SetParams(nucleus bm,nucleus t,nucleus e,nucleus r,nucleus b1,nucleus b2) {
  beam = bm;
  target = t;
  recoil = r;
  eject = e;
  break1 = b1;
  break2 = b2;
}

void Reconstruct::ResetRC() {
  resetNucleus(beam,1); resetNucleus(target,1); resetNucleus(recoil,1);
  resetNucleus(eject,1); resetNucleus(break1,1); resetNucleus(break2,1); 
}

/*Elastic()
 *Reconstructs an elastic scattering event. 
 *Returns the energy of the ejectile at the reaction point, the energy of the incoming 
 *particle, and the angle of the ejectile path. NEEDS SIM MODIFICATION
 */
/*Double_t Reconstruct::CalcElastic(Double_t &EnergyProj, Double_t &EjectTheta, Track Tr, Int_t c) {

  Double_t E_eject_si = Tr.TrEvent[c].SiEnergy;
  Double_t d_eject = Tr.TrEvent[c].PathLength;
  Double_t theta = Tr.TrEvent[c].Theta;

  //Get the inital energy of the ejectile
  Double_t E_eject_rxn = ELoss_eject->GetInitialEnergy(E_eject_si, d_eject, 0.1);

  //Calculate the energy of the proj initially, and the angle of the scattered proj
  EnergyProj = ((m_beam+m_ejectile)*(m_beam+m_ejectile)*E_eject_rxn)/
                   ((4.0*theta*theta*m_beam*m_ejectile));
  EjectTheta = asin((sin(theta)*sqrt(m_ejectile/m_beam))/(sqrt(EnergyProj/E_eject_rxn-1.0)));
  return E_eject_rxn;
}*/

/*CalcRecoil()
 *Reconstructs a recoil nucleus using the tracking information of a single ejectile, and assumed
 *knowledge of the beam and target. Should be compared to CaclRecoil_from_Qvalue, which does not
 *rely as much on the tracking to calculate the same value.
 */                            
void Reconstruct::CalcRecoil(track trEject, rcNucleus &rec) {

  TLorentzVector Eject_LV(0., 0., 0., 0.), Beam_LV(0., 0., 0., 0.), Target_LV(0., 0., 0., 0.),
                 Parent_LV(0., 0., 0., 0.), Recoil_LV(0., 0., 0., 0.);

  Float_t E_eject_si = trEject.SiE;
  Float_t d_eject = trEject.path;
  Float_t theta_eject = trEject.theta;
  Float_t phi_eject = trEject.SiPhi;
  Float_t E_eject_rxn = ELoss_eject->GetLookupEnergy(E_eject_si, -d_eject);
  Float_t E_eject_tot = E_eject_rxn + eject.mass;
  Float_t P_eject = sqrt(E_eject_tot*E_eject_tot - eject.mass*eject.mass);

  Eject_LV.SetPxPyPzE(P_eject*sin(theta_eject)*cos(phi_eject), 
                      P_eject*sin(theta_eject)*sin(phi_eject),
                      P_eject*cos(theta_eject),
                      E_eject_tot);
   
  Float_t BeamE_tot = trEject.beamKE + beam.mass;
  Float_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - beam.mass*beam.mass);
  Float_t TargetE = target.mass;

  Beam_LV.SetPxPyPzE(0.0, 0.0, BeamP_z, BeamE_tot);
  Target_LV.SetPxPyPzE(0.0,0.0,0.0, TargetE);
  
  Parent_LV = Beam_LV + Target_LV;

  Recoil_LV = Parent_LV - Eject_LV;

  rec.IntPoint = trEject.IntPoint;
  rec.beamKE = trEject.beamKE;
  rec.beamPz = Beam_LV.Pz();
  rec.theta = Recoil_LV.Theta()*180.0/TMath::Pi();
  rec.phi = Recoil_LV.Phi()*180.0/TMath::Pi();
  rec.KE = Recoil_LV.E()-Recoil_LV.M();
  rec.Ex = Recoil_LV.M() - recoil.mass;
  rec.E = Recoil_LV.E();
  rec.px = Recoil_LV.Px();
  rec.py = Recoil_LV.Py();
  rec.pz = Recoil_LV.Pz();
  rec.p = Recoil_LV.P();
  
}

/*CalcRecoil_MultiParticle()
 *This is the big boy... Reconstructs an event from 3 outgoing particles, in this specific ver,
 *1 p and 2 alpha. Then calculates the various properties of the residual nucleus from the 
 *detected particles. 
 */
void Reconstruct::CalcRecoil_MultiParticle(track trEject, track trBr1, track trBr2, 
                                           rcNucleus &rec_qval, rcNucleus &rec) 
{
  TLorentzVector Beam_LV(0.,0.,0.,0.), Target_LV(0.,0.,0.,0.),
                 Recoil_LV(0.,0.,0.,0.), eject_LV(0.,0.,0.,0.), break1_LV(0.,0.,0.,0.),
                 break2_LV(0.,0.,0.,0.), Sum_eject_LV(0.,0.,0.,0.),
                 Recoil_LV_beam_qval(0.,0.,0.,0.), Beam_LV_qval(0.,0.,0.,0.);
 
  TVector3 eject_V(0.,0.,0.), break1_V(0.,0.,0.), break2_V(0.,0.,0.), Sum_eject_V(0.,0.,0.);
 
  Target_LV.SetPxPyPzE(0.0,0.0,0.0, target.mass);
  
  Double_t BeamKE_avg = 0.0;
  
  Double_t E_eject1_rxn = ELoss_eject->GetLookupEnergy(trEject.SiE, -trEject.path);
  Double_t E_eject_tot = E_eject1_rxn + eject.mass;
  Double_t P_eject = sqrt(E_eject_tot*E_eject_tot - eject.mass*eject.mass);
  Double_t p_x = P_eject*cos(trEject.SiPhi)*sin(trEject.theta);
  Double_t p_y = P_eject*sin(trEject.SiPhi)*sin(trEject.theta);
  Double_t p_z = P_eject*cos(trEject.theta);
  eject_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
  eject_V.SetXYZ(p_x, p_y, p_z);

  Double_t E_break1_rxn = ELoss_break1->GetLookupEnergy(trBr1.SiE, -trBr1.path);
  E_eject_tot = E_break1_rxn + break1.mass;
  P_eject = sqrt(E_eject_tot*E_eject_tot - break1.mass*break1.mass);
  p_x = P_eject*cos(trBr1.SiPhi)*sin(trBr1.theta);
  p_y = P_eject*sin(trBr1.SiPhi)*sin(trBr1.theta);
  p_z = P_eject*cos(trBr1.theta);
  break1_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
  break1_V.SetXYZ(p_x, p_y, p_z);
  BeamKE_avg += trBr1.beamKE;

  Double_t E_break2_rxn = ELoss_break2->GetLookupEnergy(trBr2.SiE, -trBr2.path);
  E_eject_tot = E_break2_rxn+break2.mass;
  P_eject = sqrt(E_eject_tot*E_eject_tot - break2.mass*break2.mass);
  p_x = P_eject*cos(trBr2.SiPhi)*sin(trBr2.theta);
  p_y = P_eject*sin(trBr2.SiPhi)*sin(trBr2.theta);
  p_z = P_eject*cos(trBr2.theta);
  break2_LV.SetPxPyPzE(p_x, p_y, p_z, E_eject_tot);
  break2_V.SetXYZ(p_x, p_y, p_z);
  BeamKE_avg += trBr2.beamKE;
 
  Sum_eject_LV = eject_LV + break1_LV + break2_LV;
  Sum_eject_V = eject_V + break1_V + break2_V;
  
  Double_t Beam_KE = Sum_eject_LV.E()-target.mass-beam.mass;//Qvalue included implicitly
 
  BeamKE_avg = BeamKE_avg/2.0;
  Double_t BeamE_tot = BeamKE_avg + beam.mass; 
  Double_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - beam.mass*beam.mass);
  Double_t BeamE_tot_qval = Beam_KE + beam.mass;
 
  Double_t BeamP_z_qval = sqrt(BeamE_tot_qval*BeamE_tot_qval - beam.mass*beam.mass);
  Beam_LV_qval.SetPxPyPzE(0.0,0.0, BeamP_z_qval, BeamE_tot_qval);
  Recoil_LV_beam_qval = break1_LV+break2_LV;
 
  Beam_LV.SetPxPyPzE(0.0,0.0, BeamP_z, BeamE_tot);
  Recoil_LV = break1_LV+break2_LV;
 
  rec.beamPz = Beam_LV.Pz();
  rec_qval.beamPz = Beam_LV_qval.Pz();
  rec_qval.beamKE = Beam_KE;
  rec.beamKE = BeamKE_avg;
  rec.ejectKE = E_eject1_rxn;
  rec.break1KE = E_break1_rxn;
  rec.break2KE = E_break2_rxn;
  rec_qval.theta = Recoil_LV_beam_qval.Theta();
  rec_qval.phi = Recoil_LV_beam_qval.Phi();
  rec_qval.KE = Recoil_LV_beam_qval.E()-Recoil_LV_beam_qval.M();
  rec_qval.Ex = Recoil_LV_beam_qval.M() - recoil.mass;
  rec.theta = Recoil_LV.Theta();
  rec.phi = Recoil_LV.Phi();
  rec.KE = Recoil_LV.E() - Recoil_LV.M();
  rec.Ex = Recoil_LV.M() - recoil.mass;


  if(eject.mass != break2.mass) {
    rec.mass_sq1 = (Recoil_LV.M()*Recoil_LV.M())/1e6;
    TLorentzVector temp1 = eject_LV+break1_LV;
    TLorentzVector temp2 = eject_LV+break2_LV;
    float ex1 = temp1.M()-4667.6163636366931;
    float ex2 = temp2.M()-4667.6163636366931;
    TLorentzVector secondary_LV;
    if(ex1<ex2) {
      secondary_LV = temp1;
    } else {
      secondary_LV = temp2;
    }
    rec.mass_sq2 = (secondary_LV.M()*secondary_LV.M())/1e6;
  } else {
    TLorentzVector a2_LV = eject_LV+break2_LV;
    rec.mass_sq1 = (a2_LV.M()*a2_LV.M())/1e6;
    TLorentzVector temp1 = eject_LV+break1_LV;
    float ex1 = temp1.M()-4667.6163636366931;
    float ex2 = Recoil_LV.M()-4667.6163636366931;
    TLorentzVector secondary_LV;
    if(ex1<ex2) {
      secondary_LV = temp1;
    } else {
      secondary_LV = Recoil_LV;
    }
    rec.mass_sq2 = (secondary_LV.M()*secondary_LV.M())/1e6;
  }

  rec.E = Recoil_LV.E();
  rec.p = Recoil_LV.P();
  rec.pz = Recoil_LV.Pz();
  rec.py = Recoil_LV.Py();
  rec.px = Recoil_LV.Px();
  rec_qval.E = Recoil_LV_beam_qval.E();
  rec_qval.p = Recoil_LV_beam_qval.P();
  rec_qval.pz = Recoil_LV_beam_qval.Pz();
  rec_qval.py = Recoil_LV_beam_qval.Py();
  rec_qval.px = Recoil_LV_beam_qval.Px();
 
  rec.deltaPhi = abs(break1_LV.Phi()-break2_LV.Phi());
  rec_qval.deltaPhi = abs(break1_LV.Phi()-break2_LV.Phi());
  if(rec.deltaPhi > TMath::Pi() && rec.deltaPhi <= 2*TMath::Pi()) {
    rec.deltaPhi = 2*TMath::Pi() - rec.deltaPhi;
    rec_qval.deltaPhi = 2*TMath::Pi() - rec.deltaPhi;
  } 
  rec.deltaTheta = abs(break1_LV.Theta() - break2_LV.Theta());
  rec_qval.deltaTheta = abs(break1_LV.Theta() - break2_LV.Theta());
}

/*CalcRecoil_from_Qvalue()
 *Reconstructs a recoil nucleus from a single ejectile assuming knowledge of the beam,
 *reaction Qvalue, and target. This is mainly meant to be a comparison to CalcRecoil() to 
 *ensure that the tracking isn't super screwed up. 
 */
void Reconstruct::CalcRecoil_from_Qvalue(Double_t Qvalue, track trEject, rcNucleus &rec_qval) {

  Double_t E_eject_rxn = ELoss_eject->GetLookupEnergy(trEject.SiE, -trEject.path);
  Double_t E_eject_tot = E_eject_rxn + eject.mass;
  Double_t P_eject = sqrt(E_eject_tot*E_eject_tot - eject.mass*eject.mass);
  Double_t px = P_eject*sin(trEject.theta)*cos(trEject.SiPhi),
           py = P_eject*sin(trEject.theta)*sin(trEject.SiPhi),
           pz = P_eject*cos(trEject.theta);

  TLorentzVector Eject_LV(0.,0.,0.,0.), Recoil_LV(0.,0.,0.,0.), Beam_LV(0.,0.,0.,0.),
                 Target_LV(0.,0.,0.,0.), Parent_LV(0.,0.,0.,0.);
  Double_t a = (beam.mass-recoil.mass);
  Double_t b = (recoil.mass+eject.mass);
  Double_t c = -Qvalue*recoil.mass;
  Double_t d = 2.0*sqrt(beam.mass*eject.mass);
  Double_t cos2 = cos(trEject.theta)*cos(trEject.theta);
  Double_t a1 = a*a;
  Double_t a2 = (2.0*a*(b*E_eject_rxn+c)-d*d*E_eject_rxn*cos2);
  Double_t a3 = (b*E_eject_rxn+c)*(b*E_eject_rxn+c);
  Double_t s1 = (-a2+sqrt(a2*a2-4.0*a1*a3))/(2.0*a1);
  Double_t s2 = (-a2-sqrt(a2*a2-4.0*a1*a3))/(2.0*a1);

  if(s1<s2) swap(s1,s2);
  Double_t Beam_KE = s2;
  if(TMath::IsNaN(s2)) {
    return;
  }
  
  Eject_LV.SetPxPyPzE(px, py, pz, E_eject_tot);
  Double_t E_beam_tot = Beam_KE + beam.mass;
  Double_t pz_beam = sqrt(E_beam_tot*E_beam_tot-beam.mass*beam.mass);
  Beam_LV.SetPxPyPzE(0.,0.,pz_beam, E_beam_tot);
  Target_LV.SetPxPyPzE(0.,0.,0.,target.mass);
  Parent_LV = Beam_LV+Target_LV;
  Recoil_LV = Parent_LV-Eject_LV;

  rec_qval.beamKE = Beam_KE;
  rec_qval.Ex = Recoil_LV.M()-recoil.mass;
  rec_qval.E = Recoil_LV.E();
  rec_qval.KE = Recoil_LV.E()-Recoil_LV.M();
  rec_qval.beamPz = Beam_LV.Pz();
  rec_qval.p = Recoil_LV.P();
  rec_qval.px = Recoil_LV.Px();
  rec_qval.py = Recoil_LV.Py();
  rec_qval.pz = Recoil_LV.Pz();
  rec_qval.theta = Recoil_LV.Theta();
  rec_qval.phi = Recoil_LV.Phi();
}

