#include "calculator.h"

using namespace std;

calculator::calculator(Float_t mb, Float_t mt, Float_t mr, Float_t me, Float_t mbr1, 
                       Float_t mbr2, Float_t rxMin, Float_t rxMax) {
  target.mass = mt;
  beam.mass = mb;
  eject.mass = me;
  recoil.mass = mr;
  break1.mass = mbr1;
  break2.mass = mbr2;
  minEx = rxMin;
  maxEx = rxMax;

  random = new TRandom3();
  random->SetSeed();

}

calculator::~calculator() {
  delete beam_eloss;
  delete eject_eloss;
  delete break1_eloss;
  delete break2_eloss;
  delete random;
}

void calculator::ResetNucleus(nucleus &nuc) {
  nuc.px = -1e6;
  nuc.py = -1e6;
  nuc.pz = -1e6;
  nuc.p = -1e6;
  nuc.E = -1e6;
  nuc.KE = -1e6;
  nuc.theta = -1e6;
  nuc.phi = -1e6;
  nuc.detID = -1;
  nuc.wireID = -1;
  nuc.SX3xyz.resize(0);
  nuc.Q3xyz.resize(0);
  nuc.PCxyz.resize(0);
}

void calculator::InitializeBeam(Float_t minKE, Float_t maxKE) {
  beam.KE = random->Uniform(minKE, maxKE);
  Double_t dist_traveled = beam_eloss->GetDistance(maxKE, beam.KE, 0.1);
  Double_t z_strg = 0.0, y_strg = 0.0;
  beam_eloss->GetStraggle(dist_traveled, z_strg, y_strg);
  z_rxn = random->Gaus(dist_traveled, z_strg);
  Double_t x_dist = random->Uniform(-0.5, 0.5);//Assume beam spot of 1cm radius
  Double_t y_dist = random->Uniform(-0.5, 0.5);
  z_rxn = 55.45-z_rxn;
  if (z_rxn <0) {
    return;
  }
  x_rxn = random->Gaus(x_dist, y_strg);
  y_rxn = random->Gaus(y_dist, y_strg);
  maxBeamKE = maxKE;
  minBeamKE = minKE;
}

void calculator::Step1_CalcRxn() {
  beam.E = beam.mass+beam.KE;
  beam.p = sqrt(beam.E*beam.E-beam.mass*beam.mass);
  beam.pz = beam.p;
  beam_LV.SetPxPyPzE(0.,0.,beam.p, beam.E);
  

  target.E = target.mass;
  targ_LV.SetPxPyPzE(0.,0.,0.,target.mass);
  TLorentzVector parent_LV = beam_LV + targ_LV;
  TVector3 boost = parent_LV.BoostVector(); //use to shfit from lab to cm
  parent_LV.Boost(-boost); //send parent to cm for calc
  
  Float_t ex = random->Uniform(minEx,maxEx);
  //Matches Maria method
  Float_t E_ej_cm = (eject.mass*eject.mass - (recoil.mass+ex)*(recoil.mass+ex) +
                     parent_LV.E()*parent_LV.E())/(2.0*parent_LV.E());
  Float_t P_ej_cm = sqrt(E_ej_cm*E_ej_cm-eject.mass*eject.mass);
  Float_t theta_ej_cm = random->Uniform(0.,TMath::Pi());
  Float_t phi_ej_cm = random->Uniform(0.,2*TMath::Pi());

  eject_LV.SetPxPyPzE(P_ej_cm*sin(theta_ej_cm)*cos(phi_ej_cm),
                      P_ej_cm*sin(theta_ej_cm)*sin(phi_ej_cm),
                      P_ej_cm*cos(theta_ej_cm),
                      E_ej_cm);
  eject_LV.Boost(boost); //send ejectile to lab frame
  //fill eject nucleus with all lab frame info
  eject.E = eject_LV.E(); eject.p = sqrt(eject.E*eject.E-eject.mass*eject.mass);
  eject.px = eject_LV.Px(); eject.py = eject_LV.Py(); eject.pz = eject_LV.Pz();
  eject.theta = eject_LV.Theta(); eject.phi = eject_LV.Phi(); eject.KE = eject.E-eject.mass;
  //Lorentz vector returns phi from -pi to pi; need 0 to 2pi
  if(eject.phi<0) {
    eject.phi += 2*TMath::Pi();
  }

  parent_LV.Boost(boost); //send parent back to lab frame
  recoil_LV = parent_LV - eject_LV;

  //fill recoil nucleus
  recoil.E = recoil_LV.E(); recoil.px = recoil_LV.Px(); recoil.py = recoil_LV.Py();
  recoil.pz = recoil_LV.Pz(); recoil.theta = recoil_LV.Theta(); recoil.phi = recoil_LV.Phi();
  recoil.p = sqrt(recoil.E*recoil.E-recoil_LV.M()*recoil_LV.M()); 
  recoil.KE = recoil.E-recoil_LV.M();
  recoil.Ex = recoil_LV.M() - recoil.mass;
  if(recoil.phi<0) {
    recoil.phi += 2*TMath::Pi();
  }
}

bool calculator::Step2_CalcBreakup() {
  Float_t theta_br1_cm = random->Uniform(0.0, TMath::Pi());
  Float_t phi_br1_cm = random->Uniform(0.0, 2.0*TMath::Pi());
 
  TVector3 boost = recoil_LV.BoostVector();
  recoil_LV.Boost(-boost); //send recoil to cm frame for calc

  Double_t E_br1_cm = (break1.mass*break1.mass-break2.mass*break2.mass+
                       recoil_LV.E()*recoil_LV.E())/(2.0*recoil_LV.E());
  Double_t P_br1_cm = sqrt(E_br1_cm*E_br1_cm-break1.mass*break1.mass);

  break1_LV.SetPxPyPzE(P_br1_cm*sin(theta_br1_cm)*cos(phi_br1_cm),
                       P_br1_cm*sin(theta_br1_cm)*sin(phi_br1_cm),
                       P_br1_cm*cos(theta_br1_cm),
                       E_br1_cm);
  break1_LV.Boost(boost);//send 1st breakup particle to lab
  
  break1.E = break1_LV.E(); break1.px = break1_LV.Px(); break1.py = break1_LV.Py();
  break1.pz = break1_LV.Pz(); break1.theta = break1_LV.Theta(); break1.phi = break1_LV.Phi();
  break1.KE = break1_LV.E()-break1.mass; 
  break1.p = sqrt(break1.E*break1.E-break1.mass*break1.mass);
  if(break1.phi<0) {
    break1.phi += 2*TMath::Pi();
  }

  recoil_LV.Boost(boost);//send recoil back to lab
  break2_LV = recoil_LV-break1_LV;

  break2.E = break2_LV.E(); break2.px = break2_LV.Px(); break2.py = break2_LV.Py();
  break2.pz = break2_LV.Pz(); break2.theta = break2_LV.Theta(); break2.phi = break2_LV.Phi();
  break2.KE = break2_LV.E()-break2.mass; 
  break2.p = sqrt(break2.E*break2.E-break2.mass*break2.mass);
  if(break2.phi<0) {
    break2.phi += 2*TMath::Pi();
  }
  return true;
}

void calculator::GetPCCoords(nucleus &nuc) {
  Float_t pcr, xpc, ypc, zpc, xw, yw;
  Float_t phi_wire = 0.0, pcr_offset = 0.0, phi_woffset = 0.0;
  Float_t phi_diff = 0.0, min_diff = 380.0;
  Int_t min_index = -1;
  if(nuc.theta == -10) return;
  for (int wid=0; wid<24; wid++) {//loop over 24 wires
    if(wid == 6) continue;
    pcr = 3.846264509;//loaded pcrs
    phi_wire = (24-wid)*2.0*TMath::Pi()/24.0+TMath::Pi()/2.0;
    if (phi_wire>2.0*TMath::Pi()) phi_wire = phi_wire-2.0*TMath::Pi();

    xw = pcr*cos(phi_wire);
    yw = pcr*sin(phi_wire);
    
    //Since beam is no longer necessarily at origin, need to offset
    Float_t dxw = xw-x_rxn;
    Float_t dyw = yw-y_rxn;
    pcr_offset = sqrt(dxw*dxw+dyw*dyw);
    phi_woffset = phi_from_cart(dxw, dyw);
    phi_diff = abs(nuc.phi-phi_woffset);
    if(phi_diff<TMath::Pi()*2.0 && phi_diff>TMath::Pi()) {
       phi_diff = TMath::Pi()*2.0-phi_diff;
    }
    if(phi_diff < min_diff) {//take which ever wire has the smallest phi_diff
      min_diff = phi_diff;
      min_index = wid;
      if (nuc.theta == TMath::Pi()/2.0) zpc = z_rxn;
      else if (nuc.theta == 0.0 || nuc.theta == TMath::Pi()) zpc = -10;
      else zpc = z_rxn - pcr_offset/tan(nuc.theta);
      zpc = random->Gaus(zpc, 0.5); //Sim PC wire resolution (alpha's roughly 1cm, p's 2cm)
      xpc = xw;
      ypc = yw;
    } 
  }
  //store pc info 
  nuc.PCxyz.push_back(xpc); nuc.PCxyz.push_back(ypc); nuc.PCxyz.push_back(zpc);
  nuc.wireID = min_index;
}

bool calculator::GetSX3Coords(nucleus &nuc) {
  Float_t x1,x2,x3, y1,y2,y3, z1,z2,z3, p,q,r, L=100.0;
  bool good_hit = false, front_hit = false;  
  if(nuc.theta == -10) return false;
  for (Int_t detid = 6; detid<30; detid++) {
    if(BadDetList.IsBad(detid)) continue;//check if det is declared as bad
    for(Int_t chan = 0; chan<8; chan++) {
      if(good_hit || (chan>3 && !front_hit)) break;//done if hit was found OR there wasnt a front
      if (BadDetList.IsBad(detid, chan)) continue; //check if chan is decleared as bad

      //Need to make a line from the rxn point to some far away location, then check if line
      //intersects with plane of detector
      p = x_rxn+(L*sin(nuc.theta)*cos(nuc.phi));
      q = y_rxn+(L*sin(nuc.theta)*sin(nuc.phi));
      r = z_rxn-(L*cos(nuc.theta)); //Assumes cos(theta) = (x_rxn-r)/L everywhere
      
      //Get 3 points to define detector plane
      x1 = S3.GetX1(detid, chan);
      y1 = S3.GetY1(detid, chan);
      z1 = S3.GetZ1(detid, chan);

      x2 = S3.GetX2(detid, chan);
      y2 = S3.GetY2(detid, chan);
      z2 = S3.GetZ1(detid, chan);

      x3 = S3.GetX1(detid, chan);
      y3 = S3.GetY1(detid, chan);
      z3 = S3.GetZ2(detid, chan);

      //Make vectors
      TVector3 plane_vec1(0.,0.,0.), plane_vec2(0.,0.,0.), det_point(0.,0.,0.),
               L_point(0.,0.,0.), rxn_point(0.,0.,0.), norm(0.,0.,0.);
      plane_vec1.SetXYZ(x2-x1,y2-y1,z2-z1);
      plane_vec2.SetXYZ(x3-x1,y3-y1,z3-z1);
      norm = plane_vec1.Cross(plane_vec2);
      det_point.SetXYZ(x1,y1,z1);
      rxn_point.SetXYZ(x_rxn, y_rxn, z_rxn);
      L_point.SetXYZ(p,q,r);

      //Equation of line Vec1 = Vec_rxn + t*(Vec_distant - Vec_rxn)
      //@ the point of intersection with the detector know that n.Dot(det_point-Vec1) = 0
      //solve for t
      Double_t t = norm.Dot(det_point-rxn_point)/norm.Dot(L_point-rxn_point);
      //Now can find x,y,z on det, as long as t > 0
      if(t > 0) {
        Float_t x_si = x_rxn + t*(p-x_rxn);
        Float_t y_si = y_rxn + t*(q-y_rxn);
        Float_t z_si = z_rxn + t*(r-z_rxn);

        TVector3 hit_vec; hit_vec.SetXYZ(x_si, y_si, z_si);
 
        //also have to make sure within bounds of detector 
        Double_t t1 = (hit_vec-det_point).Dot(plane_vec1)/(plane_vec1.Mag2());
        Double_t t2 = (hit_vec-det_point).Dot(plane_vec2)/(plane_vec2.Mag2());
        if((t1>0.0 && t1<1.0) && (t2>0.0 && t2<1.0)) {
          if(chan<4) {//fronts
            front_hit = true;
          } else {//backs
            nuc.detID = detid;
            nuc.SX3xyz.push_back(x_si); nuc.SX3xyz.push_back(y_si); nuc.SX3xyz.push_back(z_si);
            good_hit = true;
            GetPCCoords(nuc);
            return good_hit;
          }
        }
      }
    }
  }
  return good_hit;
} 

bool calculator::GetQ3Coords(nucleus &nuc) {
  Float_t L;
  Float_t x1,x2, y1,y2,p,q,r;
  Float_t phi1, phi2;
  
  bool good_hit = 0;
  bool front_hit = 0;
  if(nuc.theta == -10) return false;
  for (Int_t detid = 0; detid<4; detid++) {
    if(BadDetList.IsBad(detid)) continue;
    for(Int_t chan = 0; chan<32; chan++) {
      if(chan>16 && !front_hit) break;
      if(BadDetList.IsBad(detid, chan)) continue;
      if(nuc.theta >= TMath::Pi()/2.0 || z_rxn <= Q3.GetZ(detid)) continue;
      else {
        L = (z_rxn-Q3.GetZ(detid))/cos(nuc.theta);
        p = x_rxn + L*sin(nuc.theta)*cos(nuc.phi);
        q = y_rxn + L*sin(nuc.theta)*sin(nuc.phi);
        r = Q3.GetZ(detid);
      }
      Float_t rho = sqrt(p*p+q*q);
      if (!qqq_hit(p,q) || rho<Q3.GetRInner(detid, chan) || rho>Q3.GetROuter(detid, chan))
        continue;
      
      x1 = Q3.GetX1(detid, chan);
      y1 = Q3.GetY1(detid, chan);
      x2 = Q3.GetX2(detid, chan);
      y2 = Q3.GetY2(detid, chan);
      Float_t dx1 = x1-x_rxn;
      Float_t dy1 = y1-y_rxn;
      Float_t dx2 = x2-x_rxn;
      Float_t dy2 = y2-y_rxn;
      phi1 = phi_from_cart(dx1, dy1);
      phi2 = phi_from_cart(dx2, dy2);
      Int_t detID = qqq_phi_hit(x1,y1,nuc.phi,phi1,phi2);
      if(detID < 0) continue;
      if(chan < 16) front_hit = true;
      else {
        good_hit = true;
        nuc.Q3xyz.push_back(p); nuc.Q3xyz.push_back(q); nuc.Q3xyz.push_back(r);
        nuc.detID = detID;
        GetPCCoords(nuc);
        return good_hit;
      }
    }//chan loop
  }//end detid loop
  return good_hit;
}

track calculator::GetIntPoint(nucleus nuc, LookUp *nuc_eloss) {
  track trParticle;
  track temp;
  if(nuc.SX3xyz.size() == 3) {
    trParticle.SiR = sqrt(nuc.SX3xyz[0]*nuc.SX3xyz[0]+nuc.SX3xyz[1]*nuc.SX3xyz[1]);
    trParticle.SiZ = nuc.SX3xyz[2];
    trParticle.SiPhi = phi_from_cart(nuc.SX3xyz[0], nuc.SX3xyz[1]);
    //get particle's final KE at the si detector (SiE)
    Float_t path = sqrt(pow((nuc.SX3xyz[0]-x_rxn), 2.0)+pow((nuc.SX3xyz[1]-y_rxn), 2.0)+
                        pow((nuc.SX3xyz[2]-z_rxn), 2.0));
    trParticle.SiE = nuc_eloss->GetLookupEnergy(nuc.KE, path);
    if (trParticle.SiE < 0.6) {
      trParticle.SiE = -10.0;
      return trParticle;
    }
  } else if(nuc.Q3xyz.size() == 3) {
    trParticle.SiR = sqrt(nuc.Q3xyz[0]*nuc.Q3xyz[0]+nuc.Q3xyz[1]*nuc.Q3xyz[1]);
    trParticle.SiZ = nuc.Q3xyz[2];
    trParticle.SiPhi = phi_from_cart(nuc.Q3xyz[0], nuc.Q3xyz[1]);
    //get particle's final KE at the si detector (SiE)
    Float_t path = sqrt(pow((nuc.Q3xyz[0]-x_rxn), 2.0)+pow((nuc.Q3xyz[1]-y_rxn), 2.0)+
                        pow((nuc.Q3xyz[2]-z_rxn), 2.0));
    trParticle.SiE = nuc_eloss->GetLookupEnergy(nuc.KE, path);
    if (trParticle.SiE < 0.6) {
      trParticle.SiE = -10.0;
      return trParticle;
    }
  } else {
    cout<<"Error in calculator::GetIntPoint! Both SX3xyz and Q3xyz are unset"<<endl;
    cout<<"SX size: "<<nuc.SX3xyz.size()<<" Q3 size: "<<nuc.Q3xyz.size()<<endl;
    return trParticle;
  }
  
  if(nuc.PCxyz.size() == 3) {
    trParticle.PCR = sqrt(nuc.PCxyz[0]*nuc.PCxyz[0]+nuc.PCxyz[1]*nuc.PCxyz[1]);
    trParticle.PCZ = nuc.PCxyz[2];
  } else {
    cout<<"Error in calculator::GetIntPoint! PCxyz is unset"<<endl;
    return trParticle;
  }

  //Use y = m*x+b to find intpoint from Si position and PC position assuming that 
  //rxn occurs at x=0 and y=0 (just like analysis)
  Float_t m = (trParticle.PCR-trParticle.SiR)/(trParticle.PCZ-trParticle.SiZ);
  if(m == 0.0) {
    cout<<"Error in calculator::GetIntPoint! m = 0"<<endl;
    trParticle = temp;
    return trParticle;
  }
  Float_t b = trParticle.PCR - m*trParticle.PCZ;
  trParticle.IntPoint = -b/m;
  if(trParticle.IntPoint<0.0 || trParticle.IntPoint>ana_length) {
    trParticle = temp;
    return trParticle;
  }

  //Get the theta angle from the rcnstrct int point along with the path length to detector
  Float_t zdist = trParticle.IntPoint - trParticle.SiZ;
  if(zdist>0.0) {
    trParticle.theta = atan(trParticle.SiR/zdist);
    trParticle.path = trParticle.SiR/(sin(trParticle.theta));
  } else if(zdist<0.0) {
    trParticle.theta = TMath::Pi()+atan(trParticle.SiR/zdist);
    trParticle.path = trParticle.SiR/(sin(trParticle.theta));
  } else {
    trParticle.theta = TMath::Pi()/2.0;
    trParticle.path = trParticle.SiR;
  }

  //Get reconstruct beam energy from eloss
  Float_t intp = ana_length - trParticle.IntPoint;
  if(intp>0.0 && intp<ana_length) {
    trParticle.beamKE = beam_eloss->GetLookupEnergy(maxBeamKE, intp);
    if(trParticle.beamKE<0.0 || trParticle.beamKE>maxBeamKE) {
      cout<<"Error in calculator::GetIntPoint! Illegal beamKE"<<endl;
      trParticle = temp;
      return trParticle;
    }
  } else {
    trParticle = temp;
    return trParticle;
  }
  return trParticle;
}
