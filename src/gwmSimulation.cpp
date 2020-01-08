#include "gwmSimulation.h"
#include <fstream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "nucFuncs.h"


using namespace std;

simulation::simulation() {

  be7eloss_name = "./srim/7be_in_d2_290torr.eloss";
  peloss_name = "./srim/p_in_d2_400torr.eloss";
  he4eloss_name = "./srim/4he_in_d2_400torr.eloss";
  he3eloss_name = "./srim/4he_in_d2_400torr.eloss";
}

simulation::~simulation() {
  delete be7_eloss;
  delete he4_eloss;
  delete he3_eloss;
  delete p_eloss;
  delete rootObj;
  be8file->Close();
  li5file->Close();
}

void simulation::MyFill(string name, int binsX, double lowX, double highX, double valueX) {
  TH1F *histo = (TH1F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) {
    histo->Fill(valueX);
  } else {
    TH1F *h = new TH1F(name.c_str(), name.c_str(), binsX, lowX, highX);
    h->Fill(valueX);
    rootObj->Add(h);
  }
}

void simulation::MyFill(string name, int binsX, double lowX, double highX, double valueX,
                                     int binsY, double lowY, double highY, double valueY) {
  TH2F *histo = (TH2F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) { 
    histo->Fill(valueX, valueY);
  } else {
    TH2F* h = new TH2F(name.c_str(),name.c_str(), binsX, lowX, highX, binsY, lowY, highY);
    h->Fill(valueX, valueY);
    rootObj->Add(h);
  }
}

void simulation::getCut() {
  be8file = new TFile("./cuts/be8Cut.root","READ");
  be8Cut = (TCutG*) be8file->Get("CUTG");
  be8Cut->SetName("be8Cut");
  be8Cut->SetVarX("recoil.a2_mass_sq");
  be8Cut->SetVarY("recoil.ap_mass_sq");
  rootObj->Add(be8Cut);
  li5file = new TFile("./cuts/li5Cut.root", "READ");
  li5Cut = (TCutG*) li5file->Get("CUTG");
  li5Cut->SetName("li5Cut");
  li5Cut->SetVarX("recoil.a2_mass_sq");
  li5Cut->SetVarY("recoil.ap_mass_sq");
  rootObj->Add(li5Cut);
}

void simulation::resetNuclei() {
  resetNucleus(beam,1); resetNucleus(target,1); resetNucleus(recoil,1);
  resetNucleus(ejectile,1);resetNucleus(break1,1); resetNucleus(break2,1);
  resetTrack(trEjectile); resetTrack(trBreak1); resetTrack(trBreak2);
  resetRcNucleus(rcRecoil); resetRcNucleus(rcRecoil_qval);
}

void simulation::GetReactions(string listName) {

  ifstream rxnList;
  rxnList.open(listName);
  Float_t m_targ, m_beam, m_recoil, m_eject, m_br1, m_br2, rexMax, rexMin;
  if(rxnList.is_open()) {
    cout<<"Getting reactions from "<<listName<<" ..."<<endl;
    string junk;
    rxnList>>junk>>minBeamKE>>junk>>maxBeamKE>>junk>>nBeamParticles;
    rxnList>>junk>>junk>>junk>>junk>>junk>>junk>>junk>>junk>>junk;
    while(rxnList>>junk) {
      rxnList>>m_targ;
      rxnList>>m_beam;
      rxnList>>m_recoil;
      rxnList>>m_eject;
      rxnList>>m_br1;
      rxnList>>m_br2;
      rxnList>>rexMin;
      rxnList>>rexMax;

      mt.push_back(m_targ);
      mb.push_back(m_beam);
      mr.push_back(m_recoil);
      me.push_back(m_eject);
      mbr1.push_back(m_br1);
      mbr2.push_back(m_br2);
      rxMax.push_back(rexMax);
      rxMin.push_back(rexMin);
    }
  } else {
    cout<<"Could not open reaction list!"<<endl;
    return;
  }
}

void simulation::initEloss() {
  be7_eloss = new LookUp(be7eloss_name, m_7be);
  he4_eloss = new LookUp(he4eloss_name, m_alpha);
  he3_eloss = new LookUp(he3eloss_name, m_3he);
  p_eloss = new LookUp(peloss_name, m_p);

  be7_eloss->InitializeLookupTables(30.0,200.0,0.01,0.04);
  he4_eloss->InitializeLookupTables(30.0, 900.0,0.01,0.04);
  he3_eloss->InitializeLookupTables(30.0,1200.0,0.01,0.04);
  p_eloss->InitializeLookupTables(30.0,11900.0,0.01,0.04);
}

Int_t simulation::Go2StepCalc(calculator &calc, Int_t nparticles) {
  Int_t passed = 0;
  float blarticles = nparticles;
  TLorentzVector br1_LV, br2_LV, eject_LV, a2_LV, ap_LV, temp1, temp2, beam_LV, target_LV, parent_LV;
  Float_t ap_mass, a2_mass, ex1, ex2;
  int beamBin;
  Reconstruct rc;
  for (int i=0; i<nparticles; i++) {
    cout<<"\rProgress: "<<i/blarticles*100.0<<"% "<<flush;
    resetNuclei();
    calc.ResetCalc();
    rc.ResetRC();
    detectFlag = 0;
    calc.InitializeBeam(minBeamKE, maxBeamKE);
    x_rxn = calc.x_rxn;
    y_rxn = calc.y_rxn;
    z_rxn = calc.z_rxn;
    calc.Step1_CalcRxn();
    if(!calc.Step2_CalcBreakup()) {
     continue;
    }
    
    int d1 = 0, d2 = 0, d3 = 0;
    //See if we actually detect all three particles
    if(calc.GetSX3Coords(calc.eject)) {
      d1 = 1;
    } else if (calc.GetQ3Coords(calc.eject)) {
      d1 = 1;
    } else {
      d1 = 0;
    }

    if(d1) {
      if(calc.GetSX3Coords(calc.break1)) {
        d2 = 1;
      } else if (calc.GetQ3Coords(calc.break1)) {
        d2 = 1;
      } else {
        d2 = 0;
      }
    }

    if(d1 && d2) {
      if(calc.GetSX3Coords(calc.break2)) {
        d3 = 1;
      } else if (calc.GetQ3Coords(calc.break2)) {
        d3 = 1;
      } else {
        d3 = 0;
      }
    }
    //Grab sim nuclei from calc for dalitz
    beam = calc.beam;
    target = calc.target;
    recoil = calc.recoil;
    ejectile = calc.eject;
    break1 = calc.break1;
    break2 = calc.break2;
    beamBin = ceil(beam.KE/.2);
    beam_LV.SetPxPyPzE(0., 0., beam.pz, beam.E);
    target_LV.SetPxPyPzE(0., 0., 0., target.mass);
    parent_LV = beam_LV+target_LV;
    TVector3 boost = parent_LV.BoostVector();
    parent_LV.Boost(-boost);
    Double_t Ecm = parent_LV.E()-beam.mass-target.mass;
    recoil.Ecm = Ecm;
    Double_t a1_theta_cm = -10;
    Double_t p_theta_cm = -10;
    Int_t a1_theta_bin = -10;
    Int_t p_theta_bin = -10;
    if(rxnTag == 1) {
      br1_LV.SetPxPyPzE(break1.px, break1.py, break1.pz, break1.E);
      br2_LV.SetPxPyPzE(break2.px, break2.py, break2.pz, break2.E);
      eject_LV.SetPxPyPzE(ejectile.px, ejectile.py, ejectile.pz, ejectile.E);
      temp1 = br1_LV+eject_LV;
      temp2 = br1_LV+br2_LV;
      ex1 = temp1.M()-4667.6163636366931;
      ex2 = temp2.M()-4667.6163636366931;
      if(ex1<ex2) {
        ap_LV = br1_LV+eject_LV;
        br2_LV.Boost(-boost);
        a1_theta_cm = br2_LV.Theta();
        recoil.eject_theta_cm = a1_theta_cm;
        br2_LV.Boost(boost);
      } else {
        ap_LV = br2_LV+br1_LV;
        eject_LV.Boost(-boost);
        a1_theta_cm = eject_LV.Theta();
        recoil.eject_theta_cm = a1_theta_cm;
        eject_LV.Boost(boost);
      }
      a1_theta_bin = ceil(cos(a1_theta_cm)/0.2)+5;
      a2_LV = br2_LV+eject_LV;
      ap_mass = ap_LV.M()*ap_LV.M()/1e6;
      a2_mass = a2_LV.M()*a2_LV.M()/1e6;
      recoil.a2_mass_sq = a2_mass;
      recoil.ap_mass_sq = ap_mass;
      MyFill("DPfrxn1",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);   
      MyFill("DPftot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPfrxn1_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPftot_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPftot_%ibbin_%itbin",beamBin,a1_theta_bin),90,55.55,56.00,a2_mass,
                                                                50,21.75,22.00,ap_mass);
      if(Ecm > 0.267 && Ecm < 0.468 && li5Cut->IsInside(a2_mass,ap_mass)) {
        MyFill("AngDistf_5Li_0.267_0.468",10,-1.0,1.0,cos(a1_theta_cm));
      } else if (Ecm > 1.11 && Ecm < 1.2919 && li5Cut->IsInside(a2_mass,ap_mass)) {
        MyFill("AngDistf_5Li_1.11_1.2919",10,-1.0,1.0,cos(a1_theta_cm));
      }
    } else if (rxnTag == 0) {
      br1_LV.SetPxPyPzE(break1.px, break1.py, break1.pz, break1.E);
      br2_LV.SetPxPyPzE(break2.px, break2.py, break2.pz, break2.E);
      eject_LV.SetPxPyPzE(ejectile.px, ejectile.py, ejectile.pz, ejectile.E);
      a2_LV = br1_LV+br2_LV;
      temp1 = br1_LV+eject_LV;
      temp2 = br2_LV+eject_LV;
      ex1 = temp1.M()-4667.6163636366931;
      ex2 = temp2.M()-4667.6163636366931;
      if(ex1<ex2) {
        ap_LV = br1_LV+eject_LV;
      } else {
        ap_LV = br2_LV+eject_LV;
      }
      a2_mass = a2_LV.M()*a2_LV.M()/1e6;
      ap_mass = ap_LV.M()*ap_LV.M()/1e6;
      recoil.a2_mass_sq = a2_mass;
      recoil.ap_mass_sq = ap_mass;
      MyFill("DPfrxn0",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);   
      MyFill("DPftot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      eject_LV.Boost(-boost);
      p_theta_cm = eject_LV.Theta();
      recoil.eject_theta_cm = p_theta_cm;
      eject_LV.Boost(boost);
      p_theta_bin = ceil(cos(p_theta_cm)/0.2)+5;
      MyFill(Form("DPfrxn0_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPftot_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPftot_%ibbin_%itbin",beamBin,p_theta_bin),90,55.55,56.00,a2_mass,
                                                              50,21.75,22.00,ap_mass);
      if(Ecm>1.07 && Ecm < 1.244 && be8Cut->IsInside(a2_mass, ap_mass)) {
        MyFill("AngDistf_8Be_1.07_1.244",10,-1.0,1.0,cos(p_theta_cm));
      }
    }
    
    if(rxnTag == 1 && Ecm>1.11 && Ecm<1.2919 && li5Cut->IsInside(a2_mass,ap_mass)) {
      if (d1) {
        MyFill("AngDist_eject_5Li_1.11_1.2919",20,-1.0,1.0,cos(a1_theta_cm));
      } if(d2) {
        MyFill("AngDist_br1_5Li_1.11_1.2919",20,-1.0,1.0,cos(a1_theta_cm));
      }  if(d3) {
        MyFill("AngDist_br1_5Li_1.11_1.2919",20,-1.0,1.0,cos(a1_theta_cm));
      }
    }

    //if we detected, do reconstruction
    if(d1 == 1 && d2 == 1 && d3 == 1 && ejectile.wireID != break1.wireID &&
       ejectile.wireID != break2.wireID && break1.wireID != break2.wireID) {
      trEjectile = calc.GetIntPoint(ejectile, calc.eject_eloss);
      trBreak1 = calc.GetIntPoint(break1, calc.break1_eloss);
      trBreak2 = calc.GetIntPoint(break2, calc.break2_eloss);
      if(trEjectile.IntPoint != -10 && trBreak1.IntPoint != -10 && trBreak2.IntPoint != -10) {
        passed++;
        detectFlag = 1;
        rc.SetParams(beam, target, ejectile, recoil, break1, break2);
        rc.ELoss_eject = calc.eject_eloss; rc.ELoss_break1 = calc.break1_eloss;
        rc.ELoss_break2 = calc.break2_eloss;
        rc.CalcRecoil_MultiParticle(trEjectile, trBreak1, trBreak2, rcRecoil_qval, rcRecoil);
        rc.ELoss_eject =NULL; rc.ELoss_break1 = NULL; rc.ELoss_break2 = NULL;

        if(rxnTag == 1) {
          MyFill("DPdrxn1",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill("DPdtot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill(Form("DPdrxn1_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill(Form("DPdtot_%ibin",beamBin),90,55.55,56.0,a2_mass,50,21.75,22.0,ap_mass);
          MyFill(Form("DPdtot_%ibbin_%itbin",beamBin,a1_theta_bin),90,55.55,56.00,a2_mass,
                                                                   50,21.75,22.00,ap_mass);
          if(Ecm > 0.267 && Ecm < 0.468 && li5Cut->IsInside(a2_mass,ap_mass)) {
            MyFill("AngDistd_5Li_0.267_0.468",10,-1.0,1.0,cos(a1_theta_cm));
          } else if (Ecm > 1.11 && Ecm < 1.2919 && li5Cut->IsInside(a2_mass,ap_mass)) {
            MyFill("AngDistd_5Li_1.11_1.2919",10,-1.0,1.0,cos(a1_theta_cm));
          }
        } else if (rxnTag == 0) {
          MyFill("DPdrxn0",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill("DPdtot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill(Form("DPdrxn0_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
          MyFill(Form("DPdtot_%ibin",beamBin),90,55.55,56.0,a2_mass,50,21.75,22.0,ap_mass);
          MyFill(Form("DPdtot_%ibbin_%itbin",beamBin,p_theta_bin),90,55.55,56.00,a2_mass,
                                                                  50,21.75,22.00,ap_mass);
          if(Ecm>1.07 && Ecm < 1.244 && be8Cut->IsInside(a2_mass, ap_mass)) {
            MyFill("AngDistd_8Be_1.07_1.244",10,-1.0,1.0,cos(p_theta_cm));
          }
        }
      }
    }
    //outTree->Fill();
  }
  cout<<endl;
  return passed;
}

Int_t simulation::Go1StepCalc(calculator &calc, Int_t nparticles) {
  Int_t passed = 0;
  float blarticles = (float)nparticles;
  Reconstruct rc;
  TLorentzVector beam_LV, parent_LV, target_LV, eject_LV;
  for(int i =0; i<nparticles; i++) {
    cout<<"\rProgress: "<<i/blarticles*100<<"% "<<flush;
    resetNuclei();
    rc.ResetRC();
    calc.ResetCalc();
    detectFlag=0;
    calc.InitializeBeam(minBeamKE, maxBeamKE);
    x_rxn = calc.x_rxn;
    y_rxn = calc.y_rxn;
    z_rxn = calc.z_rxn;
    calc.Step1_CalcRxn();
    int d;
    if(calc.GetSX3Coords(calc.eject)){
      d = 1;
    } else if (calc.GetQ3Coords(calc.eject)){
      d = 1;
    } else {
      d = 0;
    }
    beam = calc.beam;
    target = calc.target;
    recoil = calc.recoil;
    ejectile = calc.eject;
   
    Int_t beamBin = ceil(beam.KE/0.2); 
    beam_LV.SetPxPyPzE(0.0,0.0,beam.pz,beam.E);
    target_LV.SetPxPyPzE(0.,0.,0.,target.mass);
    parent_LV = beam_LV+target_LV;
    TVector3 boost = parent_LV.BoostVector();
    parent_LV.Boost(-boost);
    Double_t Ecm = parent_LV.E()-target.mass-beam.mass;
    recoil.Ecm = Ecm;
    eject_LV.SetPxPyPzE(ejectile.px, ejectile.py, ejectile.pz, ejectile.E);
    eject_LV.Boost(-boost);
    Double_t e_theta_cm = eject_LV.Theta();
    MyFill(Form("ExLi6f_%ibin",beamBin),100,-2.0,2.0,recoil.Ex);
    MyFill("AngDistf_6Li",10,-1.0,1.0,e_theta_cm);
    if(d==1) {
      trEjectile = calc.GetIntPoint(ejectile, calc.eject_eloss);
      if(trEjectile.IntPoint != -10) {
        detectFlag = 1;
        passed++;
        rc.SetParams(beam, target, ejectile, recoil, break1, break2);
        rc.ELoss_eject = calc.eject_eloss; rc.ELoss_break1 = NULL;
        rc.ELoss_break2 = NULL; 
        rc.CalcRecoil(trEjectile, rcRecoil);
        Double_t Q = beam.mass+target.mass-ejectile.mass-recoil.mass;
        rc.CalcRecoil_from_Qvalue(Q, trEjectile, rcRecoil_qval);
        rc.ELoss_eject = NULL;
        MyFill(Form("ExLi6d_%ibin",beamBin),100,-2.0,2.0,recoil.Ex);
        MyFill("AngDistd_6Li",10,-1.0,1.0,e_theta_cm);
      }
    }
    //outTree->Fill();
  }
  cout<<endl;
  return passed;
}

void simulation::getDalitzEff() {
  vector<Double_t> bke;
  vector<Double_t> eff;
  TGraph *beff;
  for(int i=1; i<51; i++) {
    TH2F *free = (TH2F*) rootObj->FindObject(Form("DPftot_%ibin",i));
    TH2F *detect = (TH2F*) rootObj->FindObject(Form("DPdtot_%ibin",i));
    TH2F *freer0 = (TH2F*) rootObj->FindObject(Form("DPfrxn0_%ibin",i));
    TH2F *freer1 = (TH2F*) rootObj->FindObject(Form("DPfrxn1_%ibin",i));
    TH2F *detectr1 = (TH2F*) rootObj->FindObject(Form("DPdrxn1_%ibin",i));
    TH2F *detectr0 = (TH2F*) rootObj->FindObject(Form("DPdrxn0_%ibin",i));
    TH2F *ratio = (TH2F*) detect->Clone(Form("DPeff_%ibin",i));
    TH2F *ratior1 = (TH2F*) detectr1->Clone(Form("DPeffrxn1_%ibin", i));
    TH2F *ratior0 = (TH2F*) detectr0->Clone(Form("DPeffrxn0_%ibin", i));
    ratio->Divide(free);
    ratior0->Divide(freer0);
    ratior1->Divide(freer1);
    rootObj->Add(ratio);
    rootObj->Add(ratior0);
    rootObj->Add(ratior1);
    Double_t sum = 0;
    Double_t count = 0;
    for(int j = 1; j<ratio->GetNbinsX()+1; j++) {
      for(int k = 1; k<ratio->GetNbinsY()+1; k++) {
        Double_t y_pos = ratio->GetYaxis()->GetBinLowEdge(k);
        if(y_pos>21.77 && y_pos<21.8) {
          Int_t bin = ratio->GetBin(j,k);
          Double_t value = ratio->GetBinContent(bin);
          if(value > 0 && value<1.0) {
            sum += value;
            count++;
          }
        }
      } 
    } 
    eff.push_back(sum/count);
    bke.push_back(i*0.2);
  }
  beff = new TGraph(eff.size(), &(bke[0]), &(eff[0]));
  rootObj->Add(beff);
}

void simulation::getAngDistEff() {
  TH1F *free1 = (TH1F*) rootObj->FindObject("AngDistf_5Li_0.267_0.468");
  TH1F *det1 = (TH1F*) rootObj->FindObject("AngDistd_5Li_0.267_0.468");
  TH1F *ratio1 = (TH1F*) det1->Clone("AngDiste_5Li_0.267_0.468");
  ratio1->Divide(free1);

  TH1F *free2 = (TH1F*) rootObj->FindObject("AngDistf_5Li_1.11_1.2919");
  TH1F *det2 = (TH1F*) rootObj->FindObject("AngDistd_5Li_1.11_1.2919");
  TH1F *ratio2 = (TH1F*) det2->Clone("AngDiste_5Li_0.1.11_1.2919");
  ratio2->Divide(free2);
 
  TH1F *free3 = (TH1F*) rootObj->FindObject("AngDistf_8Be_1.07_1.244"); 
  TH1F *det3 = (TH1F*) rootObj->FindObject("AngDistd_8Be_1.07_1.244"); 
  TH1F *ratio3 = (TH1F*) det3->Clone("AngDiste_8Be_1.07_1.244");
  ratio3->Divide(free3);

  TH1F *free4 = (TH1F*) rootObj->FindObject("AngDistf_6Li");
  TH1F *det4 = (TH1F*) rootObj->FindObject("AngDistd_6Li");
  TH1F *ratio4 = (TH1F*) det4->Clone("AngDiste_6Li");
  ratio4->Divide(free4);

  rootObj->Add(ratio1);
  rootObj->Add(ratio2);
  rootObj->Add(ratio3);
  rootObj->Add(ratio4);
}

void simulation::getSingleChanEff() {
  for(int i=1; i<51; i++) {
    TH1F *free = (TH1F*) rootObj->FindObject(Form("ExLi6f_%ibin",i));
    TH1F *det = (TH1F*) rootObj->FindObject(Form("ExLi6d_%ibin",i));
    if(free != NULL && det != NULL) {
      TH1F *ratio = (TH1F*) det->Clone(Form("ExLi6e_%ibin",i));
      ratio->Divide(free);
      rootObj->Add(ratio);
    }
  }
}

void simulation::getMixedEff() {
  for(int i=1; i<51; i++) {
    for(int j=0; j<20; j++) {
      TH2F *free = (TH2F*) rootObj->FindObject(Form("DPftot_%ibbin_%itbin",i,j));
      TH2F *det = (TH2F*) rootObj->FindObject(Form("DPdtot_%ibbin_%itbin",i,j));
      if(det == NULL) continue;
      TH2F *ratio = (TH2F*) det->Clone(Form("DPeff_%ibbin_%itbin",i,j));
      ratio->Divide(free);
      rootObj->Add(ratio);
    }
  }
}

void simulation::run(char *outName) {
  cout<<"Initializing eloss files... "<<endl;
  initEloss();

  string listName = "rxnList.txt";
  GetReactions(listName);
  
  cout<<"Output file: "<<outName<<endl;
  TFile *outFile = new TFile(outName, "RECREATE");
  rootObj = new THashTable();
  rootObj->SetOwner(false);

  outTree = new TTree("simTree","simTree");
  outTree->Branch("beam", &beam);
  outTree->Branch("target", &target);
  outTree->Branch("ejectile",&ejectile);
  outTree->Branch("recoil", &recoil);
  outTree->Branch("break1", &break1);
  outTree->Branch("break2", &break2);
  outTree->Branch("trEjectile", &trEjectile);
  outTree->Branch("trBreak1", &trBreak1);
  outTree->Branch("trBreak2", &trBreak2);
  outTree->Branch("rcRecoil", &rcRecoil);
  outTree->Branch("rcRecoil_qval", &rcRecoil_qval);
  outTree->Branch("x_rxn", &x_rxn, "x_rxn/F");
  outTree->Branch("y_rxn",&y_rxn, "y_rxn/F");
  outTree->Branch("z_rxn", &z_rxn, "z_rxn/F");
  outTree->Branch("rxnTag", &rxnTag, "rxnTag/I");
  outTree->Branch("detectFlag", &detectFlag, "detectFlag/I");
  getCut();

  calculator calc;
  for(unsigned int i=0; i<mb.size(); i++) {
    rxnTag = i;
    Int_t passed = 0;
    cout<<"Simulating reaction "<<rxnTag<<" ..."<<endl;
    calc.ResetCalc();
    calc.SetParameters(mb[i],mt[i],mr[i],me[i],mbr1[i],mbr2[i],rxMin[i],rxMax[i]);
    calc.beam_eloss = be7_eloss;
    switch(rxnTag) {
      case 0: 
          {//7Be(d,p)->2a 8Be
            calc.eject_eloss = p_eloss;
            calc.break1_eloss = he4_eloss;
            calc.break2_eloss = he4_eloss;
            passed = Go2StepCalc(calc, nBeamParticles);
            Float_t eff = ((Float_t) passed)/((Float_t) nBeamParticles);
            cout<<"Efficiency: "<<eff<<endl;
            break;
          }
      case 1:
          {//7Be(d,a)->p+a 0.0 MeV 5Li
            calc.eject_eloss = he4_eloss;
            calc.break1_eloss = p_eloss;
            calc.break2_eloss = he4_eloss;
            passed = Go2StepCalc(calc, nBeamParticles);
            getDalitzEff();
            getMixedEff();
            Float_t eff = ((Float_t) passed)/((Float_t) nBeamParticles);
            cout<<"Efficiency: "<<eff<<endl;
            break;
          }
      case 2:
          {//7Be(d,3He)6Li
            calc.eject_eloss = he3_eloss;
            passed = Go1StepCalc(calc, nBeamParticles);
            getSingleChanEff();
            getAngDistEff();
            Float_t eff = ((Float_t) passed)/((Float_t) nBeamParticles);
            cout<<"Efficiency: "<<eff<<endl;
            break;
          }
      default:
          {
            cout<<"Unexpected reaction! Check rxnTag switch in simulation::run() and add ";
            cout<<"a case."<<endl;
          }
    }
    calc.beam_eloss = NULL;
    calc.eject_eloss = NULL;
    calc.break1_eloss = NULL;
    calc.break2_eloss = NULL;
  }
  outFile->cd();
  //outTree->Write(outTree->GetName(), TObject::kOverwrite);
  rootObj->Write();
  rootObj->Clear();
  outFile->Close();
}
