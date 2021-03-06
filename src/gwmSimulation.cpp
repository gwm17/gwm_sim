#include "gwmSimulation.h"
#include <fstream>
#include <TLorentzVector.h>


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

void simulation::resetNucleus(nucleus &nuc) {
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

void simulation::resetNucleus(track &tr) {
  tr.IntPoint = -10;
  tr.SiR = -10;
  tr.PCR = -10;
  tr.SiZ = -10;
  tr.PCZ = -10;
  tr.theta = -10;
  tr.SiPhi = -10;
  tr.path = -10;
  tr.beamKE = -10;
  tr.SiE = -10;
}

void simulation::resetNucleus(rcNucleus &rc) {
  rc.px = -10;
  rc.py = -10;
  rc.pz = -10;
  rc.p = -10;
  rc.E = -10;
  rc.KE = -10;
  rc.theta = -10;
  rc.phi = -10;
  rc.Ex = -10;
  rc.beamKE = -10;
  rc.beamPz = -10;
  rc.IntPoint = -10;
  rc.deltaPhi  = -10;
  rc.deltaTheta = -10;
  rc.mass_sq1 = -10;
  rc.mass_sq2 = -10;
  rc.ejectKE = -10;
  rc.break1KE = -10;
  rc.break2KE = -10;
}

void simulation::resetNuclei() {
  resetNucleus(beam); resetNucleus(target); resetNucleus(recoil); resetNucleus(ejectile);
  resetNucleus(break1); resetNucleus(break2); resetNucleus(trEjectile);
  resetNucleus(trBreak1); resetNucleus(trBreak2); resetNucleus(rcRecoil);
  resetNucleus(rcRecoil_qval);
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
  TLorentzVector br1_LV, br2_LV, eject_LV, a2_LV, ap_LV, temp1, temp2;
  Float_t ap_mass, a2_mass, ex1, ex2;
  int beamBin;

  for (int i=0; i<nparticles; i++) {
    cout<<"\rProgress: "<<i/blarticles*100.0<<"% "<<flush;
    resetNuclei();
    calc.ResetNucleus(calc.beam);
    calc.ResetNucleus(calc.target);
    calc.ResetNucleus(calc.recoil);
    calc.ResetNucleus(calc.eject);
    calc.ResetNucleus(calc.break1);
    calc.ResetNucleus(calc.break2);
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

    if(calc.GetSX3Coords(calc.break1)) {
      d2 = 1;
    } else if (calc.GetQ3Coords(calc.break1)) {
      d2 = 1;
    } else {
      d2 = 0;
    }

    if(calc.GetSX3Coords(calc.break2)) {
      d3 = 1;
    } else if (calc.GetQ3Coords(calc.break2)) {
      d3 = 1;
    } else {
      d3 = 0;
    }
    //Grab sim nuclei from calc
    beam = calc.beam;
    target = calc.target;
    recoil = calc.recoil;
    ejectile = calc.eject;
    break1 = calc.break1;
    break2 = calc.break2;
    
    beamBin = floor(beam.KE/.2);
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
      } else {
        ap_LV = br2_LV+br1_LV;
      }
      a2_LV = br2_LV+eject_LV;
      ap_mass = ap_LV.M()*ap_LV.M()/1e6;
      a2_mass = a2_LV.M()*a2_LV.M()/1e6;
      MyFill("DPfrxn1",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);   
      MyFill("DPftot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
      MyFill(Form("DPftot_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
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
      MyFill("DPfrxn0",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);   
      MyFill("DPftot",90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);   
      MyFill(Form("DPftot_%ibin",beamBin),90,55.55,56.00,a2_mass,50,21.75,22.00,ap_mass);
    }
 
    //if we detected, do reconstruction
    if(d1 == 1 && d2 == 1 && d3 == 1) {
      trEjectile = calc.GetIntPoint(ejectile, calc.eject_eloss);
      trBreak1 = calc.GetIntPoint(break1, calc.break1_eloss);
      trBreak2 = calc.GetIntPoint(break2, calc.break2_eloss);
      if(trEjectile.IntPoint != -10 && trBreak1.IntPoint != -10 && trBreak2.IntPoint != -10) {
        passed++;
        detectFlag = 1;
        Reconstruct rc(beam, target, ejectile, recoil, break1, break2);
        rc.ELoss_eject = calc.eject_eloss; rc.ELoss_break1 = calc.break1_eloss;
        rc.ELoss_break2 = calc.break2_eloss;
        rc.CalcRecoil_MultiParticle(trEjectile, trBreak1, trBreak2, rcRecoil_qval, rcRecoil);
        rc.ELoss_eject =NULL; rc.ELoss_break1 = NULL; rc.ELoss_break2 = NULL;

        if(rxnTag == 1) {
          MyFill("DPdrxn1",90,55.55,56.00,rcRecoil.mass_sq1,50,21.75,22.00,rcRecoil.mass_sq2);
          MyFill("DPdtot",90,55.55,56.00,rcRecoil.mass_sq1,50,21.75,22.00,rcRecoil.mass_sq2);
          MyFill(Form("DPdtot_%ibin",beamBin),90,55.55,56.0,rcRecoil.mass_sq1,
                                              50,21.75,22.0,rcRecoil.mass_sq2);
        } else if (rxnTag == 0) {
          MyFill("DPdrxn0",90,55.55,56.00,rcRecoil.mass_sq1,50,21.75,22.00,rcRecoil.mass_sq2);
          MyFill("DPdtot",90,55.55,56.00,rcRecoil.mass_sq1,50,21.75,22.00,rcRecoil.mass_sq2);
          MyFill(Form("DPdtot_%ibin",beamBin),90,55.55,56.0,rcRecoil.mass_sq1,
                                              50,21.75,22.0,rcRecoil.mass_sq2);
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
  for(int i =0; i<nparticles; i++) {
    cout<<"\rProgress: "<<i/blarticles*100<<"% "<<flush;
    resetNuclei();
    calc.ResetNucleus(calc.beam);
    calc.ResetNucleus(calc.target);
    calc.ResetNucleus(calc.recoil);
    calc.ResetNucleus(calc.eject);
    calc.ResetNucleus(calc.break1);
    calc.ResetNucleus(calc.break2);
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
    if(d==1) {
      trEjectile = calc.GetIntPoint(ejectile, calc.eject_eloss);
      if(trEjectile.IntPoint != -10) {
        detectFlag = 1;
        passed++;
        Reconstruct rc(beam, target, ejectile, recoil, break1, break2);
        rc.ELoss_eject = calc.eject_eloss; rc.ELoss_break1 = NULL;
        rc.ELoss_break2 = NULL; 
        rc.CalcRecoil(trEjectile, rcRecoil);
        Double_t Q = beam.mass+target.mass-ejectile.mass-recoil.mass;
        rc.CalcRecoil_from_Qvalue(Q, trEjectile, rcRecoil_qval);
        rc.ELoss_eject = NULL;
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
    TH2F *ratio = new TH2F(Form("DPeff_%ibin",i),Form("DPeff_%ibin",i), 90,55.55,56.0,50,21.75,22.00);
    ratio->Divide(detect,free);
    rootObj->Add(ratio);
    Double_t sum = 0;
    Double_t count = 0;
    for(int j = 1; j<ratio->GetNbinsX()+1; j++) {
      for(int k = 1; k<ratio->GetNbinsY()+1; k++) {
        Int_t bin = ratio->GetBin(j,k);
        Double_t value = ratio->GetBinContent(bin);
        if(value > 0 && value<1.0) {
          sum += value;
          count++;
        }
      } 
    } 
    eff.push_back(sum/count);
    bke.push_back(i*0.2);
  }
  beff = new TGraph(eff.size(), &(bke[0]), &(eff[0]));
  rootObj->Add(beff);
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

  /*outTree = new TTree("simTree","simTree");
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
  outTree->Branch("detectFlag", &detectFlag, "detectFlag/I");*/
  for(unsigned int i=0; i<mb.size(); i++) {
    rxnTag = i;
    Int_t passed = 0;
    cout<<"Simulating reaction "<<rxnTag<<" ..."<<endl;
    calculator calc(mb[i],mt[i],mr[i],me[i],mbr1[i],mbr2[i],rxMin[i],rxMax[i]);
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
            Float_t eff = ((Float_t) passed)/((Float_t) nBeamParticles);
            cout<<"Efficiency: "<<eff<<endl;
            break;
          }
      case 2:
          {//7Be(d,3He)6Li
            calc.eject_eloss = he3_eloss;
            passed = Go1StepCalc(calc, nBeamParticles);
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
  //outTree->Write(outTree->GetName(), TObject::kOverwrite);
  rootObj->Write();
  rootObj->Clear();
  outFile->Close();
}
