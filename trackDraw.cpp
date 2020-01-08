#include <TROOT.h>
#include <TFile.h>
#include <TApplication.h>
#include <vector>
#include <TTree.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TCutG.h>
#include <TMath.h>
#include "nucleus.h"
#include "SX3.h"
#include "QQQ.h"
#include "BadDetectorList.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
  TApplication *app = new TApplication("app", &argc, argv);
  TFile *data = new TFile("rootfiles/test48.root");
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Draw();
  TTree *tree = (TTree*) data->Get("simTree");
  Int_t nentries = tree->GetEntries();
  TCutG *be8Cut = (TCutG*) data->Get("li5Cut");

  nucleus *eject=0, *recoil=0;
  Int_t df, rt;
  Float_t x_rxn, y_rxn, z_rxn;
  SX3 S3; QQQ Q3;
  BadDetectorList blist;

  tree->SetBranchAddress("ejectile",&eject);
  tree->SetBranchAddress("recoil",&recoil);
  tree->SetBranchAddress("detectFlag",&df);
  tree->SetBranchAddress("rxnTag",&rt);
  tree->SetBranchAddress("x_rxn",&x_rxn);
  tree->SetBranchAddress("y_rxn",&y_rxn);
  tree->SetBranchAddress("z_rxn",&z_rxn);

  vector<Float_t> xpnts, ypnts, zpnts;
  TPolyLine3D *baseline = new TPolyLine3D(2);
  TPolyLine3D *Q3line = new TPolyLine3D(2);
  baseline->SetPoint(0,0,0,55.45);
  baseline->SetPoint(1,0,0,0);
  Q3line->SetPoint(0,-9.9,0,0);
  Q3line->SetPoint(1,9.9,0,0);
  baseline->Draw();
  Q3line->Draw("same");
  vector<TPolyLine3D*> lines;
  int lineCount=0;
  Double_t L = 15;
  for(Int_t i=0; i<nentries; i++) {
    tree->GetEntry(i);
    xpnts.resize(0);
    ypnts.resize(0);
    zpnts.resize(0);
    if(rt == 1 && be8Cut->IsInside(recoil->a2_mass_sq, recoil->ap_mass_sq) && 
       recoil->Ecm>1.11 && recoil->Ecm<1.2919 && cos(eject->theta_cm)>-.2 &&
       cos(eject->theta_cm<0.0)) {
      lines.push_back(new TPolyLine3D(2));
      if(eject->detID < 4 && df) {
        xpnts.push_back(eject->Q3xyz[0]);
        ypnts.push_back(eject->Q3xyz[1]);
        zpnts.push_back(eject->Q3xyz[2]);
      } else if(df) {
        xpnts.push_back(eject->SX3xyz[0]);
        ypnts.push_back(eject->SX3xyz[1]);
        zpnts.push_back(eject->SX3xyz[2]);
      } else {
        Double_t p = x_rxn+(L*sin(eject->theta)*cos(eject->phi));
        Double_t q = y_rxn+(L*sin(eject->theta)*sin(eject->phi));
        Double_t r = z_rxn+(L*cos(eject->theta));
        xpnts.push_back(p);
        ypnts.push_back(q);
        zpnts.push_back(r);
      }
      xpnts.push_back(x_rxn);
      ypnts.push_back(y_rxn);
      zpnts.push_back(z_rxn);
      lines[lineCount]->SetPoint(0, xpnts[0], ypnts[0], zpnts[0]);
      lines[lineCount]->SetPoint(1, xpnts[1], ypnts[1], zpnts[1]);
      if(df) lines[lineCount]->SetLineColor(kRed);
      else if(!df) lines[lineCount]->SetLineColor(kBlue);
      lines[lineCount]->Draw("same");
      lineCount++;
      if(lineCount > 100) break;
    }
  }
  for(int i=0; i<30; i++) {
    if(i>5 && !blist.IsBad(i)) {
      Double_t x1, x2, x3, x4;
      Double_t y1, y2, y3, y4;
      Double_t z1, z2, z3, z4;
      x1 = S3.GetX1(i);
      y1 = S3.GetY1(i);
      z1 = S3.GetZ1(i);

      x2 = S3.GetX1(i);
      y2 = S3.GetY1(i);
      z2 = S3.GetZ2(i);

      x3 = S3.GetX2(i);
      y3 = S3.GetY2(i);
      z3 = S3.GetZ2(i);

      x4 = S3.GetX2(i);
      y4 = S3.GetY2(i);
      z4 = S3.GetZ1(i);

      TPolyLine3D *l = new TPolyLine3D(5);
      l->SetPoint(0,x1,y1,z1); l->SetPoint(1,x2,y2,z2); l->SetPoint(2,x3,y3,z3);
      l->SetPoint(3,x4,y4,z4); l->SetPoint(4,x1,y1,z1);
      lines.push_back(l);
      l->Draw();
    }
  }
  while(c1->WaitPrimitive()) {};
  delete baseline;
  delete Q3line;
  for(unsigned int i=0; i<lines.size(); i++) {
    delete lines[i];
  }
  cout<<"Number of tracks: "<<lineCount<<endl;
  data->Close();
}
