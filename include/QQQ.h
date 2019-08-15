#ifndef __QQQ_h
#define __QQQ_h

#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TString.h>
#include <TMath.h>

class QQQ {

 public:

  QQQ();

  Float_t GetX1(Int_t detid);
  Float_t GetY1(Int_t detid);
  Float_t GetX2(Int_t detid);
  Float_t GetY2(Int_t detid);
  Float_t GetZ(Int_t detid);
  Float_t GetTh1(Int_t detid);
  Float_t GetPh1(Int_t detid);
  Float_t GetTh2(Int_t detid);
  Float_t GetPh2(Int_t detid);
  Float_t GetRInner(Int_t detid);
  Float_t GetROuter(Int_t detid);

  Float_t GetX1(Int_t detid, Int_t chnum);
  Float_t GetY1(Int_t detid, Int_t chnum);
  Float_t GetX2(Int_t detid, Int_t chnum);
  Float_t GetY2(Int_t detid, Int_t chnum);
  Float_t GetZ(Int_t detid, Int_t chnum);
  Float_t GetTh1(Int_t detid, Int_t chnum);
  Float_t GetPh1(Int_t detid, Int_t chnum);
  Float_t GetTh2(Int_t detid, Int_t chnum);
  Float_t GetPh2(Int_t detid, Int_t chnum);
  Float_t GetRInner(Int_t detid, Int_t chnum);
  Float_t GetROuter(Int_t detid, Int_t chnum);

 private:

  //radii (cm)
  const Float_t R_INNER,
                R_OUTER;
  //***this is CYLINDRICAL R, perp. to beam axis; i.e. R^2 = X^2 + Y^2 (NO Z^2)

  Float_t ringwidth; //cm
  Float_t wedgeinnerwidth;
  Float_t wedgeouterwidth;

  //whole detector coords:
  Float_t X1[4], Y1[4], X2[4], Y2[4]; //coords for two outer corners of each det
  Float_t Th1[4], Ph1[4], Th2[4], Ph2[4]; //corresponding theta and phi
  Float_t Z[4]; //z position, relative to most-downstream dets

  //individual channel coords:
  //0-15 are fronts, 16-31 are backs
  Float_t X1c[4][32], Y1c[4][32], X2c[4][32], Y2c[4][32]; //outer corners
  Float_t Th1c[4][32], Ph1c[4][32], Th2c[4][32], Ph2c[4][32];
  Float_t Rc_inner[32], Rc_outer[32]; //radii

};

/*

QQQ layout:

fronts (rings):
 R_OUTER -- 0
 THETA2  --
	 --
	 --
	 --
	 --
	 --
	 --    ^
	 --    |
	 --    | R & theta
	 --    |
	 --
	 --
	 --
	 --
 R_INNER -- 15
 THETA1
                

backs (wedges/strips):
      15||||||||||||||||0
     PHI2  <--- CCW   PHI1

Here, X1,Y1 ~ PHI1 and X2,Y2 ~ PHI2, unlike the SX3s which are backwards (see SX3.h)
*/

#endif
