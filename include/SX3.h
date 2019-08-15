#ifndef __SX3_h
#define __SX3_h

#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TString.h>
#include <TMath.h>

class SX3 {

 public: 

  SX3();

  Float_t GetX1(Int_t detid);
  Float_t GetY1(Int_t detid);
  Float_t GetX2(Int_t detid);
  Float_t GetY2(Int_t detid);
  Float_t GetZ1(Int_t detid);
  Float_t GetZ2(Int_t detid);
  Float_t GetR1(Int_t detid);
  Float_t GetR2(Int_t detid);
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
  Float_t GetZ1(Int_t detid, Int_t chnum);
  Float_t GetZ2(Int_t detid, Int_t chnum);
  Float_t GetR1(Int_t detid, Int_t chnum);
  Float_t GetR2(Int_t detid, Int_t chnum);
  Float_t GetTh1(Int_t detid, Int_t chnum);
  Float_t GetPh1(Int_t detid, Int_t chnum);
  Float_t GetTh2(Int_t detid, Int_t chnum);
  Float_t GetPh2(Int_t detid, Int_t chnum);

 private:

  const Float_t DETLENGTH, DETWIDTH; //cm

  Float_t frontlength,
          frontwidth,
          backlength,
          backwidth;

  //whole detector coords (24 total, ids 6-29):
  Float_t X1[24], Y1[24], X2[24], Y2[24], Z1[24], Z2[24]; //cartesian
  Float_t R1[24], R2[24], Th1[24], Ph1[24], Th2[24], Ph2[24]; //polar
  //***this is CYLINDRICAL R, perp. to beam axis; i.e. R^2 = X^2 + Y^2 (NO Z^2)

  //individual channel coords (8 per detector):
  Float_t X1c[24][8], Y1c[24][8], X2c[24][8], Y2c[24][8], Z1c[24][8], Z2c[24][8];
  Float_t R1c[24][8], R2c[24][8], Th1c[24][8], Ph1c[24][8], Th2c[24][8], Ph2c[24][8];

};

/*

Here's a diagram of the SX3s:

front:                        back:
                                X2  3     2     1     0
  =========================      *======================== PHI1 & R2   |
  |                       |4     |     |     |     |     |             | PHI (CCW)
  |-----------------------|      |     |     |     |     |             v
  |                       |5     |     |     |     |     |
  |-----------------------|      |     |     |     |     |        <------BEAM-------
  |                       |6     |     |     |     |     |
  |-----------------------|      |     |     |     |     |
  |                       |7     |     |     |     |     |
  =========================      *======================== PHI2 & R1
  Z1                      Z2    X1
THETA1                  THETA2

Unfortunately, X1/Y1 ~ PHI2, and X2/Y2 ~ PHI1 in the current implementation

*/

#endif
