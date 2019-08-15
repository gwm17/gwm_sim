#include "SX3.h"

using namespace std;

SX3::SX3() : DETLENGTH(7.5), DETWIDTH(4.0) {

  frontlength = DETLENGTH;
  frontwidth = DETWIDTH/4;
  backlength = DETLENGTH/4;
  backwidth = DETWIDTH;

  for (int i=0; i<24; i++) {
    X1[i] = 0;
    X2[i] = 0;
    Y1[i] = 0;
    Y2[i] = 0;
    Z1[i] = 0;
    Z2[i] = 0;
    R1[i] = 0;
    R2[i] = 0;
    Th1[i] = 0;
    Ph1[i] = 0;
    Th2[i] = 0;
    Ph2[i] = 0;
    for (int j=0; j<8; j++) {
      X1c[i][j] = 0;
      X2c[i][j] = 0;
      Y1c[i][j] = 0;
      Y2c[i][j] = 0;
      Z1c[i][j] = 0;
      Z2c[i][j] = 0;
      R1c[i][j] = 0;
      R2c[i][j] = 0;
      Th1c[i][j] = 0;
      Ph1c[i][j] = 0;
      Th2c[i][j] = 0;
      Ph2c[i][j] = 0;
    }
  }

  TString headers = "",
          comment = "";

  ifstream wcfile;
  wcfile.open("WorldCoordinates2018Dec4");
  for (int i=0; i<7; i++) //skip first seven lines
    for (int j=0; j<7; j++)
      wcfile >> headers;
  for (int i=0; i<24; i++)
    wcfile >> comment >> Z1[i] >> X1[i] >> X2[i] >> Y1[i] >> Y2[i] >> comment;
  wcfile.close();

  for (int det=0; det<24; det++) {

    Z2[det] = Z1[det] + DETLENGTH;
    R1[det] = TMath::Sqrt(X1[det]*X1[det] + Y1[det]*Y1[det]);
    R2[det] = TMath::Sqrt(X2[det]*X2[det] + Y2[det]*Y2[det]);
    Ph1[det] = TMath::ATan2(Y2[det],X2[det]); //unfortunate, but current
    Ph2[det] = TMath::ATan2(Y1[det],X1[det]);
    Th1[det] = TMath::ATan2(R1[det],Z1[det]);
    Th2[det] = TMath::ATan2(R2[det],Z2[det]);

    if (Ph1[det] < 0) Ph1[det] += 2*TMath::Pi();
    if (Ph2[det] < 0) Ph2[det] += 2*TMath::Pi();

    Float_t phidet = Ph1[det] - Ph2[det]; //whole angle spanned by det

    for (int ch=0; ch<8; ch++) {

      if (ch < 4) { //backs
	Z1c[det][ch] = Z1[det] + ch*backlength;
	Z2c[det][ch] = Z1[det] + (ch+1)*backlength;
	X1c[det][ch] = X1[det];
	X2c[det][ch] = X2[det];
	Y1c[det][ch] = Y1[det];
	Y2c[det][ch] = Y2[det];
	R1c[det][ch] = R1[det];
	R2c[det][ch] = R2[det];
	Ph1c[det][ch] = Ph1[det];
	Ph2c[det][ch] = Ph2[det];;
	Th1c[det][ch] = TMath::ATan2(R1c[det][ch],Z1c[det][ch]);
	Th2c[det][ch] = TMath::ATan2(R2c[det][ch],Z2c[det][ch]);

	if (Ph1c[det][ch] < 0) Ph1c[det][ch] += 2*TMath::Pi();
	if (Ph2c[det][ch] < 0) Ph2c[det][ch] += 2*TMath::Pi();

      }
      else { //fronts
	Z1c[det][ch] = Z1[det];
	Z2c[det][ch] = Z2[det];
	Th1c[det][ch] = Th1[det];
	Th2c[det][ch] = Th2[det];
	Ph1c[det][ch] = Ph2[det] + (ch-4+1)*phidet/4;
	Ph2c[det][ch] = Ph2[det] + (ch-4)*phidet/4;

	if (Ph1c[det][ch] < 0) Ph1c[det][ch] += 2*TMath::Pi();
	if (Ph2c[det][ch] < 0) Ph2c[det][ch] += 2*TMath::Pi();

      }

    }

    //have to do more clever calcs for the rest of the front strip params
    R2c[det][4] = R2[det]; //outer radii
    R1c[det][7] = R1[det];
    R2c[det][6] = R1[det]*TMath::Cos(phidet/2); //middle radius (right triangle)
    R1c[det][5] = R2c[det][6];
    //now for other two radii (equal):
    R1c[det][4] = R2c[det][6]/TMath::Cos(phidet/4);
    R2c[det][5] = R1c[det][4];
    R1c[det][6] = R1c[det][4];
    R2c[det][7] = R1c[det][4];

    for (int ch=4; ch<8; ch++) {
      X1c[det][ch] = R1c[det][ch]*TMath::Cos(Ph2c[det][ch]);
      Y1c[det][ch] = R1c[det][ch]*TMath::Sin(Ph2c[det][ch]);
      X2c[det][ch] = R2c[det][ch]*TMath::Cos(Ph1c[det][ch]);
      Y2c[det][ch] = R2c[det][ch]*TMath::Sin(Ph1c[det][ch]);
    }

  }

}

Float_t SX3::GetX1(Int_t detid) {return (detid >= 6 && detid < 30 ? X1[detid-6] : 0.0);}
Float_t SX3::GetY1(Int_t detid) {return (detid >= 6 && detid < 30 ? Y1[detid-6] : 0.0);}
Float_t SX3::GetX2(Int_t detid) {return (detid >= 6 && detid < 30 ? X2[detid-6] : 0.0);}
Float_t SX3::GetY2(Int_t detid) {return (detid >= 6 && detid < 30 ? Y2[detid-6] : 0.0);}
Float_t SX3::GetZ1(Int_t detid) {return (detid >= 6 && detid < 30 ? Z1[detid-6] : 0.0);}
Float_t SX3::GetZ2(Int_t detid) {return (detid >= 6 && detid < 30 ? Z2[detid-6] : 0.0);}
Float_t SX3::GetR1(Int_t detid) {return (detid >= 6 && detid < 30 ? R1[detid-6] : 0.0);}
Float_t SX3::GetR2(Int_t detid) {return (detid >= 6 && detid < 30 ? R2[detid-6] : 0.0);}
Float_t SX3::GetTh1(Int_t detid) {return (detid >= 6 && detid < 30 ? Th1[detid-6] : 0.0);}
Float_t SX3::GetPh1(Int_t detid) {return (detid >= 6 && detid < 30 ? Ph1[detid-6] : 0.0);}
Float_t SX3::GetTh2(Int_t detid) {return (detid >= 6 && detid < 30 ? Th1[detid-6] : 0.0);}
Float_t SX3::GetPh2(Int_t detid) {return (detid >= 6 && detid < 30 ? Ph1[detid-6] : 0.0);}

Float_t SX3::GetX1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? X1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetY1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Y1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetX2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? X2c[detid-6][chnum] : 0.0);}
Float_t SX3::GetY2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Y2c[detid-6][chnum] : 0.0);}
Float_t SX3::GetZ1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Z1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetZ2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Z2c[detid-6][chnum] : 0.0);}
Float_t SX3::GetR1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? R1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetR2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? R2c[detid-6][chnum] : 0.0);}
Float_t SX3::GetTh1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Th1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetPh1(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Ph1c[detid-6][chnum] : 0.0);}
Float_t SX3::GetTh2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Th2c[detid-6][chnum] : 0.0);}
Float_t SX3::GetPh2(Int_t detid, Int_t chnum)
{return (detid >= 6 && detid < 30 && chnum >= 0 && chnum < 8 ? Ph2c[detid-6][chnum] : 0.0);}
