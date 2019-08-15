#include "QQQ.h"

using namespace std;

QQQ::QQQ() : R_INNER(5.01), R_OUTER(9.90) {

  ringwidth = (R_OUTER - R_INNER)/16;
  wedgeinnerwidth = R_INNER*(TMath::Pi()/2)/16;
  wedgeouterwidth = R_OUTER*(TMath::Pi()/2)/16;

  for (int i=0; i<4; i++) {
    X1[i] = 0;
    X2[i] = 0;
    Y1[i] = 0;
    Y2[i] = 0;
    Z[i] = 0;
    Th1[i] = 0;
    Ph1[i] = 0;
    Th2[i] = 0;
    Ph2[i] = 0;
    for (int j=0; j<32; j++) {
      X1c[i][j] = 0;
      X2c[i][j] = 0;
      Y1c[i][j] = 0;
      Y2c[i][j] = 0;
      Th1c[i][j] = 0;
      Ph1c[i][j] = 0;
      Th2c[i][j] = 0;
      Ph2c[i][j] = 0;
      Rc_inner[j] = 0;
      Rc_outer[j] = 0;
    }
  }

  TString headers = "",
          comment = "";

  ifstream wcfile;
  wcfile.open("WorldCoordinates2018Dec4");
  for (int i=0; i<7; i++) wcfile >> headers;
  for (int i=0; i<4; i++)
    wcfile >> comment >> Z[i] >> X1[i] >> X2[i] >> Y1[i] >> Y2[i] >> comment;
  wcfile.close();

  for (int det=0; det<4; det++) {

    Ph1[det] = TMath::ATan2(Y1[det],X1[det]);
    Ph2[det] = TMath::ATan2(Y2[det],X2[det]);
    Th1[det] = TMath::ATan2(R_INNER,Z[det]);
    Th2[det] = TMath::ATan2(R_OUTER,Z[det]);

    if (Ph1[det] < 0) Ph1[det] += 2*TMath::Pi();
    if (Ph2[det] < 0) Ph2[det] += 2*TMath::Pi();

    for (int ch=0; ch<32; ch++) {

      if (ch < 16) {//front rings
	Rc_outer[ch] = R_OUTER - ch*ringwidth;
	Rc_inner[ch] = R_OUTER - (ch+1)*ringwidth;
	if (TMath::Abs(X1[det]) > TMath::Abs(X2[det])) {
	  X1c[det][ch] = X1[det] - ch*ringwidth;
	  Y1c[det][ch] = Y1[det];
	  X2c[det][ch] = X2[det];
	  Y2c[det][ch] = Y2[det] - ch*ringwidth;
	}
	else {
	  X1c[det][ch] = X1[det];
	  Y1c[det][ch] = Y1[det] - ch*ringwidth;
	  X2c[det][ch] = X2[det] - ch*ringwidth;
	  Y2c[det][ch] = Y2[det];
	}
	Ph1c[det][ch] = TMath::ATan2(Y1c[det][ch],X1c[det][ch]);
	Ph2c[det][ch] = TMath::ATan2(Y2c[det][ch],X2c[det][ch]);
	Th1c[det][ch] = TMath::ATan2(Rc_outer[ch],Z[det]);
	Th2c[det][ch] = TMath::ATan2(Rc_outer[ch],Z[det]);

	if (Ph1c[det][ch] < 0) Ph1c[det][ch] += 2*TMath::Pi();
	if (Ph2c[det][ch] < 0) Ph2c[det][ch] += 2*TMath::Pi();

      }

      else {//back wedges
	Int_t chshift = ch-16; //shift back to zero for the moment
	Rc_inner[ch] = R_INNER;
	Rc_outer[ch] = R_OUTER;
	Th1c[det][ch] = Th1[det];
	Th2c[det][ch] = Th2[det];
	Ph1c[det][ch] = Ph1[det] + chshift*(TMath::Pi()/2)/16;
	Ph2c[det][ch] = Ph1[det] + (chshift+1)*(TMath::Pi()/2)/16;
	X1c[det][ch] = R_OUTER*TMath::Cos(Ph1c[det][ch]);
	Y1c[det][ch] = R_OUTER*TMath::Sin(Ph1c[det][ch]);
	X2c[det][ch] = R_OUTER*TMath::Cos(Ph2c[det][ch]);
	Y2c[det][ch] = R_OUTER*TMath::Sin(Ph2c[det][ch]);

	if (Ph1c[det][ch] < 0) Ph1c[det][ch] += 2*TMath::Pi();
	if (Ph2c[det][ch] < 0) Ph2c[det][ch] += 2*TMath::Pi();
	
      }

    }

  }

}

Float_t QQQ::GetX1(Int_t detid) {return (detid >= 0 && detid < 4 ? X1[detid] : 0.0);}
Float_t QQQ::GetY1(Int_t detid) {return (detid >= 0 && detid < 4 ? Y1[detid] : 0.0);}
Float_t QQQ::GetX2(Int_t detid) {return (detid >= 0 && detid < 4 ? X2[detid] : 0.0);}
Float_t QQQ::GetY2(Int_t detid) {return (detid >= 0 && detid < 4 ? Y2[detid] : 0.0);}
Float_t QQQ::GetZ(Int_t detid) {return (detid >= 0 && detid < 4 ? Z[detid] : 0.0);}
Float_t QQQ::GetTh1(Int_t detid) {return (detid >= 0 && detid < 4 ? Th1[detid] : 0.0);}
Float_t QQQ::GetPh1(Int_t detid) {return (detid >= 0 && detid < 4 ? Ph1[detid] : 0.0);}
Float_t QQQ::GetTh2(Int_t detid) {return (detid >= 0 && detid < 4 ? Th2[detid] : 0.0);}
Float_t QQQ::GetPh2(Int_t detid) {return (detid >= 0 && detid < 4 ? Ph2[detid] : 0.0);}
Float_t QQQ::GetRInner(Int_t detid) {return (detid >= 0 && detid < 4 ? R_INNER : 0.0);}
Float_t QQQ::GetROuter(Int_t detid) {return (detid >= 0 && detid < 4 ? R_OUTER : 0.0);}

Float_t QQQ::GetX1(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? X1c[detid][chnum] : 0.0);}
Float_t QQQ::GetY1(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Y1c[detid][chnum] : 0.0);}
Float_t QQQ::GetX2(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? X2c[detid][chnum] : 0.0);}
Float_t QQQ::GetY2(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Y2c[detid][chnum] : 0.0);}
Float_t QQQ::GetZ(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Z[detid] : 0.0);}
Float_t QQQ::GetTh1(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Th1c[detid][chnum] : 0.0);}
Float_t QQQ::GetPh1(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Ph1c[detid][chnum] : 0.0);}
Float_t QQQ::GetTh2(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Th2c[detid][chnum] : 0.0);}
Float_t QQQ::GetPh2(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Ph2c[detid][chnum] : 0.0);}
Float_t QQQ::GetRInner(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Rc_inner[chnum] : 0.0);}
Float_t QQQ::GetROuter(Int_t detid, Int_t chnum)
{return (detid >= 0 && detid < 4 && chnum >= 0 && chnum < 32 ? Rc_outer[chnum] : 0.0);}
