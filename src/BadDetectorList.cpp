#include "BadDetectorList.h"

using namespace std;

BadDetectorList::BadDetectorList() {

  ifstream inList;
  inList.open("badDetectors.txt");

  if (!inList) {
    
    cerr << "*** WARNING: badDetectors.txt NOT FOUND! Initializing to zero... ***\n";
    total = -1; //error value
    DetList = new Int_t*[1];
    DetList[0] = new Int_t[2];
    DetList[0][0] = 0;
    DetList[0][1] = 0;

  }

  else {

    string junk;
    
    inList >> junk; //for "Total:" part of header
    inList >> total;

    inList >> junk >> junk; //for "DetID  ChNum" part of header

    DetList = new Int_t*[total];

    for (Int_t i=0; i<total; i++) {
      
      DetList[i] = new Int_t[2];
      inList >> DetList[i][0]; //DetID
      inList >> DetList[i][1]; //ChNum

    }

  }

}

BadDetectorList::~BadDetectorList() {

  for (Int_t i=0; i<total; i++)
    delete [] DetList[i];

  delete [] DetList;

}

Int_t BadDetectorList::GetTotal() {
  
  return total;

}

Bool_t BadDetectorList::IsBad(Int_t iDetID) {

  Int_t numbad = 0;

  for (Int_t i=0; i<total; i++)
    if (iDetID == DetList[i][0])
      numbad++; //count how many bad instances of this detector we have

  if (iDetID > 0 && iDetID < 3) {//QQQ
    if (numbad == 32)
      return true;
  }
  else //SX3
    if (numbad == 12)
      return true;
  
  //otherwise it's fine ("not bad")
  return false;

}

Bool_t BadDetectorList::IsBad(Int_t iDetID, Int_t iChNum) {

  for (Int_t i=0; i<total; i++) {

    //check if given detector/channel combo is "bad"
    if (iDetID == DetList[i][0] && iChNum == DetList[i][1])
      return true;

  }

  //otherwise it's fine ("not bad")
  return false;

}
