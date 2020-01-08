#include "dalitz.h"

MakeEff::MakeEff() {
  depth = 1;
}

MakeEff::~MakeEff() {

}

bool MakeEff::checkBinValid(Int_t gbin) {
  if(hfree->GetBinContent(gbin) != 0) {
    return true;
  } else {
    return false;
  }
}

Double_t MakeEff::interpBinContent(Int_t binx, Int_t biny, Int_t gbin) {
  if(checkBinValid(gbin)) {
    Double_t sum = 0;
    Double_t count = 0;
    for(int i = binx-depth; i<binx+depth+1; i++) {
      for(int j = biny-depth; j<biny+depth+1; j++) {
        if(i <= 0 || i> heff->GetNbinsX()) continue;
        if(j <= 0 || j> heff->GetNbinsY()) continue;
        if((i<binx+depth && i>binx-depth) && (j<biny+depth && j>biny-depth)) continue;
        Int_t this_bin = heff->GetBin(i,j);
        Double_t value = heff->GetBinContent(this_bin);
        auto iter = interp_set.find(this_bin);
        if((iter == interp_set.end()) || (value == 0 && checkBinValid(this_bin))) {
          continue;
        } else {
          sum += value;
          count++;
        }
      }
    }
    if(sum == 0) {
      depth++;
      return interpBinContent(binx, biny, gbin);
    } else {
      return sum/count;
    }
  } else {
    return 0;
  }
}

TH2F* MakeEff::makeEffHisto(TH2F *free, TH2F *detect, int Ebin) {
  hfree = free;
  hdetect = detect;

  heff = (TH2F*) hdetect->Clone(Form("DPeff_%ibin", Ebin));
  heff->Divide(hfree);
  for(int i=0; i<heff->GetNbinsX()+1; i++) {
    for(int j=0; j<heff->GetNbinsY()+1; j++) {
      Int_t this_bin = heff->GetBin(i,j);
      interp_set.insert(this_bin);
      Double_t value = heff->GetBinContent(this_bin);
      if(value == 0) {
        depth = 1;
        Double_t interp_value = interpBinContent(i,j,this_bin);
        heff->SetBinContent(this_bin, interp_value);
      }
    }
  }
  return heff;
}
