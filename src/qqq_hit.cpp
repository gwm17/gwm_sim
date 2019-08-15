#include "qqq_hit.h"

bool qqq_hit(Float_t Xq, Float_t Yq) {

  Float_t rho = sqrt(Xq*Xq + Yq*Yq);

  //takes care of the gaps 
  if(Xq > -0.2505 && Xq < 0.2505)
    return false;
  if(Yq > -0.2505 && Yq < 0.2505)
    return false;
  //within rhos of qqqs
  if(rho < 5.01 || rho > 9.90)
    return false;

  return true;

}

Int_t qqq_phi_hit(Float_t Xd, Float_t Yd, Float_t phi, Float_t phi1, Float_t phi2) {

  if(Xd > 0 && Yd > 0 && phi>=phi1 && phi<=phi2)
    return 0;
  if(Xd <  0 && Yd > 0 && phi>=phi1 && phi<=phi2)
    return 3;
   if(Xd <  0 && Yd < 0 && phi>=phi1 && phi<=phi2)
    return 2;
   if(Xd >  0 && Yd < 0 && phi>=phi1 && phi<phi2)
    return 1;

   return -1;
}
