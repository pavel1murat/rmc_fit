//-----------------------------------------------------------------------------
// response of the TRIUMF detector as described in '1992 PhysRevC
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TMath.h"
//-----------------------------------------------------------------------------
double rmc_data::triumf_response_1992(double E, double Er) {

  double prob;

  //  printf("E, Er : %12.5e %12.5e\n",E,Er);

  double ps0[3] = { -0.5836,  0.0352 ,  0.       };
  double ps1[3] = { -5.879 ,  0.1653 , -5.149e-4 };
  double ps2[3] = {  1.596 , -0.03859,  3.883e-4 };
  double ps3[3] = {-47.80  ,  1.010  , -4.406e-3 };

  double sig0 = ps0[0] + ps0[1]*E + ps0[2]*E*E;
  double sig1 = ps1[0] + ps1[1]*E + ps1[2]*E*E;
  double sig2 = ps2[0] + ps2[1]*E + ps2[2]*E*E;
  double sig3 = ps3[0] + ps3[1]*E + ps3[2]*E*E;

  // printf(" sig0, sig1, sig2, sig3: %12.4e %12.4e %12.4e %12.4e\n",
  // 	 sig0, sig1, sig2, sig3);
  

  double pa [4] = {  3.259e-4, -4.120e-4, 1.015e-5, -4.050e-8};
  double pfa[3] = { -0.1337  , 2.828e-3, -9.701e-6};

  double a  = pa[0] + pa[1]*E + pa[2]*E*E + pa[3]*E*E*E;
  double fa = pfa[0] + pfa[1]*E + pfa[2]*E*E;

  // printf("a,fa : %12.5f %12.5f\n", a, fa);

  double pe0[2][3] = { 22.73 , 0.1995,  5.993e-3,
		       -1.161, 0.9481,  1.724e-4 };
  
  double e0;
  if (E < 60) e0   = pe0[0][0] + pe0[0][1]*E + pe0[0][2]*E*E;
  else        e0   = pe0[1][0] + pe0[1][1]*E + pe0[1][2]*E*E;

  double  e1 = e0-sig0*sig0/sig1;
  double  e2 = e0+sig0*sig0/sig2;

  double  e3 = 1.068 + 0.7507*E;

  // printf("e0, e1, e2, e3: %12.5e %12.5e %12.5e %12.5e\n",e0,e1,e2,e3); 

  double x, ba;
  
  double x3  = (Er-e3)/sig3;

  // printf ("x3 : %12.5f\n",x3);

  if (Er < e1) {

    double s01 = sig0/sig1;
    ba   = TMath::Exp(-s01*s01/2);
    x    = (e1-Er)/sig1;
    prob = a*(ba*TMath::Exp(-x/2.) + fa*TMath::Exp(-x3*x3/2));
  }
  else if (Er < e2) {
    x    = (Er-e0)/sig0;
    prob = a*(TMath::Exp(-x*x/2.) + fa*TMath::Exp(-x3*x3/2));
  }
  else {
    double s02 = sig0/sig2;
    double ca  = TMath::Exp(-s02*s02/2);

    x    = (Er-e2)/sig2;
    prob = a*(ca*TMath::Exp(-x/2.) + fa*TMath::Exp(-x3*x3/2));
  }

  return prob;
}
