//-----------------------------------------------------------------------------
//  polynomial fits from Phys Rev C v58 i3 p 1767 (1998):
//  ----------------------------------------------------
//       a0        a1     a2            a3
// a    56.1      62.5    -0.826       0.0
// A    9.41e-4  2.61e-3  -0.270e-2    0.835e-3
// E0   54.5      57.7    -0.315       0.0
// s0   2.03      10.4     1.25       -0.428
// s2   0.786     0.508   0.425       -0.164
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include  "TMath.h"

//-----------------------------------------------------------------------------
double rmc_data::triumf_response_1998(double E, double Er) {

  double prob;

  double y  = (E-60.)/60.;

  //  printf(">>> triumf_response_1998: E, Er : %12.5e %12.5e %12.5e\n",E,Er,y);
//-----------------------------------------------------------------------------
// E0 parameterization
//-----------------------------------------------------------------------------
  double pe0[4] = { 54.5 ,     57.7,    -0.315,      0.0 };
  double e0     = pe0[0] + pe0[1]*y + pe0[2]*y*y + pe0[3]*y*y*y;

  double ps0[4] = { 2.03 ,     10.4,     1.25,    -0.428 };
  double sig0   = ps0[0] + ps0[1]*y + ps0[2]*y*y + ps0[3]*y*y*y;

  double ps2[4] = { 0.786,    0.508,     0.425,   -0.164 };
  double sig2   = ps2[0] + ps2[1]*y + ps2[2]*y*y + ps2[3]*y*y*y;

  double  e1    = e0-sig0;
  double  e2    = e0+sig0*sig0/sig2;

  //  printf(" sig0, sig2: %12.4e %12.4e\n", sig0, sig2);

  double pa [4] = { 9.41e-4,  2.61e-3,  -0.270e-2, 0.835e-3};
  double a      = pa[0] + pa[1]*y + pa[2]*y*y + pa[3]*y*y*y;

  //  printf("a : %12.5f\n", a);

  //  printf("e0, e1, e2: %12.5e %12.5e %12.5e\n",e0,e1,e2); 

  double palp[4] = { 56.1 ,   62.5,    -0.826,     0.0};
  double alp     = palp[0] + palp[1]*y + palp[2]*y*y + palp[3]*y*y*y;

  double x;
  
  if      (Er < e1) {
    double x0   = (alp-e1)/(alp-37.);
    double beta = a*TMath::Exp(-0.5)*x0/TMath::Log(x0);
    x           = (alp - Er)/(alp-37);
    prob        = beta*TMath::Log(x)/x;
  }
  else if (Er < e2) {
    x    = (Er-e0)/sig0;
    prob = a*TMath::Exp(-x*x/2.);
  }
  else {
    x    = (Er-e0)/sig2;
    prob = a*TMath::Exp(-x/2.);
    //    printf(" emoe: x, Er, e0, sig2: %12.5e %12.5e %12.5e %12.5e\n", x,Er,e0,sig2);
  }

  //  double slope = 0.015;
  double delta = 1; // -(E-70)*slope;
  
  prob = prob*delta;

  if (prob < 0) prob = 0;
  
  return prob;
}


//-----------------------------------------------------------------------------
double rmc_data::fun_triumf_response_1998(double* X, double* P) {
  double E  = X[0];
  double Er = X[1];

  double f = triumf_response_1998(E,Er);
  return f;
}
