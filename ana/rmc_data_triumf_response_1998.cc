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
double rmc_data::r98_sigma0(double E) {
  double ps0[4] { 2.03 ,     10.4,     1.25,    -0.428 };
  
  double y      = (E-60.)/60.;
  double sig0   = ps0[0] + ps0[1]*y + ps0[2]*y*y + ps0[3]*y*y*y;
  return sig0;
}

//-----------------------------------------------------------------------------
double rmc_data::r98_sigma2(double E) {
  double ps2[4] { 0.786,    0.508,     0.425,   -0.164 };
  
  double y      = (E-60.)/60.;
  double sig2   = ps2[0] + ps2[1]*y + ps2[2]*y*y + ps2[3]*y*y*y;
  return sig2;
}

//-----------------------------------------------------------------------------
double rmc_data::r98_e0(double E) {
  double pe0[4] { 54.5 ,     57.7,    -0.315,      0.0 };
  
  double y      = (E-60.)/60.;
  double e0     = pe0[0] + pe0[1]*y + pe0[2]*y*y + pe0[3]*y*y*y;
  return e0;
}

//-----------------------------------------------------------------------------
double rmc_data::r98_e1(double E) {
  double e1 = rmc_data::r98_e0(E)-rmc_data::r98_sigma0(E);
  return e1;
}

//-----------------------------------------------------------------------------
double rmc_data::r98_e2(double E) {
  double s0 = rmc_data::r98_sigma0(E);
  double s2 = rmc_data::r98_sigma2(E);
  double e2 = rmc_data::r98_e0(E)+s0*s0/s2;
  return e2;
}

//-----------------------------------------------------------------------------
double rmc_data::r98_a(double E) {
  double pa [4] { 9.41e-4,  2.61e-3,  -0.270e-2, 0.835e-3};
  double y      = (E-60.)/60.;
  double a      = pa[0] + pa[1]*y + pa[2]*y*y + pa[3]*y*y*y;

  return a;
}


//-----------------------------------------------------------------------------
double rmc_data::r98_alp(double E) {
  double palp[4] { 56.1 ,   62.5,    -0.826,     0.0};

  double y      = (E-60.)/60.;
  double alp    = palp[0] + palp[1]*y + palp[2]*y*y + palp[3]*y*y*y;
  return alp;
}


//-----------------------------------------------------------------------------
double rmc_data::triumf_response_1998(double E, double Er) {

  double prob;

  //  double y  = (E-60.)/60.;

  //  printf(">>> triumf_response_1998: E, Er : %12.5e %12.5e %12.5e\n",E,Er,y);
//-----------------------------------------------------------------------------
// E0 parameterization
//-----------------------------------------------------------------------------
  double e0     = rmc_data::r98_e0    (E);
  double sig0   = rmc_data::r98_sigma0(E);
  double sig2   = rmc_data::r98_sigma2(E);

  double  e1    = rmc_data::r98_e1(E);
  double  e2    = rmc_data::r98_e2(E);

  //  printf(" sig0, sig2: %12.4e %12.4e\n", sig0, sig2);

  double a      = rmc_data::r98_a(E);

  //  printf("a : %12.5f\n", a);

  //  printf("e0, e1, e2: %12.5e %12.5e %12.5e\n",e0,e1,e2); 

  double alp     = rmc_data::r98_alp(E);

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


//-----------------------------------------------------------------------------
void rmc_data::plot_r98(const char* Parameter, double EMin, double EMax, int NBins) {
  TString par(Parameter);
  
  TH1F* hist = new TH1F(Form("r98_%s",Parameter),Form("r98: %s",Parameter),NBins,EMin,EMax);

  double bin = (EMax-EMin)/NBins;
    
  for (int i=1; i<=NBins; i++) {
    double e = EMin+(i-0.5)*bin;
    double f (0);

    if      (par == "sigma0"  ) f = rmc_data::r98_sigma0  (e);
    else if (par == "sigma2"  ) f = rmc_data::r98_sigma2  (e);
    else if (par == "a"       ) f = rmc_data::r98_a       (e);
    else if (par == "alp"     ) f = rmc_data::r98_alp     (e);
    else if (par == "e0"      ) f = rmc_data::r98_e0      (e);
    else if (par == "e1"      ) f = rmc_data::r98_e1      (e);
    else if (par == "e2"      ) f = rmc_data::r98_e2      (e);
    
    hist->SetBinContent(i,f);
  }

  hist->Draw();
}
