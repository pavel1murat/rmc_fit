//-----------------------------------------------------------------------------
// response of the TRIUMF detector as described in '1992 PhysRevC
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TMath.h"


//-----------------------------------------------------------------------------
double rmc_data::r92_sigma0(double E) {
  double c[3]{ -0.5836,  0.0352 ,  0. };

  double s = c[0] + c[1]*E + c[2]*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_sigma1(double E) {
  double c[3]{ -5.879 ,  0.1653 , -5.149e-4 };

  double s = c[0] + c[1]*E + c[2]*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_sigma2(double E) {
  double c[3]{  1.596 , -0.03859,  3.883e-4 };
  double s = c[0] + c[1]*E + c[2]*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_sigma3(double E) {
  double c[3]{-47.80  ,  1.010  , -4.406e-3 };
  double s = c[0] + c[1]*E + c[2]*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_a(double E) {
  double c[4]{  3.259e-4, -4.120e-4, 1.015e-5, -4.050e-8};
  double s = c[0] + c[1]*E + c[2]*E*E+c[3]*E*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_f_over_a(double E) {
  double c[4]{ -0.1337  , 2.828e-3, -9.701e-6};
  double s = c[0] + c[1]*E + c[2]*E*E;
  return s;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_e0(double E) {
  
  double pe0[2][3] = { 22.73 , 0.1995,  5.993e-3,
		       -1.161, 0.9481,  1.724e-4 };

  double e0;
  if (E < 60) e0   = pe0[0][0] + pe0[0][1]*E + pe0[0][2]*E*E;
  else        e0   = pe0[1][0] + pe0[1][1]*E + pe0[1][2]*E*E;

  return e0;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_e1(double E) {

  double e0 = r92_e0    (E);
  double s0 = r92_sigma0(E);
  double s1 = r92_sigma1(E);
  
  return e0-s0*s0/s1;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_e2(double E) {

  double e0 = r92_e0    (E);
  double s0 = r92_sigma0(E);
  double s1 = r92_sigma2(E);
  
  return e0+s0*s0/s1;
}

//-----------------------------------------------------------------------------
double rmc_data::r92_e3(double E) {
  return 1.068 + 0.7507*E;
}

//-----------------------------------------------------------------------------
double rmc_data::triumf_response_1992(double E, double Er) {

  double prob;

  //  printf("E, Er : %12.5e %12.5e\n",E,Er);

  double sig0 = rmc_data::r92_sigma0(E);
  double sig1 = rmc_data::r92_sigma1(E);
  double sig2 = rmc_data::r92_sigma2(E);
  double sig3 = rmc_data::r92_sigma3(E);

  // printf(" sig0, sig1, sig2, sig3: %12.4e %12.4e %12.4e %12.4e\n",
  // 	 sig0, sig1, sig2, sig3);
  
  double a    = rmc_data::r92_a(E);
  double fa   = rmc_data::r92_f_over_a(E);

  double e0   = rmc_data::r92_e0(E);

  double e1   = rmc_data::r92_e1(E);
  double e2   = rmc_data::r92_e2(E);

  double e3   = rmc_data::r92_e3(E);

  // printf("a,fa : %12.5f %12.5f\n", a, fa);
  // printf("e0, e1, e2, e3: %12.5e %12.5e %12.5e %12.5e\n",e0,e1,e2,e3); 

  double x, ba;
  
  double x3  = (Er-e3)/sig3;

  // printf ("x3 : %12.5f\n",x3);

  if (Er < e1) {

    double s01 = sig0/sig1;
    ba   = TMath::Exp(-s01*s01/2);
    x    = (e1-Er)/sig1;
    //    prob = a*(ba*TMath::Exp(-x/2.) + fa*TMath::Exp(-x3*x3/2));
    prob = a*(ba*TMath::Exp(-x) + fa*TMath::Exp(-x3*x3/2));
  }
  else if (Er < e2) {
    x    = (Er-e0)/sig0;
    prob = a*(TMath::Exp(-x*x/2.) + fa*TMath::Exp(-x3*x3/2));
  }
  else {
    double s02 = sig0/sig2;
    double ca  = TMath::Exp(-s02*s02/2);

    x    = (Er-e2)/sig2;
    prob = a*(ca*TMath::Exp(-x) + fa*TMath::Exp(-x3*x3/2));
  }

  if (prob < 0) prob = 0;
  
  return prob;
}


//-----------------------------------------------------------------------------
void rmc_data::plot_r92(const char* Parameter, double EMin, double EMax, int NBins) {
  TString par(Parameter);
  
  TH1F* hist = new TH1F(Form("r92_%s",Parameter),Form("r92: %s",Parameter),NBins,EMin,EMax);

  double bin = (EMax-EMin)/NBins;
    
  for (int i=1; i<=NBins; i++) {
    double e = EMin+(i-0.5)*bin;
    double f (0);

    if      (par == "sigma0"  ) f = rmc_data::r92_sigma0  (e);
    else if (par == "sigma1"  ) f = rmc_data::r92_sigma1  (e);
    else if (par == "sigma2"  ) f = rmc_data::r92_sigma2  (e);
    else if (par == "sigma3"  ) f = rmc_data::r92_sigma3  (e);
    else if (par == "a"       ) f = rmc_data::r92_a       (e);
    else if (par == "f_over_a") f = rmc_data::r92_f_over_a(e);
    else if (par == "e0"      ) f = rmc_data::r92_e0      (e);
    else if (par == "e1"      ) f = rmc_data::r92_e1      (e);
    else if (par == "e2"      ) f = rmc_data::r92_e2      (e);
    else if (par == "e3"      ) f = rmc_data::r92_e3      (e);
    
    hist->SetBinContent(i,f);
  }

  hist->Draw();
}
