#ifndef __rmc_fit__
#define __rmc_fit__

#include "TObject.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"

#include "rmc_data.hh"

class rmc_fit : public TNamed {
public:
					// data members
  TGraph*          fGraphChi2;
  TF1*             fFun;
  int              fYear;
  TString          fTarget;		// cached name, convenience
  rmc_data*        fRmcData;
  rmc_data::Data_t fData;
  int              fNCanvases;
  int              fNFunctions;
//  -----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  rmc_fit(const char* Name = "rmc_fit");
  ~rmc_fit();

  // static double response_1992(double E, double Er);
  // static double response_1998(double E, double Er);

  //  static double fun_response_1998    (double* X, double* P);

  static double fit_pol2             (double* X, double* P);
  static double fun_closure          (double* X, double* P);

  int           get_response_hist    (double Response(double,double), double E, TH1F* Hist);

  int           get_smeared_closure_spectrum(double KMax, double Response(double,double),
					     TF1** Func, int Debug = 0);

					// fit functions
  
  void          fit            (int Year, const char* Target, double KMax = 90.,
				double MinFitE = 57., double MaxFitE = 100.);

  void          scan_kmax_range(int Year, const char* Target, double EMin, double EMax,
				int NSteps=10, double MinFitE = 57, double MaxFitE=100);

					// test functions

  void          test0(int Year, const char* Target, double E);
  void          test1(int Year, const char* Target, double E, double EMin, double EMax, int NBins=100);
  void          test2(int Year, const char* Target, double E);
  int           test3(int Year, const char* Target, double KMax);

  ClassDef(rmc_fit,0)
};
#endif
