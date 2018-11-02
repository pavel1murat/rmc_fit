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

  static double fit_pol2             (double* X, double* P);
  static double fun_closure          (double* X, double* P);

  int           get_response_hist    (double Response(double,double), double E, TH1F* Hist);

  int           get_convoluted_closure_spectrum(double KMax                   ,
						double Response(double,double),
						TF1**  Func                   ,
						int    Debug = 0              );

					// fit functions
  
  void          fit  (int Year, const char* Target, const char* ResponseModel,
		      double KMax = 90., double MinFitE = 57., double MaxFitE = 100.);

  void          nfit (int Year, const char* Target, const char* ResponseModel,
		      double KMax = 90., double MinFitE = 57., double MaxFitE = 100.);
  
  void          scan (int Year, const char* Target, const char* ResponseModel,
		      double EMin, double EMax,
		      int NSteps=10, double MinFitE = 57, double MaxFitE=100);
  
  void          nscan(int Year, const char* Target, const char* ResponseModel,
		      double EMin, double EMax,
		      int NSteps=10, double MinFitE = 57, double MaxFitE=100);
  
  void          nchi2(const TH1F* Hist, const TF1* Func, double* Chi2, int* NDof);
//-----------------------------------------------------------------------------
// test functions
//-----------------------------------------------------------------------------
  void          test0(const char* ResponseModel, double E);

  void          plot_response_function(const char* ResponseModel, double E,
				       double EMin=0, double EMax=300, int NBins=3000);
  
  int           plot_convoluted_closure_approximation_spectrum(const char* ResponseModel, double KMax);

  ClassDef(rmc_fit,0)
};
#endif
