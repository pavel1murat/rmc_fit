///////////////////////////////////////////////////////////////////////////////
// T.P.Gorringe et al Phys Rev C v58 i3 p 1767 (1998 ) . Armstrong
// data for 1mm converter
///////////////////////////////////////////////////////////////////////////////
#include "murat/plot/smooth.hh"
#include "ana/rmc_fit.hh"
#include "TCanvas.h"
#include "Getline.h"

ClassImp(rmc_fit)
//-----------------------------------------------------------------------------
rmc_fit::rmc_fit(const char* Name): TNamed(Name,Name) {
  fRmcData = new rmc_data();
  fData.fHist = nullptr;
  fData.fResp = nullptr;
  fNCanvases  = 0;
  fNFunctions = 0;
}

//-----------------------------------------------------------------------------
rmc_fit::~rmc_fit() {
  delete fRmcData;
  if (fData.fHist) delete fData.fHist;
}

//-----------------------------------------------------------------------------
double rmc_fit::fit_pol2(double* X, double* P) {

  double dx = X[0]-P[2];

  double f = P[0]+P[1]*dx*dx;

  return f;
}

//-----------------------------------------------------------------------------
// closure approximation:
//
// P[0] : normalization constant
// P[1] : kMax
//-----------------------------------------------------------------------------
double rmc_fit::fun_closure(double* X, double* P) {
  double f{0};

  double kmax = P[1];

  double x = X[0]/kmax;
  
  if (x <= 1) {
    f = P[0]*(1-2*x+2*x*x)*x*(1-x)*(1-x);
  }

  return f;
}

//-----------------------------------------------------------------------------
// histogram is assumed to be created by the caller
// data are also assumed to be retrieved
//-----------------------------------------------------------------------------
int rmc_fit::get_response_hist(double Response(double,double), double E, TH1F* Hist) {

  Hist->Reset();

  int nbins = Hist->GetNbinsX();

  double e, p;
  for (int i=1; i<=nbins; i++) {
    e = Hist->GetBinCenter(i);
    p = Response(E,e);

    Hist->SetBinContent(i,p);
    //    printf(">>> get_response_hist: e,E,p : %12.4e %12.4e %12.4e \n",e,E,p);
  }
  return 0;
}

//-----------------------------------------------------------------------------
// convolute closure approximation spectrum with the given 'KMax' with the
// data response function, return resulting TF1
// 'Func' has just one parameter: overall normalization
//-----------------------------------------------------------------------------
int rmc_fit::get_convoluted_closure_spectrum(double KMax,double Response(double,double),
					  TF1** Func, int Debug) {
  
  fNFunctions++;

  TF1* f_closure = new TF1(Form("f_closure_%s_%i",GetName(),fNFunctions),fun_closure,0,100,2);
  f_closure->SetParameter(0,10);
  f_closure->SetParameter(1,KMax);
  f_closure->SetParError(1,0.);
//-----------------------------------------------------------------------------
// for given kMax and response, define the expected spectrum
//-----------------------------------------------------------------------------    
  while (TObject* o = gROOT->FindObject("h_resp")) delete o;
  TH1F* h_resp = new TH1F("h_resp","h_resp",3000,0,300);

  while (TObject* o = gROOT->FindObject("h_sp")) delete o;
  
  int nbins_spectrum{1100};
  
  TH1F* h_spectrum = new TH1F("h_sp","h_sp",nbins_spectrum,0,110);

  int    nj  = 10000;
  double qnj = nj;

  for (int i=1; i<=1000; i++) {
    double e = (i-1/2)*0.1;

    if (e < 40) continue;

    double w = f_closure->Eval(e);

    get_response_hist(Response,e,h_resp);
    double eff = h_resp->Integral();

    for (int j=0; j<nj; j++) {
      double r = h_resp->GetRandom();
      h_spectrum->Fill(r,eff*w/qnj);
    }
  }

  h_spectrum->Rebin(10);
  smooth* smf;

  smf = new smooth(h_spectrum,0,100);
  *Func = smf->GetFunc();

  if (Debug != 0) h_spectrum->Draw("e,p");

  return 0;
}

//-----------------------------------------------------------------------------
// testing procedures
// test3: convolute the closure approximation spectrum with given 'KMax' value
//        with a response function corresponding to the data defined by Year and Target
//        *doesn't make much sense* ...
//-----------------------------------------------------------------------------
int rmc_fit::plot_convoluted_closure_approximation_spectrum(const  char* ResponseModel,
							    double KMax               ) {
  TF1* f;
  double (*resp)(double,double) {NULL};
  
  int debug(1);

					// use Year / Target to get the response function
  
  fRmcData->get_response_function(ResponseModel,&resp);

					// convolute closure approximation with the response function
  
  get_convoluted_closure_spectrum(KMax,resp,&f,debug);

  f->SetLineWidth(1);
  f->Draw("same");

  return 0;
}

//-----------------------------------------------------------------------------
// print detector response for a given photon (e+e-) energy 'E' +/- 10 MeVx
//-----------------------------------------------------------------------------
void rmc_fit::test0(const char* ResponseModel, double E) {

  double (*resp)(double,double) {NULL};

  fRmcData->get_response_function(ResponseModel,&resp);

  double Er = E-10;

  double p  = resp(E,Er);

  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);

  Er        = E;
  p         = resp(E,Er);
  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);

  Er        = E+10;
  p         = resp(E,Er);
  printf(" E, Er, p = %12.5e %12.5e %12.5e\n",E,Er,p);
}

//-----------------------------------------------------------------------------
// test1: detector response for a given photon (e+e-) energy 'E'
// assume 'Al' data are always there
//-----------------------------------------------------------------------------
void rmc_fit::plot_response_function(const char* ResponseModel,
				     double E, double EMin, double EMax, int NBins) {

  double (*resp)(double,double) {NULL};

  fRmcData->get_response_function(ResponseModel,&resp);

  TH1F* h1 = new TH1F(Form("h_f_%s_%i",GetName(),fNFunctions),"h",NBins,EMin,EMax);
  fNFunctions++;

  get_response_hist(resp,E,h1);
  h1->Draw();
  
  printf(" integral: %12.5e\n",h1->Integral()*h1->GetBinWidth(0));
}

//-----------------------------------------------------------------------------
// test3: Armstrong'1999 data and '1998 Gorringe data,
// reasonable initial approximations below
//
// test3(1,86,2.5)
// test3(2,85,2.4)
// test3(3,89,3.2)
// test3(4,87,3.2)
// test3(5,86,3. )
//-----------------------------------------------------------------------------
void rmc_fit::fit(int         Year         ,
		  const char* Target       ,
		  const char* ResponseModel,
		  double      KMax         ,
		  double      MinFitE      ,
		  double      MaxFitE      ) {
  TF1* f;
  
  fRmcData->get_experimental_data(Year,Target,ResponseModel,&fData);
  
  fData.fHist->GetXaxis()->SetRangeUser(50,99.9);

  get_convoluted_closure_spectrum(KMax,fData.fResp,&f);
  
  double hist_int = fData.fHist->Integral()*fData.fHist->GetBinWidth(1);
  double fun_int  = f->Integral(55,100,1.e-5);
  double cnorm    = hist_int/fun_int/1.5;

  //    printf("hist_int, fun_int, cnorm: %12.5e %12.5e %12.5e\n",hist_int,fun_int,cnorm);
    
  f->SetParameter(0,cnorm);
    
  //    debugging
  //    fFun->Draw("same");
  //    c_fit->Modified();
  //    c_fit->Update();
  //    Getline("hit <RETURN>");

  fNCanvases++;
  TCanvas* c = new TCanvas(Form("c_%s_%i",GetName(),fNCanvases),"fit",1200,700);

  //  fData.fHist->Fit(f,"QL","p,e",MinFitE,MaxFitE);
  fData.fHist->Fit(f,"Q","p,e",MinFitE,MaxFitE);

  c->Modified();
  c->Update();
    
  double chi2    = f->GetChisquare();
  int    ndof    = f->GetNDF();
  double chi2dof = chi2/(ndof+1.e-12);
  double anorm   = f->GetParameter(0);
  
  printf("kMax, chi2, ndof, chi2/ndof, anorm = %12.5f %12.5e %3i %12.5e %12.5e\n",
	 KMax, chi2, ndof, chi2dof, anorm);

  for(int i=57; i<=100; i++) {
    double x = i+0.5;
    double dat = fData.fHist->GetBinContent(i+1);
    double fv = f->Eval(x);
    printf("x, dat, fv : %12.5f %12.5f %12.5f \n",x,dat,fv);
  }
}


//-----------------------------------------------------------------------------
void rmc_fit::nchi2(const TH1F* Hist, const TF1* Func, double* Chi2, int* NPt) {
  *Chi2 = 0;
  *NPt  = 0;
  
  int nb = Hist->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    float x = Hist->GetBinCenter(i);
    if (x < 57) continue;
    float y  = Hist->GetBinContent(i);
    float fy = Func->Eval(x);
					// skip empty bins, assume gaussian statistics
    if (y > 0) {
      double dy = fy-y;
      *Chi2 += dy*dy/y;
      *NPt  += 1;
    }
  }
}

//-----------------------------------------------------------------------------
// normalize the closure approximation to the number of events in the histogram
// and calculate the chi2...
//-----------------------------------------------------------------------------
void rmc_fit::nfit(int         Year         ,
		   const char* Target       ,
		   const char* ResponseModel,
		   double      KMax         ,
		   double      MinFitE      ,
		   double      MaxFitE      ) {
  TF1* f;
  
  fRmcData->get_experimental_data(Year,Target,ResponseModel,&fData);
  
  fData.fHist->GetXaxis()->SetRangeUser(50,99.9);

  get_convoluted_closure_spectrum(KMax,fData.fResp,&f);
  
  double hist_int = fData.fHist->Integral();
  double fun_int  = f->Integral(57,100,1.e-5);

  double cnorm    = hist_int/fun_int;
  printf("hist_int, fun_int, cnorm, p0: %12.5e %12.5e %12.5e %12.5e\n",
	 hist_int,fun_int,cnorm,f->GetParameter(0));
    
  f->SetParameter(0,cnorm/f->GetParameter(0));
    
  //    debugging
  //    fFun->Draw("same");
  //    c_fit->Modified();
  //    c_fit->Update();
  //    Getline("hit <RETURN>");

  fNCanvases++;
  TCanvas* c = new TCanvas(Form("c_%s_%i",GetName(),fNCanvases),"fit",1200,700);

  fData.fHist->Draw();
  f->Draw("sames");

  c->Modified();
  c->Update();

  double chi2;
  int ndof;

  nchi2(fData.fHist, f, &chi2, &ndof);
    
  printf("kMax, chi2, ndof, chi2/ndof: %12.5f %12.5e %3i %12.5e\n",
	 KMax, chi2, ndof, chi2/(ndof-1)); 
}


//-----------------------------------------------------------------------------
// for a fiven spectrum, perform a chi2 fit with fixed kMax and normalization
// varied as the only fit parameter
// scan range of kMax values [KMax1, KMax2] in 'NSteps' steps
// NSteps = 0: perform fit for only one value of kMax (KMax1)
//-----------------------------------------------------------------------------
void rmc_fit::scan(int Year, const char* Target, const char* ResponseModel,
		   double KMax1, double KMax2, int NSteps,
		   double MinFitE, double MaxFitE) {

  double* kmax    = new double [NSteps+1];
  double* chi2    = new double [NSteps+1];
  int*    ndof    = new int    [NSteps+1];
  double* chi2dof = new double [NSteps+1];
  double* anorm   = new double [NSteps+1];
//-----------------------------------------------------------------------------
// retrieve the experimental spectrum to fit
//-----------------------------------------------------------------------------
  fRmcData->get_experimental_data(Year,Target,ResponseModel,&fData);
  
  fData.fHist->GetXaxis()->SetRangeUser(50,99.9);

  double hist_int = fData.fHist->Integral()*fData.fHist->GetBinWidth(1);

  fNCanvases++;
  TCanvas* c = new TCanvas(Form("c_%s_%i",GetName(),fNCanvases),"fit",1800,700);
  c->Divide(2,1);
  c->cd(1);

  double step(0);
  if (NSteps > 0) step = (KMax2-KMax1)/NSteps;

  int imin        = -1;
  double chi2_min = 1.e6;

  TF1**   f = new TF1* [NSteps+1];
  
  for (int i=0; i<NSteps+1; i++) {
    kmax[i] = KMax1+i*step;
//-----------------------------------------------------------------------------
// get closure approximation spectrum convoluted with the detector response
//-----------------------------------------------------------------------------
    get_convoluted_closure_spectrum(kmax[i],fData.fResp,f+i);
//-----------------------------------------------------------------------------
// perform the fit
//-----------------------------------------------------------------------------
    fData.fHist->Draw("e,p");

    double fun_int   = f[i]->Integral(55,100,1.e-5);

    double cnorm = hist_int/fun_int/1.5;

    //    printf("hist_int, fun_int, cnorm: %12.5e %12.5e %12.5e\n",hist_int,fun_int,cnorm);
    
    f[i]->SetParameter(0,cnorm);
    
    //    debugging
    //    fFun->Draw("same");
    //    c_fit->Modified();
    //    c_fit->Update();
    //    Getline("hit <RETURN>");
    
    //    fData.fHist->Fit(f[i],"QL","",MinFitE,MaxFitE);
    fData.fHist->Fit(f[i],"Q","",MinFitE,MaxFitE);
    
    //    c_fit->Modified();
    //    c_fit->Update();

    chi2   [i] = f[i]->GetChisquare();
    ndof   [i] = f[i]->GetNDF();
    chi2dof[i] = chi2[i]/(ndof[i]+1.e-12);
    anorm  [i] = f[i]->GetParameter(0);
    
    printf("kMax, chi2, ndof, chi2/ndof, anorm = %12.5f %12.5e %3i %12.5e %12.5e\n",
    	   kmax[i], chi2[i], ndof[i], chi2dof[i], anorm[i]);

    if (chi2[i] < chi2_min) {
      chi2_min = chi2[i];
      imin     = i;
    }
    
    //    Getline("hit <RETURN>");
  }
//-----------------------------------------------------------------------------
// at this point have chi2 vs kmax, fit
//-----------------------------------------------------------------------------
  fGraphChi2 = new TGraph(NSteps+1,kmax,chi2);
  
  fGraphChi2->SetMarkerStyle(20);
  fGraphChi2->GetXaxis()->SetTitle("kmax, MeV");
  fGraphChi2->SetMarkerSize(1);
  fGraphChi2->SetTitle("chi2 vs kMax");
//-----------------------------------------------------------------------------
// don't need to fit the chi2 vs kMax in the whole range tested,
//-----------------------------------------------------------------------------
  TF1* f_pol2 = new TF1("fit_pol2",fit_pol2,KMax1,KMax2,3);
  f_pol2->SetParameter(0,chi2[imin]);
  f_pol2->SetParameter(1, 1);
  f_pol2->SetParameter(2,KMax1+(KMax2-KMax1)/NSteps*imin);

  c->cd(2);
  fGraphChi2->Fit(f_pol2,"","",KMax1,KMax2);
  fGraphChi2->Draw("AP");
//-----------------------------------------------------------------------------
// finally, plot the best fit
//-----------------------------------------------------------------------------
  c->cd(1);
  double kmax_best = f_pol2->GetParameter(2);

  TF1* fbest;
  get_convoluted_closure_spectrum(kmax_best,fData.fResp,&fbest);

  fbest->SetParameter(0,fData.fHist->Integral()/fbest->Integral(55,100,1.e-5)/1.5);
    
  //  fData.fHist->Fit(fbest,"L","",MinFitE,MaxFitE);
  fData.fHist->Fit(fbest,"","",MinFitE,MaxFitE);

  fData.fHist->Draw("ep");
//-----------------------------------------------------------------------------
// an attempt to clean up in the end
//-----------------------------------------------------------------------------
  delete [] kmax;
  delete [] chi2;
  delete [] anorm;
  delete [] ndof;
  delete [] chi2dof;
}

//-----------------------------------------------------------------------------
void rmc_fit::nscan(int Year, const char* Target, const char* ResponseModel,
		    double KMax1, double KMax2, int NSteps,
		    double MinFitE, double MaxFitE) {

  double* kmax    = new double [NSteps+1];
  double* chi2    = new double [NSteps+1];
  int*    ndof    = new int    [NSteps+1];
  double* chi2dof = new double [NSteps+1];
  double* anorm   = new double [NSteps+1];
//-----------------------------------------------------------------------------
// retrieve the experimental spectrum to fit
//-----------------------------------------------------------------------------
  fRmcData->get_experimental_data(Year,Target,ResponseModel,&fData);
  
  fData.fHist->GetXaxis()->SetRangeUser(50,99.9);

  double hist_int = fData.fHist->Integral()*fData.fHist->GetBinWidth(1);

  fNCanvases++;
  TCanvas* c = new TCanvas(Form("c_%s_%i",GetName(),fNCanvases),"fit",1800,700);
  c->Divide(2,1);
  c->cd(1);

  double step(0);
  if (NSteps > 0) step = (KMax2-KMax1)/NSteps;

  int imin        = -1;
  double chi2_min = 1.e6;

  TF1**   f = new TF1* [NSteps+1];
  
  for (int i=0; i<NSteps+1; i++) {
    kmax[i] = KMax1+i*step;
//-----------------------------------------------------------------------------
// get closure approximation spectrum convoluted with the detector response
//-----------------------------------------------------------------------------
    get_convoluted_closure_spectrum(kmax[i],fData.fResp,f+i);
//-----------------------------------------------------------------------------
// perform the fit
//-----------------------------------------------------------------------------
    fData.fHist->Draw("e,p");

    double fun_int   = f[i]->Integral(57,100,1.e-5);

    double cnorm = hist_int/fun_int/f[i]->GetParameter(0);

    //    printf("hist_int, fun_int, cnorm: %12.5e %12.5e %12.5e\n",hist_int,fun_int,cnorm);
    
    f[i]->SetParameter(0,cnorm);
    
    //    debugging
    //    fFun->Draw("same");
    //    c_fit->Modified();
    //    c_fit->Update();
    //    Getline("hit <RETURN>");
    
    nchi2(fData.fHist,f[i], &chi2[i], &ndof[i]);
    
    //    c_fit->Modified();
    //    c_fit->Update();

    ndof   [i] -= 1;
    chi2dof[i] = chi2[i]/(ndof[i]+1.e-12);
    anorm  [i] = f[i]->GetParameter(0);
    
    printf("kMax, chi2, ndof, chi2/ndof, anorm = %12.5f %12.5e %3i %12.5e %12.5e\n",
    	   kmax[i], chi2[i], ndof[i], chi2dof[i], anorm[i]);

    if (chi2[i] < chi2_min) {
      chi2_min = chi2[i];
      imin     = i;
    }
    
    //    Getline("hit <RETURN>");
  }
//-----------------------------------------------------------------------------
// at this point have chi2 vs kmax, fit
//-----------------------------------------------------------------------------
  fGraphChi2 = new TGraph(NSteps+1,kmax,chi2);
  
  fGraphChi2->SetMarkerStyle(20);
  fGraphChi2->GetXaxis()->SetTitle("kmax, MeV");
  fGraphChi2->SetMarkerSize(1);
  fGraphChi2->SetTitle("chi2 vs kMax");
//-----------------------------------------------------------------------------
// don't need to fit the chi2 vs kMax in the whole range tested,
//-----------------------------------------------------------------------------
  TF1* f_pol2 = new TF1("fit_pol2",fit_pol2,KMax1,KMax2,3);
  f_pol2->SetParameter(0,chi2[imin]);
  f_pol2->SetParameter(1, 1);
  f_pol2->SetParameter(2,KMax1+(KMax2-KMax1)/NSteps*imin);

  c->cd(2);
  fGraphChi2->Fit(f_pol2,"","",KMax1,KMax2);
  fGraphChi2->Draw("AP");
//-----------------------------------------------------------------------------
// finally, plot the best fit
//-----------------------------------------------------------------------------
  c->cd(1);
  double kmax_best = f_pol2->GetParameter(2);

  TF1* fbest;
  get_convoluted_closure_spectrum(kmax_best,fData.fResp,&fbest);

  double best_norm = fData.fHist->Integral()/fbest->Integral(55,100,1.e-5)/fbest->GetParameter(0);
  fbest->SetParameter(0,best_norm);
    
  //  fData.fHist->Fit(fbest,"L","",MinFitE,MaxFitE);
  fData.fHist->Fit(fbest,"","",MinFitE,MaxFitE);

  fData.fHist->Draw("ep");
//-----------------------------------------------------------------------------
// an attempt to clean up in the end
//-----------------------------------------------------------------------------
  delete [] kmax;
  delete [] chi2;
  delete [] anorm;
  delete [] ndof;
  delete [] chi2dof;
}

