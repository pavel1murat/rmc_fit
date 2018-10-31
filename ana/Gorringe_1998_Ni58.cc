//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetGorringe_1998_Ni58(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5, 70., 58.5, 78., 59.5, 59., 60.5, 56., 61.5, 63.,
    62.5, 67., 63.5, 68., 64.5, 60., 65.5, 51., 66.5, 49.,
    67.5, 44., 68.5, 32., 69.5, 49., 70.5, 37., 71.5, 39.,
    72.5, 29., 73.5, 25., 74.5, 24., 75.5, 19., 76.5, 24.,
    77.5, 15., 78.5, 13., 79.5, 13., 80.5, 16., 81.5,  8.,
    82.5,  5., 83.5, 11., 84.5,  1., 85.5,  4., 86.5,  2.,
    87.5,  2., 88.5,  1., 89.5,  2., 90.5,  2., 91.5,  0,
    -1 };

  const char hist_name[] = "h_Gorringe_1998_Ni58";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998), Table II, ^{58}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni58_2MeV(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(50);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    59.0, 137.0, 61.0, 119.0, 63.0, 135.0, 65.0, 111.0, 67.0,  93.0,
    69.0,  81.0, 71.0,  76.0, 73.0,  54.0, 75.0,  43.0, 77.0,  39.0,
    79.0,  26.0, 81.0,  24.0, 83.0,  16.0, 85.0,   5.0, 87.0,   4.0,
    89.0,   3.0, 91.0,   2.0, 93.0,   0.0, 95.0,   0.0, 97.0,   0.0,    
    99.0,   0.0,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Gorringe_1998_Ni58","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998), Table II, ^{58}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}
