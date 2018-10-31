//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni60(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,101., 58.5, 98., 59.5,102., 60.5,101., 61.5, 94.,
    62.5,107., 63.5, 91., 64.5,106., 65.5, 62., 66.5, 82.,
    67.5, 61., 68.5, 65., 69.5, 52., 70.5, 54., 71.5, 46.,
    72.5, 58., 73.5, 38., 74.5, 27., 75.5, 27., 76.5, 24.,
    77.5, 23., 78.5, 16., 79.5,  7., 80.5,  7., 81.5,  6.,
    82.5, 10., 83.5,  8., 84.5,  5., 85.5,  4., 86.5,  2.,
    87.5,  1., 88.5,  0., 89.5,  2., 90.5,  1., 91.5,  0.,
    92.5,  0., 93.5,  2., 94.5,  0., 95.5,  0., 96.5,  0.,
    97.5,  1., 98.5,  0., 99.5,  0.,
    -1 };

  const char hist_name[] = "h_Gorringe_1998_Ni60";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width + 1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998)  data from table II, ^{60}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni60_2MeV(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(50);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    59.0, 190.7, 61.0, 185.9, 63.0, 188.8, 65.0, 160.2, 67.0, 136.3,
    69.0, 111.5, 71.0,  95.3, 73.0,  91.5, 75.0,  51.5, 77.0,  44.8,
    79.0,  21.9, 81.0,  12.4, 83.0,  17.2, 85.0,   8.6, 87.0,   2.9,
    89.0,   1.9, 91.0,   1.0, 93.0,   1.9, 95.0,   0.0, 97.0,   1.0,
    99.0,   0.0,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Gorringe_1998_Ni60","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width + 1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998)  data from table II, ^{60}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}
