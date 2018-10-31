//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni62(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,62., 58.5,50., 59.5,68., 60.5,54., 61.5,62.,
    62.5,52., 63.5,44., 64.5,57., 65.5,58., 66.5,50.,
    67.5,43., 68.5,31., 69.5,32., 70.5,30., 71.5,20.,
    72.5,29., 73.5,15., 74.5,16., 75.5,12., 76.5,11.,
    77.5,15., 78.5, 6., 79.5, 5., 80.5, 6., 81.5, 7.,
    82.5, 2., 83.5, 1., 84.5, 4., 85.5, 2., 86.5, 1.,
    87.5, 0., 88.5, 1., 89.5, 1., 90.5, 0., 91.5, 0.,
    -1 };

  const char hist_name[] = "h_Gorringe_1998_Ni62";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width + 1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998) data from table II, ^{62}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni62_2MeV(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(50);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    59.0, 118.0, 61.0, 116.0, 63.0,  96.0, 65.0, 115.0, 67.0,  93.0,
    69.0,  63.0, 71.0,  50.0, 73.0,  44.0, 75.0,  28.0, 77.0,  26.0,
    79.0,  11.0, 81.0,  13.0, 83.0,   3.0, 85.0,   6.0, 87.0,   1.0,
    89.0,   2.0, 91.0,   0.0, 93.0,   0.0, 95.0,   0.0, 97.0,   0.0,
    99.0,   0.0,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Gorringe_1998_Ni62","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width + 1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998) data from table II, ^{62}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}
