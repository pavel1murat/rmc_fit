///////////////////////////////////////////////////////////////////////////////
// PHYSICAL REVIEW C VOLUME 32, NUMBER 5, page 1506 NOVEMBER 1985
// this is a unfolded data plot, thus assume a delta-function response 
///////////////////////////////////////////////////////////////////////////////
#include  "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetFrischknecht_1985_Ca40_table_IV(Data_t* Data) {
  // PRC v38 i5 p1996
  // data from Table II
  // data in the format (x1,y1,ey1,x2,y2,ey2,....,xn,yn,eyn,-1)
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    55.5, 3.06, 0.31, 56.5, 2.46, 0.26, 57.5, 1.93, 0.22, 58.5, 1.59, 0.20, 59.5, 1.46, 0.19,
    60.5, 1.41, 0.17, 61.5, 1.37, 0.17, 62.5, 1.05, 0.16, 63.5, 0.89, 0.14, 64.5, 0.84, 0.13,    
    65.5, 0.87, 0.12, 66.5, 0.97, 0.12, 67.5, 0.93, 0.12, 68.5, 0.85, 0.11, 69.5, 0.60, 0.10,     
    70.5, 0.57, 0.09, 71.5, 0.68, 0.09, 72.5, 0.56, 0.09, 73.5, 0.46, 0.08, 74.5, 0.55, 0.09,
    75.5, 0.37, 0.08, 76.5, 0.49, 0.08, 77.5, 0.34, 0.07, 78.5, 0.26, 0.06, 79.5, 0.31, 0.07,
    80.5, 0.15, 0.06, 81.5, 0.32, 0.06, 82.5, 0.12, 0.05, 83.5, 0.17, 0.05, 84.5, 0.17, 0.05,
    85.5, 0.21, 0.05, 86.5, 0.10, 0.04, 87.5, 0.05, 0.04, 88.5, 0.06, 0.04, 89.5, 0.03, 0.04,
    90.5, 0.00, 0.04, 91.5, 0.06, 0.04, 92.5,-0.01, 0.04, 93.5,-0.08, 0.04, 94.5, 0.01, 0.04,
    95.5, 0.03, 0.04, 96.5, 0.04, 0.04, 97.5,-0.03, 0.04, 98.5,-0.02, 0.03, 99.5, 0.01, 0.04,
    -1 };

  const char hist_name[] = "h_Frischknecht_1988_O16";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[3*np] > 0; np++) {
    int bin = (data[3*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[3*np+1]);
    Data->fHist->SetBinError  (bin,data[3*np+2]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PRC v32 i5 p1506 table IV, ^{40}Ca, unfolded");

  Data->fResp = rmc_data::delta_function_response;

  return 0;
}
