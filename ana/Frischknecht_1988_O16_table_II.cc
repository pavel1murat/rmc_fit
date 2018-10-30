///////////////////////////////////////////////////////////////////////////////
// this is a unfolded data plot, thus assume a delta-function response 
///////////////////////////////////////////////////////////////////////////////
#include  "ana/rmc_data.hh"
#include "TAxis.h"
int rmc_data::GetFrischknecht_1988_O16_table_II(Data_t* Data) {
  // PRC v38 i5 p1996
  // data from Table II
  // data in the format (x1,y1,ey1,x2,y2,ey2,....,xn,yn,eyn,-1)
  
  int   nbins(50);
  float xmin(1), xmax(101);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {

    58, 4.20, 0.42, 60, 2.32, 0.39, 62, 2.53, 0.31, 64, 1.89, 0.25, 66, 1.72, 0.24,
    68, 1.60, 0.21, 70, 1.48, 0.21, 72, 1.30, 0.17, 74, 0.86, 0.16, 76, 0.60, 0.15,
    78, 0.53, 0.15, 80, 0.45, 0.14, 82, 0.40, 0.13, 84, 0.34, 0.11, 86, 0.27, 0.10,
    88, 0.13, 0.12, 90, 0.19, 0.09, 92, 0.27, 0.07, 94, 0.07, 0.11, 96, 0.03, 0.11,
    98, 0.02, 0.09,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Frischknecht_1988_O16","",nbins,xmin,xmax);

  for (int np=0; data[3*np] > 0; np++) {
    int bin = (data[3*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[3*np+1]);
    Data->fHist->SetBinError  (bin,data[3*np+2]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PRC v38 i5 p1996 table II, ^{16}O");

  Data->fResp = rmc_data::delta_function_response;

  return 0;
}
