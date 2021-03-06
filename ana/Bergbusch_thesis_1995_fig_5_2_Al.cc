///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetBergbusch_thesis_1995_fig_5_2_Al(Data_t* Data) {
  // Bergbusch PhD thesis (1993)
  // digitized fig 6 : RMC spectrum on Al27
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = { // v2
    57.5,221, 58.5,223, 59.5,215, 60.5,172, 61.5,182, 
    62.5,160, 63.5,163, 64.5,164, 65.5,160, 66.5,142,
    67.5,140, 68.5,120, 69.5, 96, 70.5,101, 71.5, 78,
    72.5, 71, 73.5, 68, 74.5, 52, 75.5, 56, 76.5, 32,
    77.5, 32, 78.5, 27, 79.5, 19, 80.5, 18, 81.5, 12,
    82.5, 16, 83.5,  7, 84.5,  2, 85.5,  5, 86.5,  6, 
    87.5,  5, 88.5,  3, 89.5,  2, 91.5,  1, 96.5,  1,
    -1 };

  const char hist_name[] = "h_bergbusch_thesis_1995_Al";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Bergbusch MS thesis (1995), figure 5.2, ^{27}Al");

  Data->fResp = triumf_response_1998;

  return 0;
}
