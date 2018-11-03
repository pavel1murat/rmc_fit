///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetBergbusch_thesis_1995_fig_5_8_Zr(Data_t* Data) {
  // Bergbusch PhD thesis (1995)
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,82., 58.5,80., 59.5,56., 60.5,63., 61.5,59.,
    62.5,47., 63.5,59., 64.5,54., 65.5,50., 66.5,39.,
    67.5,20., 68.5,37., 69.5,31., 70.5,33., 71.5,28.,
    72.5,25., 73.5,22., 74.5,18., 75.5,10., 76.5,13.,
    77.5,12., 78.5, 8., 79.5, 5., 80.5,10., 81.5, 2.,
    82.5, 1., 83.5, 2., 84.5, 2., 86.5, 1., 87.5, 2.,
    88.5, 1., 89.5, 1., 92.5, 1.,
    -1 };

  const char hist_name[] = "h_bergbusch_thesis_1995_Zr";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Bergbusch MS thesis (1995), figure 5.8, Zr");

  Data->fResp = triumf_response_1998;

  return 0;
}
