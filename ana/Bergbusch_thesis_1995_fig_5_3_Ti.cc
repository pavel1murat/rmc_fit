///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetBergbusch_thesis_1995_fig_5_3_Ti(Data_t* Data) {
  // Bergbusch PhD thesis (1995)
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,134., 58.5,110., 59.5,114., 60.5,152., 61.5,109.,
    62.5, 95., 63.5,105., 64.5, 79., 65.5, 83., 66.5,101.,
    67.5, 70., 68.5, 70., 69.5, 50., 70.5, 52., 71.5, 62.,
    72.5, 42., 73.5, 29., 74.5, 37., 75.5, 29., 76.5, 20.,
    77.5, 22., 78.5, 17., 79.5, 19., 80.5,  9., 81.5,  4.,
    82.5,  7., 83.5,  6., 84.5,  2., 85.5,  3., 86.5,  4.,
    87.5,  2., 92.5,  1., 96.5,  1., 99.5,  1.,
    -1 };

  const char hist_name[] = "h_bergbusch_thesis_1995_Ti";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Bergbusch MS thesis (1995), figure 5.3, Ti");

  Data->fResp = triumf_response_1998;

  return 0;
}
