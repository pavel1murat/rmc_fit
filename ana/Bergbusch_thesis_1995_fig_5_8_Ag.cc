///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetBergbusch_thesis_1995_fig_5_8_Ag(Data_t* Data) {
  // Bergbusch PhD thesis (1995)
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,62., 58.5,57., 59.5,73., 60.5,43., 61.5,60.,
    62.5,39., 63.5,44., 64.5,43., 65.5,42., 66.5,42.,
    67.5,27., 68.5,28., 69.5,34., 70.5,22., 71.5,16.,
    72.5,19., 73.5,16., 74.5,13., 75.5,14., 76.5, 8.,
    77.5,10., 78.5, 7., 79.5, 5., 80.5, 2., 81.5, 3.,
    82.5, 6., 83.5, 1., 84.5, 2., 85.5, 1., 90.5, 1.,
    94.5, 1., 96.5, 1.,
    -1 };

  const char hist_name[] = "h_bergbusch_thesis_1995_Ar";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Bergbusch MS thesis (1995), figure 5.8, Ar");

  Data->fResp = triumf_response_1998;

  return 0;
}
