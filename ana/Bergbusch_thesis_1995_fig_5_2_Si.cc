///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

int rmc_data::GetBergbusch_thesis_1995_fig_5_2_Si(Data_t* Data) {
  // Bergbusch PhD thesis (1993)
  // digitized fig 6 : RMC spectrum on Si
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5,216., 58.5,225., 59.5,231., 60.5,211., 61.5,205.,
    62.5,195., 63.5,182., 64.5,193., 65.5,153., 66.5,145.,
    67.5,115., 68.5,135., 69.5,112., 70.5 ,95., 71.5,110.,
    72.5, 95., 73.5, 81., 74.5, 57., 75.5, 52., 76.5, 43.,
    77.5, 39., 78.5, 45., 79.5,  7., 80.5, 18., 81.5, 17.,
    82.5,  9., 83.5, 10., 84.5,  3., 85.5,  1., 86.5,  3.,
    87.5,  4., 88.5,  3., 96.5,  1.,
    -1 };

  const char hist_name[] = "h_bergbusch_thesis_1995_Si";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Bergbusch MS thesis (1995), figure 5.2, Si");

  Data->fResp = triumf_response_1998;

  return 0;
}
