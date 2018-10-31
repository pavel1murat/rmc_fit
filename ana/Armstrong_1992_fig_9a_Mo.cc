//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include  "TROOT.h"

int rmc_data::GetArmstrong_1992_fig_9a_Mo(Data_t* Data) {
  // Phys Rev C v46 i3 p1094 (1992)
  // digitized fig 6a : RMC spectrum on Al
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);

  float bin_width = (xmax-xmin)/nbins;

  double data [] = {
    57.5, 78., 58.5, 79., 59.5, 87., 60.5, 96., 61.5,103.,
    62.5, 83., 63.5, 78., 64.5, 84., 65.5, 78., 66.5, 79.,
    67.5, 69., 68.5, 81., 69.5, 62., 70.5, 69., 71.5, 48.,
    72.5, 50., 73.5, 54., 74.5, 35., 75.5, 32., 76.5, 27.,
    77.5, 26., 78.5, 25., 79.5, 19., 80.5, 12., 81.5, 14.,
    82.5,  9., 83.5,  9., 84.5,  9., 85.5,  2., 86.5,  2.,
    87.5,  4., 88.5,  1., 89.5,  3., 90.5,  0., 91.5,  0.,
    92.5,  1., 93.5,  0., 94.5,  0., 95.5,  0., 96.5,  0.,
    97.5,  0., 98.5,  0., 99.5,  0.,
    -1 };

  const char hist_name[] = "h_Armstrong_1992_fig_9a_Mo";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Phys Rev C 46 1094 (1992), figure 9a, Mo");

  Data->fResp = triumf_response_1992;

  return 0;
}
