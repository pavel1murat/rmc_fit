//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"

int rmc_data::GetArmstrong_1992_fig_9c_Pb(Data_t* Data) {
  // Phys Rev C v46 i3 p1094 (1992)
  // digitized fig 6a : RMC spectrum on Pb
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);

  float bin_width = (xmax-xmin)/nbins;

  double data [] = {
    57.5, 75., 58.5, 78., 59.5, 83., 60.5, 81., 61.5, 71.,
    62.5, 66., 63.5, 60., 64.5, 72., 65.5, 62., 66.5, 52.,
    67.5, 46., 68.5, 40., 69.5, 43., 70.5, 25., 71.5, 29.,
    72.5, 20., 73.5, 15., 74.5, 13., 75.5, 15., 76.5, 10.,
    77.5,  6., 78.5,  8., 79.5,  2., 80.5,  2., 81.5,  4.,
    82.5,  4., 83.5,  1., 84.5,  1., 85.5,  1., 86.5,  0.,
    87.5,  1., 88.5,  1., 89.5,  0., 90.5,  0., 91.5,  0.,
    92.5,  0., 93.5,  0., 94.5,  0., 95.5,  1., 96.5,  0.,
    97.5,  0., 98.5,  1., 99.5,  0.,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Armstrong_1992_fig_9c_Pb","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Phys Rev C 46 1094 (1992), figure 9c, Pb");

  Data->fResp = triumf_response_1992;

  return 0;
}
