//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"

int rmc_data::GetArmstrong_1992_fig_9b_Sn(Data_t* Data) {
  // Phys Rev C v46 i3 p1094 (1992)
  // digitized fig 6a : RMC spectrum on Al
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);

  float bin_width = (xmax-xmin)/nbins;

  double data [] = {
    57.5,89., 58.5,89., 59.5,92., 60.5,85., 61.5,87., 62.5,78.,
    63.5,77., 64.5,92., 65.5,76., 66.5,58., 67.5,46., 68.5,61.,
    69.5,55., 70.5,42., 71.5,38., 72.5,38., 73.5,28., 74.5,24.,
    75.5,24., 76.5,20., 77.5,21., 78.5,14., 79.5,12., 80.5, 7.,
    81.5, 9., 82.5, 9., 83.5, 6., 84.5, 3., 85.5, 2., 86.5, 1.,
    87.5, 2., 88.5, 0., 89.5, 0., 90.5, 0., 91.5, 0., 92.5, 1.,
    93.5, 0., 94.5, 0., 95.5, 0., 96.5, 0., 97.5, 1., 98.5, 0.,
    99.5, 0.,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Armstrong_1992_fig_9b_Sn","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1 ;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Phys Rev C 46 1094 (1992), figure 9b, Sn");

  Data->fResp = triumf_response_1992;

  return 0;
}
