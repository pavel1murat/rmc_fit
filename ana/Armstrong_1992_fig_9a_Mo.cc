//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"

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
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Armstrong_1992_fig_9a_Mo","",nbins,xmin,xmax);

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
