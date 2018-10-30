//

#include "ana/rmc_data.hh"
#include "TAxis.h"
//-----------------------------------------------------------------------------
int rmc_data::GetArmstrong_1999_fig_14_Ti(Data_t* Data) {
  // Phys Rev C v46 i3 p1094 (1999)
  // digitized fig 14 : RMC spectrum on Ti
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(100);
  float xmin(0), xmax(100);

  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5, 112, 58.5, 140, 59.5, 104, 60.5, 113, 61.5, 147,
    62.5, 117, 63.5,  95, 64.5, 101, 65.5,  82, 66.5,  82,
    67.5,  99, 68.5,  79, 69.5,  69, 70.5,  49, 71.5,  56,
    72.5,  55, 73.5,  58, 74.5,  27, 75.5,  34, 76.5,  34,
    77.5,  27, 78.5,  20, 79.5,  22, 80.5,  12, 81.5,  16,
    82.5,   8, 83.5,   4, 84.5,   6, 85.5,   6, 86.5,   3,
    87.5,   2, 88.5,   5, 89.5,   1, 94.5,   1, 98.5,   1,
    -1 };

  if (Data->fHist) delete (Data->fHist);
  Data->fHist = new TH1F("h_Armstrong_1999_fig_14_Ti","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Phys Rev C  59 2853-2864 (1999), figure 14, Ti");

  Data->fResp = triumf_response_1998;

  return 0;
}
