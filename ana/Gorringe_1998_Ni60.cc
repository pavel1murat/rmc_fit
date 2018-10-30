//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"

//-----------------------------------------------------------------------------
int rmc_data::GetGorringe_1998_Ni60(Data_t* Data) {
  // PhysRevC v58 1767  (1998)  data from table II
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1
  
  int   nbins(50);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    59.0, 190.7, 61.0, 185.9, 63.0, 188.8, 65.0, 160.2, 67.0, 136.3,
    69.0, 111.5, 71.0,  95.3, 73.0,  91.5, 75.0,  51.5, 77.0,  44.8,
    79.0,  21.9, 81.0,  12.4, 83.0,  17.2, 85.0,   8.6, 87.0,   2.9,
    89.0,   1.9, 91.0,   1.0, 93.0,   1.9, 95.0,   0.0, 97.0,   1.0,
    99.0,   0.0,
    -1 };

  if (Data->fHist) delete Data->fHist;
  Data->fHist = new TH1F("h_Gorringe_1998_Ni60","",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width + 1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("PhysRevC v58 1767  (1998)  data from table II, ^{60}Ni");

  Data->fResp = triumf_response_1998;

  return 0;
}
