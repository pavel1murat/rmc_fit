//-----------------------------------------------------------------------------
// Option = 'D' or 'd' - draw...
//-----------------------------------------------------------------------------
#include "ana/rmc_data.hh"
#include "TAxis.h"
#include "TROOT.h"

//-----------------------------------------------------------------------------
int rmc_data::GetArmstrong_1992_fig_8b_Ca40(Data_t* Data) {
  // Phys Rev C v46 i3 p1094 (1992)
  // digitized fig 6a : RMC spectrum on Al
  // data in the format (x1,y1,x2,y2,....,xn,yn,-1)
  // assuming figure plots number of events, the Y digitizations
  // were rounded down to the closest integer. Estimated errors less than 0.1

  int   nbins(100);
  float xmin(0), xmax(100);
  float bin_width = (xmax-xmin)/nbins;
  
  double data [] = {
    57.5, 112., 58.5, 121., 59.5, 112., 60.5, 107., 61.5, 125.,
    62.5, 123., 63.5, 127., 64.5, 103., 65.5, 143., 66.5, 120.,
    67.5, 112., 68.5, 126., 69.5, 106., 70.5, 102., 71.5, 101.,
    72.5, 115., 73.5,  86., 74.5,  78., 75.5,  61., 76.5,  61.,
    77.5,  58., 78.5,  42., 79.5,  31., 80.5,  36., 81.5,  23.,
    82.5,  22., 83.5,  21., 84.5,  14., 85.5,  22., 86.5,  16.,
    87.5,   6., 88.5,  14., 89.5,  10., 90.5,   5., 91.5,   5.,
    92.5,   2., 93.5,   0., 94.5,   2., 95.5,   1., 96.5,   0.,
    97.5,   3., 98.5,   3., 99.5,   2.,
    -1 };

  const char hist_name[] = "h_Armstrong_1992_fig_8b_Ca40";
  while (TObject* o = gROOT->FindObject(hist_name)) delete o;
  Data->fHist = new TH1F(hist_name,"",nbins,xmin,xmax);

  for (int np=0; data[2*np] > 0; np++) {
    int bin = (data[2*np]-xmin)/bin_width +1;
    Data->fHist->SetBinContent(bin,data[2*np+1]);
  }

  Data->fHist->SetMarkerStyle(20);
  Data->fHist->GetXaxis()->SetTitle("E_{#gamma}, MeV");
  Data->fHist->SetMarkerSize(1);
  Data->fHist->SetTitle("Phys Rev C 46 1094 (1992), figure 8b, Ca40");

  Data->fResp = triumf_response_1992;

  return 0;
}


// //-----------------------------------------------------------------------------
// //
// //-----------------------------------------------------------------------------
// TGraphErrors* Armstrong_1992_fig_6a_Al27_fit_kmax_90(const char* Option) {

//   double data[] = {
//     57.4902,62.6385, 58.1979,64.2429, 59.0688,66.2482, 59.9949,67.9879, 61.0859,69.0637,
//     62.0684,69.7391, 62.7237,70.0119, 63.2703,69.8842, 63.9812,69.6251, 64.8022,68.8346,
//     65.6232,68.0441, 66.4451,66.7212, 67.3223,64.9995, 67.9811,63.1425, 68.9690,60.6233,
//     69.6829,58.5006, 70.5613,56.1134, 71.2770,52.9258, 72.1560,50.1392, 73.0361,46.6871,
//     73.8068,43.3670, 74.6871,39.7818, 75.6776,35.6652, 76.6135,31.5481, 77.5492,27.5641,
//     78.4842,23.9794, 79.4745,19.9959, 80.5195,16.0130, 81.3985,13.2264, 82.1677,10.7050,
//     82.8815, 8.7154, 83.4304, 7.2566, 84.1436, 5.5333, 84.9655, 4.2104, 85.9516, 2.7560,
//     86.9918, 1.5684, 87.7033, 0.9099, 88.6876,0.52043, 90,      0.0,
//     -1.};
  
//   // double data [] = { // v1
//   //   57.4967, 62.5299, 58.6196,65.4436, 59.9318,67.7487, 61.1201,69.2890, 61.9969,69.7571,
//   //   63.1251, 69.9226, 64.2548,69.3249, 64.9456,68.7221, 66.0768,67.3610, 67.0843,65.0825,
//   //   68.3432, 62.5014, 69.5415,58.8510, 70.2358,56.4162, 71.1815,53.6790, 72.1282,50.4837,
//   //   73.0132, 46.6770, 73.9610,42.8710, 75.1625,37.5413, 76.2987,33.5848, 77.2477,29.1682,
//   //   78.3859, 24.1430, 79.7755,18.8154, 80.9747,14.7070, 82.2371,10.2939, 83.4963, 7.5602,
//   //   84.6307,  4.5198, 86.3276, 2.4018, 87.7092, 1.1962, 89.4040, 0.1469, 90.0   , 0.0   ,
//   //   -1.};

//   double x[100], y[100], ex[100], ey[100];
//   int np{0};
  
//   for (; data[2*np] > 0; np++) {
//     x [np] = data[2*np  ];
//     y [np] = data[2*np+1];
//     ex[np] = 0.01;
//     ey[np] = 0;
//   }

//   TGraphErrors* gr = new TGraphErrors(np,x,y,ex,ey);

//   gr->SetMarkerStyle(20);
//   gr->GetXaxis()->SetTitle("E_{#gamma}, MeV");
//   gr->SetMarkerSize(1);
//   gr->SetTitle("Armstrong PRD_1992 Al fit kMax=90 MeV");

//   TString opt(Option);
//   opt.ToUpper();
//   if (opt.Index('D') >= 0) {
//     TCanvas* c = new TCanvas("c_armstrong_1999_AL_fit_kmax_90","1992 Al fit kmax=90",1200,800);
//     gr->Draw("AP");
//   }

//   return gr;
// }
