//


#include "ana/rmc_data.hh"
#include "TMath.h"

ClassImp(rmc_data)

//-----------------------------------------------------------------------------
rmc_data::rmc_data() {
}

//-----------------------------------------------------------------------------
rmc_data::~rmc_data() {
}


//-----------------------------------------------------------------------------
// typical resolution scale is of the order of a few MeV,
// so use 0.01 MeV as an approximation of a perfect resolution
// function integrates to 1.
//-----------------------------------------------------------------------------
double rmc_data::delta_function_response(double E, double Er) {
  double sigma(0.01); // in MeV

  double dx = (Er-E)/sigma;
  double f = 1/(sigma*sqrt(2*M_PI))*TMath::Exp(-dx*dx/2);
  return f;
}

//-----------------------------------------------------------------------------
int rmc_data::get_response_function(const char* Response, double (**F)(double,double)) {
  int     rc (0);
  TString response(Response);
  
  if      (response == "1992") *F = triumf_response_1992;
  else if (response == "1998") *F = triumf_response_1998;
  else    rc = -1;

  return rc;
}

//-----------------------------------------------------------------------------
// get experimental data
//-----------------------------------------------------------------------------
int rmc_data::get_experimental_data(int Year, const char* Target, const char* Response, rmc_data::Data_t* Data) {

  TString target(Target), response(Response);
;

  if (Year == 1995) {
    if      (target == "Al" ) GetBergbusch_thesis_1995_fig_5_2_Al (Data);
    else if (target == "O16") GetBergbusch_thesis_1995_fig_5_3_O16(Data);
  }
  else if (Year == 1998) {
    if      (target == "Ni58"     ) GetGorringe_1998_Ni58     (Data);
    if      (target == "Ni58_2MeV") GetGorringe_1998_Ni58_2MeV(Data);
    else if (target == "Ni60"     ) GetGorringe_1998_Ni60     (Data);
    if      (target == "Ni60_2MeV") GetGorringe_1998_Ni60_2MeV(Data);
    else if (target == "Ni62"     ) GetGorringe_1998_Ni62     (Data);
    if      (target == "Ni62_2MeV") GetGorringe_1998_Ni62_2MeV(Data);
  }
  else if (Year == 1999) {
    if (target == "O16") GetArmstrong_1999_fig_6_O16(Data);
    if (target == "Ti" ) GetArmstrong_1999_fig_14_Ti(Data);
  }
  else if (Year == 1988) {
    if (target == "O16") GetFrischknecht_1988_O16_table_II(Data);
  }
  else if (Year == 1992) {
    if      (target == "Al"  ) GetArmstrong_1992_fig_6a_Al  (Data);
    else if (target == "Si"  ) GetArmstrong_1992_fig_6b_Si  (Data);
    else if (target == "Ca40") GetArmstrong_1992_fig_8b_Ca40(Data);
    else if (target == "Mo"  ) GetArmstrong_1992_fig_9a_Mo  (Data);
    else if (target == "Sn"  ) GetArmstrong_1992_fig_9b_Sn  (Data);
    else if (target == "Pb"  ) GetArmstrong_1992_fig_9c_Pb  (Data);
  }
  else if (Year == 1991) {
    if      (target == "Ca40") GetArmstrong_1991_fig_18_Ca40(Data);
  }
//-----------------------------------------------------------------------------
// force response, if it was requested explicitly
//-----------------------------------------------------------------------------
  if      (response == "1992") Data->fResp = triumf_response_1992;
  else if (response == "1998") Data->fResp = triumf_response_1998;
  
  return 0;
}

