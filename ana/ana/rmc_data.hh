//

#ifndef __triumf_data_1995_1999__
#define __triumf_data_1995_1999__

#include "TH1F.h"

class rmc_data : public TObject {
public:

  struct Data_t {
    TH1F*  fHist;                     // data histogram
    double  (*fResp)(double,double);  // response function: prob(Etrue,Emeas)
  };
  
  rmc_data();
  ~rmc_data();

  int  get_experimental_data(int Year, const char* Target, Data_t* Data);
//-----------------------------------------------------------------------------
// response parameterization
//-----------------------------------------------------------------------------
  static double  delta_function_response(double E, double Er);
  static double  triumf_response_1992   (double E, double Er);
  static double  triumf_response_1998   (double E, double Er);
//-----------------------------------------------------------------------------
// same response parameterization, different format - for TF1's and ROOT fits
//-----------------------------------------------------------------------------
  static double fun_triumf_response_1992(double* X, double* P);
  static double fun_triumf_response_1998(double* X, double* P);
//-----------------------------------------------------------------------------
// Armstrong'1999 PRC - 2 spectra - O16 and Ti
//-----------------------------------------------------------------------------
  int GetArmstrong_1999_fig_6_O16(Data_t* Data);
  int GetArmstrong_1999_fig_14_Ti(Data_t* Data);
//-----------------------------------------------------------------------------
// Gorringe'1998 PRC - 3 spectra for Ni isotopes, parameterization in the paper
// tables - 2 MeV/bin, plots - 1 MeV/bin
//-----------------------------------------------------------------------------
  int GetGorringe_1998_Ni58      (Data_t* Data);
  int GetGorringe_1998_Ni58_2MeV (Data_t* Data);
  int GetGorringe_1998_Ni60      (Data_t* Data);
  int GetGorringe_1998_Ni60_2MeV (Data_t* Data);
  int GetGorringe_1998_Ni62      (Data_t* Data);
  int GetGorringe_1998_Ni62_2MeV (Data_t* Data);
//-----------------------------------------------------------------------------
// Bergbusch MS thesis'1995 - several spectra, response, strictly speaking,
// unknown. Can only hope it should be close to the 1998 response, as the spectra
// look very different from 1992 spectra
//-----------------------------------------------------------------------------
  int GetBergbusch_thesis_1995_fig_5_2_Al (Data_t* Data);
  int GetBergbusch_thesis_1995_fig_5_2_Si (Data_t* Data); // ** not defined yet
  int GetBergbusch_thesis_1995_fig_5_3_O16(Data_t* Data);
  int GetBergbusch_thesis_1995_fig_5_3_Ti (Data_t* Data); // ** not defined yet
  int GetBergbusch_thesis_1995_fig_5_8_Ag (Data_t* Data); // ** not defined yet
  int GetBergbusch_thesis_1995_fig_5_8_Zr (Data_t* Data); // ** not defined yet
//-----------------------------------------------------------------------------
// Armstrong'1992 PRC - 5 spectra, come with the parameterization from the paper
//-----------------------------------------------------------------------------
  int GetArmstrong_1992_fig_6a_Al  (Data_t* Data);
  int GetArmstrong_1992_fig_6b_Si  (Data_t* Data);
  int GetArmstrong_1992_fig_8b_Ca40(Data_t* Data);
  int GetArmstrong_1992_fig_9a_Mo  (Data_t* Data); // ** not defined yet
  int GetArmstrong_1992_fig_9b_Sn  (Data_t* Data); // ** not defined yet
  int GetArmstrong_1992_fig_9c_Pb  (Data_t* Data); // ** not defined yet
//-----------------------------------------------------------------------------
// Armstrong'1991 PRC - Ca40 plus low statistics O16 and C12
//-----------------------------------------------------------------------------
  int GetArmstrong_1991_fig_18_Ca40(Data_t* Data);
//-----------------------------------------------------------------------------
// GetFrischknecht'1988 PRC - O16
//-----------------------------------------------------------------------------
  int GetFrischknecht_1988_O16_table_II (Data_t* Data);
//-----------------------------------------------------------------------------
// GetFrischknecht'1985 PRC - Ca40
//-----------------------------------------------------------------------------
  int GetFrischknecht_1985_Ca40_table_IV(Data_t* Data);
  
  ClassDef(rmc_data,0)
};
#endif
