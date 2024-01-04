#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::export]]
List derivatives2_cpp(NumericVector i_sr, NumericVector i_meth, NumericVector qhat, NumericVector yield, 
                               NumericVector ks, NumericVector cum_inhib_sr, NumericVector cum_inhib_meth, 
                               NumericVector xa_fresh, NumericVector decay_rate, NumericVector conc_fresh, NumericVector alpha, 
                               NumericVector COD_conv, NumericVector lnA, NumericVector E,
                               CharacterVector names_conc_fresh, CharacterVector names_alpha, CharacterVector names_COD_conv,
                               CharacterVector names_lnA, CharacterVector names_E,
                               double slurry_prod_rate, double ks_SO4, 
                               double scale_ks, double scale_yield, double scale_xa_fresh, 
                               double NH3_emis_rate_pit, double NH3_emis_rate_floor, double H2S_emis_rate,
                               
                               NumericVector xa, 
                               double slurry_mass, 
                               double xa_dead, 
                               double RFd, 
                               double iNDF, 
                               double VSd, 
                               double starch, 
                               double CP, 
                               double CF, 
                               double VFA, 
                               double urea, 
                               double TAN, 
                               double sulfate, 
                               double sulfide, 
                               double VSd_A,
                               double VSnd_A, 
                               double CH4_A_emis_cum, 
                               double NH3_emis_cum,
                               double CH4_emis_cum, 
                               double CO2_emis_cum, 
                               double COD_conv_cum,
                               double COD_conv_cum_meth, 
                               double COD_conv_cum_respir,
                               double COD_conv_cum_sr,
                               
                               double area, double sub_respir, double respiration, 
                               double temp_K, double R, double CH4_VS, double rain, 
                               double evap) {
  
  int n = qhat.size();
  NumericVector rut(n);
  
  for (int i:i_sr){
    rut[i] = ((qhat[i] * VFA / (slurry_mass) * xa[i] / (slurry_mass) / (scale_ks * ks[i] + VFA / (slurry_mass))) * (slurry_mass) *
      (sulfate / (slurry_mass)) / (ks_SO4 + sulfate / (slurry_mass))) * cum_inhib_sr[i];
  } 
  
  for(int i:i_meth){
    rut[i] = ((qhat[i] * VFA / (slurry_mass) * xa[i] / (slurry_mass)) / (scale_ks * ks[i] + VFA / (slurry_mass)) * 
      (slurry_mass)) * cum_inhib_meth[i];
  }
  
  if (any(rut < 0).is_true()) {
    stop("In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)");
  }
  
  NumericVector rutsr = rut[i_sr];
  if (rutsr.size() == 0){
    rutsr = 0;
  }
  
  double sum_decay = 0.0;
  for (int i = 0; i < decay_rate.size(); i++) {
    sum_decay += decay_rate[i] * xa[i];
  }
  
  int c_xa_dead;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "xa_dead") {
      c_xa_dead = i;
    }
  }
  
  int a_xa_dead;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "xa_dead") {
     a_xa_dead = i;
    }
  }
  
  int c_RFd;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "RFd") {
      c_RFd = i;
    }
  }
  
  int a_RFd;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "RFd") {
      a_RFd = i;
    }
  }
  
  int c_iNDF;
  
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "iNDF") {
      c_iNDF = i;
    }
  }
  
  int c_VSd;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "VSd") {
      c_VSd = i;
    }
  }
  
  int a_VSd;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "VSd") {
      a_VSd = i;
    }
  }
  
  int c_starch;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "starch") {
     c_starch = i;
    }
  }
  
  int a_starch;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "starch") {
      a_starch = i;
    }
  }
  
  int c_CP;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "CP") {
      c_CP = i;
    }
  }
  
  int a_CP;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "CP") {
      a_CP = i;
    }
  }
  
  int c_CF;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "CF") {
      c_CF = i;
    }
  }
  
  int a_CF;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "CF") {
      a_CF = i;
    }
  }
  
  int c_VFA;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "VFA") {
      c_VFA = i;
    }
  }
  
  double sum_rut = std::accumulate(rut.begin(), rut.end(), 0.0);
  
  int c_TAN;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "TAN") {
      c_TAN = i;
    }
  }
  
  int c_urea;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "urea") {
      c_urea = i;
    }
  }   
  
  int a_urea;
  for (int i = 0; i < names_alpha.size(); i++) {
    if (names_alpha[i] == "urea") {
      a_urea = i;
    }
  }
  
  int N_CP;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "N_CP") {
      N_CP = i;
    }
  }
  
  int c_SO4;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "SO4") {
      c_SO4 = i;
    }
  }
  
  int COD_conv_S;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "S") {
      COD_conv_S = i;
    }
  }
  
  double sum_rutsr = std::accumulate(rutsr.begin(), rutsr.end(), 0.0);
  
  int c_S2;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "S2") {
      c_S2 = i;
    }
  }
  
  int lnA_VSd_A;
  for (int i = 0; i < names_lnA.size(); i++) {
    if (names_lnA[i] == "VSd_A") {
      lnA_VSd_A = i;
    }
  }
  
  int E_VSd_A;
  for (int i = 0; i < names_E.size(); i++) {
    if (names_E[i] == "VSd_A") {
      E_VSd_A = i;
    }
  }
  
  int c_VSd_A;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "VSd_A") {
      c_VSd_A = i;
    }
  }
  
  int c_VSnd_A;
  for (int i = 0; i < names_conc_fresh.size(); i++) {
    if (names_conc_fresh[i] == "VSnd_A") {
      c_VSnd_A = i;
    }
  }
  
  // works until here
  NumericVector rutmeth = rut[i_meth];
  double sum_rutmeth = std::accumulate(rutmeth.begin(), rutmeth.end(), 0.0);
  
  int COD_conv_CH4;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "CH4") {
      COD_conv_CH4 = i;
    }
  }
  
  int COD_conv_CO2_anaer;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "CO2_anaer") {
      COD_conv_CO2_anaer = i;
    }
  }
  
  int COD_conv_CO2_sr;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "CO2_sr") {
      COD_conv_CO2_sr = i;
    }
  }
  
  int COD_conv_CO2_aer;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "CO2_aer") {
      COD_conv_CO2_aer = i;
    }
  } 
  
  int COD_conv_CO2_ureo;
  for (int i = 0; i < names_COD_conv.size(); i++) {
    if (names_COD_conv[i] == "CO2_ureo") {
      COD_conv_CO2_ureo = i;
    }
  }
  
  int nn = 23 + n;  
  
  NumericVector y(nn);
  NumericVector xa_n(n);
  
  for (int i = 0; i < n; i++) {
    xa_n[i] = scale_yield * yield[i] * rut[i] + scale_xa_fresh * xa_fresh[i] * slurry_prod_rate - decay_rate[i] * xa[i];
    y[i] = xa_n[i];
  }
  y[n] = slurry_prod_rate + (rain - evap) * area; // slurry_mass
  y[n+1] = slurry_prod_rate * conc_fresh[c_xa_dead] - alpha[a_xa_dead] * xa_dead + sum_decay; //xa_dead
  y[n+2] = slurry_prod_rate * conc_fresh[c_RFd] - alpha[a_RFd] * RFd - respiration * RFd/sub_respir; // RFd
  y[n+3] = slurry_prod_rate * conc_fresh[c_iNDF]; // iNDF
  y[n+4] = slurry_prod_rate * conc_fresh[c_VSd] - alpha[a_VSd] * VSd - respiration * VSd/sub_respir; //VSd
  
  y[n+5] = slurry_prod_rate * conc_fresh[c_starch] - alpha[a_starch] * starch - respiration * starch/sub_respir; //starch
  y[n+6] = slurry_prod_rate * conc_fresh[c_CP] - alpha[a_CP] * CP - respiration * CP/sub_respir;//CP
  
  y[n+7] = slurry_prod_rate * conc_fresh[c_CF] - alpha[a_CF] * CF - respiration * CF/sub_respir;//CF
  y[n+8] = alpha[a_xa_dead] * xa_dead + alpha[a_starch] * starch + alpha[a_CP] * CP + alpha[a_CF] * CF + alpha[a_RFd] * RFd + alpha[a_VSd] * VSd - sum_rut + slurry_prod_rate * conc_fresh[c_VFA]; //VFA
  
  y[n+9] = slurry_prod_rate * conc_fresh[c_urea] - alpha[a_urea] * urea; // urea
  y[n+10] = slurry_prod_rate * conc_fresh[c_TAN] + alpha[a_urea] * urea + alpha[a_CP] * CP * COD_conv[N_CP] + respiration * CP/sub_respir * COD_conv[N_CP] - NH3_emis_rate_pit - NH3_emis_rate_floor; //TAN
  y[n+11] = slurry_prod_rate * conc_fresh[c_SO4] - sum_rutsr * COD_conv[COD_conv_S]; //sulfate
  y[n+12] = slurry_prod_rate * conc_fresh[c_S2] + sum_rutsr * COD_conv[COD_conv_S] - H2S_emis_rate; //sulfide   
  
  y[n+13] = - VSd_A * ((exp(lnA[lnA_VSd_A] - E[E_VSd_A] / (R * temp_K))) * 24 / 1000 * CH4_VS) + slurry_prod_rate * conc_fresh[c_VSd_A];//VSd_A
  y[n+14] = slurry_prod_rate * conc_fresh[c_VSnd_A];  //VSnd_A          
  y[n+15] = VSd_A * (exp(lnA[lnA_VSd_A] - E[E_VSd_A] / (R * temp_K))) * 24 / 1000; //CH4_A_emis_cum
  y[n+16] = NH3_emis_rate_pit + NH3_emis_rate_floor; //NH3_emis_cum
  y[n+17] = sum_rutmeth * COD_conv[COD_conv_CH4]; //CH4_emis_cum
  y[n+18] = sum_rutmeth * COD_conv[COD_conv_CO2_anaer] + sum_rutsr * COD_conv[COD_conv_CO2_sr] + respiration * COD_conv[COD_conv_CO2_aer] + alpha[a_urea] * urea * COD_conv[COD_conv_CO2_ureo]; //CO2_emis_cum
  y[n+19] = sum_rutmeth + respiration + sum_rutsr; //COD_emis_cum
  y[n+20] = sum_rutmeth; // COD_conv_cum_meth
  y[n+21] = respiration; // COD_conv_cum_respir
  y[n+22] = sum_rutsr; // COD_conv_cum_sr                                                                                                                                        
  
  List res = Rcpp::List::create(Rcpp::Named("derivatives") = y,
                                        Rcpp::Named("rut") = rut);
    return res;
}




