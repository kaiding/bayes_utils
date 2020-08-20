#include <Rcpp.h>
using namespace Rcpp;

//Simulation-based functions
// Predictive distribution after interim n patients each arm
double mu_n_MC_Normal_cpp(double n, double N, double xbar, double s, double mustar, double sstar){
  double num = mustar/std::pow(sstar,2) + n*xbar/std::pow(s,2);
  double den = 1/std::pow(sstar,2) + n/std::pow(s,2);
  return(num/den);
}

double s_n_MC_Normal_cpp(double n, double N, double s, double sstar){
  double sq = std::pow(s,2)/(N-n) + 1/(1/std::pow(sstar,2) + n/std::pow(s,2));
  return(std::sqrt(sq));
}

// Predictive distribution after final sample size N each arm
double mu_N_MC_Normal_cpp(double n, double N, double xbar_n, double x_bar_N_n, 
                      double s, double mustar, double sstar){
  double num = mustar/std::pow(sstar,2) + (n*xbar_n + (N-n)*x_bar_N_n)/std::pow(s,2);
  double den = 1/std::pow(sstar,2) + N/std::pow(s,2);
  return(num/den);
}

double var_N_MC_Normal_cpp(double N, double s, double sstar){
  return(1/(1/std::pow(sstar,2) + N/std::pow(s,2)));
}

// Posterior probabiliy at final sample size N each arm
double mu_final_MC_Normal_cpp(double N, double x_bar, double s, double mustar, double sstar){
  double num = mustar/std::pow(sstar,2) + (N*x_bar)/std::pow(s,2);
  double den = 1/std::pow(sstar,2) + N/std::pow(s,2);
  return(num/den);
}

double var_final_MC_Normal_cpp(double N, double s, double sstar){
  return(1/(1/std::pow(sstar,2) + N/std::pow(s,2)));
}

double Post_Normal_cpp(double N, double x1bar, double x2bar, double s,
                          double m1star, double m2star, double sstar,
                          double delta){
  double mu1 = mu_final_MC_Normal_cpp(N, x1bar, s, m1star, sstar);
  double var1 = var_final_MC_Normal_cpp(N, s, sstar);
  double mu2 = mu_final_MC_Normal_cpp(N, x2bar, s, m2star, sstar);
  double var2 = var_final_MC_Normal_cpp(N, s, sstar);
  double res = 1 - R::pnorm(0, mu1 - mu2 - delta, sqrt(var1 + var2), 1, 0);
  return(res);
}

// Vector subset mean
NumericVector mean_sub_cpp(NumericVector n, NumericVector x){
  NumericVector tmp(n.length());
  for (int i = 0; i < n.length(); ++i){
    tmp[i] = mean(head(x, n[i]));
  }
  return(tmp);
}

// Vector subset standard deviation
NumericVector sd_sub_cpp(NumericVector n, NumericVector x){
  NumericVector tmp(n.length());
  for (int i = 0; i < n.length(); ++i){
    tmp[i] = sd(head(x, n[i]));
  }
  return(tmp);
}

// Pooled SD
double pooledsd_cpp(double sd1, double sd2, double n1, double n2){
  return(std::sqrt(((n1-1)*sd1*sd1 + (n2-1)*sd2*sd2)/(n1-1+n2-1)));
}

// Vector subset pooled SD
NumericVector psd_sub_cpp(NumericVector n, NumericVector x1, NumericVector x2){
  NumericVector res(n.length());
  for (int i = 0; i < n.length(); ++i){
    res[i] = pooledsd_cpp(x1[i], x2[i], n[i], n[i]);
  }
  return(res);
}

// Interim decision based on power/probability, eg., result=1 for futiliy stopping
double interim(double p, double tau){
  double res;
  if (p < tau) res = 1;
  else res = 0;
  return(res);
}


// [[Rcpp::export]]
// Predictive Power
double Ppred_MC_Normal_cpp(double n1, double n2, double N, 
                       double x1bar, double x2bar, double s, 
                       double m1star, double m2star, double sstar, 
                       double alpha, double nsim){
  double tmp = 0;
  double mu1 = mu_n_MC_Normal_cpp(n1, N, x1bar, s, m1star, sstar);
  double s1 = s_n_MC_Normal_cpp(n1, N, s, sstar);
  double mu2 = mu_n_MC_Normal_cpp(n2, N, x2bar, s, m2star, sstar);
  double s2 = s_n_MC_Normal_cpp(n2, N, s, sstar);
  for (int i =0; i < nsim; ++i){
    double x1bar_N_n = rnorm(1, mu1, s1)[0];
    double x2bar_N_n = rnorm(1, mu2, s2)[0];
    double x1bar_N = (n1*x1bar + (N - n1)*x1bar_N_n)/N;
    double x2bar_N = (n2*x2bar + (N - n2)*x2bar_N_n)/N;
    double ZN = std::sqrt(N/2)*(x1bar_N - x2bar_N)/s;
    double p_tmp = 1 - R::pnorm(ZN, 0, 1, 1, 0);
    if(p_tmp < alpha) tmp += 1;
  }
  return(tmp/nsim);
}

// [[Rcpp::export]]
// Predictive Probability
double Pbayes_MC_Normal_cpp(double n1, double n2, double N, 
                        double x1bar, double x2bar, double s, 
                        double m1star, double m2star, double sstar, 
                        double neta, double delta, double nsim){
  double tmp = 0;
  double mu1 = mu_n_MC_Normal_cpp(n1, N, x1bar, s, m1star, sstar);
  double s1 = s_n_MC_Normal_cpp(n1, N, s, sstar);
  double mu2 = mu_n_MC_Normal_cpp(n2, N, x2bar, s, m2star, sstar);
  double s2 = s_n_MC_Normal_cpp(n2, N, s, sstar);
  for (int i = 0; i < nsim; ++i){
    double x1bar_N_n = rnorm(1, mu1, s1)[0];
    double x2bar_N_n = rnorm(1, mu2, s2)[0];
    double mu1_N_tmp = mu_N_MC_Normal_cpp(n1, N, x1bar, x1bar_N_n, s, m1star, sstar);
    double var1_N_tmp = var_N_MC_Normal_cpp(N, s, sstar);
    double mu2_N_tmp = mu_N_MC_Normal_cpp(n2, N, x2bar, x2bar_N_n, s, m2star, sstar);
    double var2_N_tmp = var_N_MC_Normal_cpp(N, s, sstar);
    double p_tmp = 1 - R::pnorm(0, mu1_N_tmp - mu2_N_tmp - delta, sqrt(var1_N_tmp + var2_N_tmp),
                                1, 0);
    if(p_tmp > neta) tmp += 1;
  }
  return(tmp/nsim);
}

//Analytical functions
//Parameters in predictive probability
double an_pb_NORMAL_cpp(double n, double sstar, double s){
  return(1/(1 + n*std::pow((sstar/s),2)));
}

double bn_pb_NORMAL_cpp(double n, double sstar, double s, 
                        double m1star, double m2star, double N){
  return(std::sqrt((n*N)/2)*(m1star-m2star)/(1 + n*std::pow((sstar/s),2)));
}

double cn_pb_NORMAL_cpp(double n, double sstar, double s){
  return(1 - (1/(1 + (1/n)*pow((s/sstar),2))));
}

// Parameters in predictive power
double an_pp_NORMAL_cpp(double n, double N, double sstar, double s){
  return((n-N)/N + ((N-n)/N)*an_pb_NORMAL_cpp(1/n, s, sstar));
}

double bn_pp_NORMAL_cpp(double n, double N, double sstar, double s, double m1star, double m2star){
  return(std::sqrt(n/(2*N))*(N-n)*(m1star-m2star)*an_pb_NORMAL_cpp(n, sstar, s));
}

double cn_pp_NORMAL_cpp(double n, double sstar, double s){
  return(an_pb_NORMAL_cpp(n, sstar, s)); 
}

// [[Rcpp::export]]
// Predictive power - analytical
double Ppred_Normal_Analytical_cpp(double n1, double n2, double N, 
                        double x1bar, double x2bar, double s, 
                        double m1star, double m2star, double sstar, 
                        double alpha){
  double n = (n1 + n2)/2;
  double Zn = sqrt(n/2)*(x1bar-x2bar)/s;
  double an = an_pp_NORMAL_cpp(n,N,sstar,s);
  double bn = bn_pp_NORMAL_cpp(n,N,sstar,s,m1star,m2star);
  double cn = cn_pp_NORMAL_cpp(n,sstar,s);
  double num = std::sqrt(N)*Zn*(1+an) + bn/s - std::sqrt(n)*R::qnorm(1-alpha, 0,1,1,0);
  double den = std::sqrt(((N-n)/N)*((N-n)*(1-cn)+n));
  return(R::pnorm(num/den, 0, 1, 1, 0));
}


// [[Rcpp::export]]
// Predictive probability - analytical
double Pbayes_Normal_Analytical_cpp(double n1, double n2, double N, 
                     double x1bar, double x2bar, double s, 
                     double m1star, double m2star, double sstar, 
                     double neta, double delta){
  double n = (n1 + n2)/2;
  double Zn = std::sqrt(n/2)*(x1bar-x2bar)/s;
  double an = an_pb_NORMAL_cpp(n,sstar,s);
  double bn = bn_pb_NORMAL_cpp(n,sstar,s,m1star,m2star,N);
  double Zneta = R::qnorm(neta, 0, 1, 1, 0);
  double cN = cn_pb_NORMAL_cpp(N,sstar,s);
  double aN = an_pb_NORMAL_cpp(N,sstar,s);
  double num = std::sqrt(N)*Zn*(1-an) + bn/s - 
    (delta/s)*std::sqrt((n*N)/2) - std::sqrt(n)*Zneta*(1-cN);
  double den = std::sqrt(((N-n)/N)*std::pow((1-aN),2)*((N-n)*(1-an)+n));
  return(R::pnorm(num/den, 0, 1, 1, 0));
}

//Wrapper for predictive prob/power using either analytical or simulation methods
// [[Rcpp::export]]
// Predictive power
double Ppred_Normal_cpp(double n1, double n2, double N, 
                           double x1bar, double x2bar, double s, 
                           double m1star, double m2star, double sstar, 
                           double alpha, double nsim, int method){
  double tmp;
  //method = 1 for analytical
  if (method == 1) tmp = Ppred_Normal_Analytical_cpp(n1, n2, N, x1bar, x2bar, s, 
                                                   m1star, m2star, sstar, alpha);
  //else (by default) for simulation
  else tmp = Ppred_MC_Normal_cpp(n1, n2, N, x1bar, x2bar, s, m1star, m2star, sstar, alpha, nsim);
  
  return(tmp);
}

// [[Rcpp::export]]
// Predictive probability
double Pbayes_Normal_cpp(double n1, double n2, double N, 
                            double x1bar, double x2bar, double s, 
                            double m1star, double m2star, double sstar, 
                            double neta, double delta, double nsim, int method){
  double tmp;
  //method = 1 for analytical
  if (method == 1) tmp = Pbayes_Normal_Analytical_cpp(n1, n2, N, x1bar, x2bar, s, 
                                                      m1star, m2star, sstar,
                                                      neta, delta);
  //else (by default) for simulation 
  else tmp = Pbayes_MC_Normal_cpp(n1, n2, N, x1bar, x2bar, s, m1star, m2star, sstar,
                                  neta, delta, nsim);
  
  return(tmp);
}

// [[Rcpp::export]]
// Conditional Power
double Pcond_Normal_cpp(double n, double alpha, double N, 
                        double x1bar, double x2bar, double s, double es, double MU1, double MU2, double sigma){
  double Zn = std::sqrt(n/2)*(x1bar - x2bar)/s;
  double theta = 0;
  if (es == 1){
    theta = MU1 - MU2;
  }
  if (es == 2){
    theta = x1bar - x2bar;
  }
  if (es == 3){
    theta = 0;
  }
  return(R::pnorm((sqrt(n)*Zn+(N-n)*theta/(sigma*sqrt(2))-sqrt(N)*R::qnorm(1-alpha, 0, 1, 1, 0))/sqrt(N-n), 0, 1, 1, 0));
}

// [[Rcpp::export]]
// Interim monitoring
List IA_Normal_cpp(double n1, double n2, double N, double nsim_p, int methodpb, int methodpp,
                   double delta, double neta, double alpha, 
                   double x1bar, double x2bar, double s1, double s2,
                   NumericVector m1star, NumericVector m2star, NumericVector sstar,
                   double MU1, double MU2, double sigma){
  
  int nr = m1star.length();
  NumericVector prob_BP(nr);
  NumericVector prob_PP(nr);
  NumericVector prob_CP(3);
  
  double psd = pooledsd_cpp(s1, s2, n1, n2);
  double n = (n1 + n2)/2;
  
  for (int es =1; es <4; es++){
    prob_CP[es-1] = Pcond_Normal_cpp(n, alpha, N, x1bar, x2bar, psd, es, MU1, MU2, sigma);
  }
  
  for (int i=0; i < nr; i++) {
    prob_BP[i] = Pbayes_Normal_cpp(n1, n2, N, x1bar, x2bar, psd,
                                   m1star[i], m2star[i], sstar[i],
                                   neta, delta, nsim_p, methodpb);
    
    prob_PP[i] = Ppred_Normal_cpp(n1, n2, N, x1bar, x2bar, psd,
                                  m1star[i], m2star[i], sstar[i], 
                                  alpha, nsim_p, methodpp);
  }
  
  List output = List::create(Named("PredProb") = prob_BP, _["PredPower"] = prob_PP, 
                             _["CondPower"] = prob_CP);
  return(output);
  
}



// [[Rcpp::export]]
// Operating Characteristics
List OC_Normal_cpp(double nsim, double mu1star, double mu2star, double sdstar,
                            NumericVector n, double sdcom, double nsim_p, int methodpb, int methodpp,
                            double delta, double neta, double alpha, double tau,
                            NumericVector m1star, NumericVector m2star, NumericVector sstar,
                            double MU1, double MU2, double sigma){
  int nr = m1star.length();
  int nc = n.length();
  double N = n[nc-1];
  
 // Rcout << "Interim+final: " << nc << "\n";
  
  NumericMatrix res_BP(nr, nc+2);
  NumericMatrix res_PP(nr, nc+2);
  NumericMatrix res_CP1(nr, nc+2);
  NumericMatrix res_CP2(nr, nc+2);
  NumericMatrix res_CP3(nr, nc+2);
  
  double z_final;
  z_final = R::qnorm(1-alpha, 0, 1, 1, 0);
  
  for(int k = 0; k < nsim; ++k){
    double mu1 = rnorm(1, mu1star, sdstar)[0];
    double mu2 = rnorm(1, mu2star, sdstar)[0];
    NumericVector x1 = rnorm(N, mu1, sdcom);
    NumericVector x2 = rnorm(N, mu2, sdcom);
    NumericVector x1bar = mean_sub_cpp(n, x1);
    NumericVector x2bar = mean_sub_cpp(n, x2);
    NumericVector sd1 = sd_sub_cpp(n, x1);
    NumericVector sd2 = sd_sub_cpp(n, x2);
    NumericVector psd = psd_sub_cpp(n, sd1, sd2);
    
    double ZN;
    ZN = std::sqrt(N/2)*(x1bar[nc-1] - x2bar[nc-1])/psd[nc-1];
    
    for(int j = 0; j < nr; ++j){
      NumericVector ind_BP(nr);
      NumericVector ind_PP(nr);
      NumericVector ind_CP1(nr);
      NumericVector ind_CP2(nr);
      NumericVector ind_CP3(nr);
      
      
      for(int i = 0; i < nc; ++i){
        double prob_BP;
        double prob_PP;
        double prob_CP1;
        double prob_CP2;
        double prob_CP3;
        //Rcout << "sim: "<< k+1 << "\n";
        //Rcout << "scenario: "<< j+1 << "\n";
        //Rcout << "look: "<< i+1 << "\n";
        //Rcout << "Current Ind_CP: " << ind_CP[j] << "\n";
        
        if (i == 0) {
          
          prob_BP = Pbayes_Normal_cpp(n[i], n[i], N,
                                              x1bar[i], x2bar[i], psd[i],
                                              m1star[j], m2star[j], sstar[j],
                                              neta, delta, nsim_p, methodpb);
          
          prob_PP = Ppred_Normal_cpp(n[i], n[i], N,
                                              x1bar[i], x2bar[i], psd[i],
                                              m1star[j], m2star[j], sstar[j],
                                              alpha, nsim_p, methodpp);
          
          prob_CP1 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 1, MU1, MU2, sigma);
          prob_CP2 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 2, MU1, MU2, sigma);
          prob_CP3 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 3, MU1, MU2, sigma);
          
          ind_BP[j] = interim(prob_BP, tau);
          ind_PP[j] = interim(prob_PP, tau);
          ind_CP1[j] = interim(prob_CP1, tau);
          ind_CP2[j] = interim(prob_CP2, tau);
          ind_CP3[j] = interim(prob_CP3, tau);
          
          res_BP(j, i) += ind_BP[j];
          res_PP(j, i) += ind_PP[j];
          res_CP1(j, i) += ind_CP1[j];
          res_CP2(j, i) += ind_CP2[j];
          res_CP3(j, i) += ind_CP3[j];
          
          res_BP(j, nc+1) += ind_BP[j]*n[i];
          res_PP(j, nc+1) += ind_PP[j]*n[i];
          res_CP1(j, nc+1) += ind_CP1[j]*n[i];
          res_CP2(j, nc+1) += ind_CP2[j]*n[i];
          res_CP3(j, nc+1) += ind_CP3[j]*n[i];
        }
        
        else if (i < nc - 1) {
          if(ind_BP[j] == 0)  {
            prob_BP = Pbayes_Normal_cpp(n[i], n[i], N,
                                           x1bar[i], x2bar[i], psd[i],
                                           m1star[j], m2star[j], sstar[j],
                                           neta, delta, nsim_p, methodpb);
            ind_BP[j] = interim(prob_BP, tau);
            res_BP(j, i) += ind_BP[j];
            res_BP(j, nc+1) += ind_BP[j]*n[i];
          }
          
          if(ind_PP[j] == 0) {
            prob_PP = Ppred_Normal_cpp(n[i], n[i], N,
                                          x1bar[i], x2bar[i], psd[i],
                                          m1star[j], m2star[j], sstar[j],
                                          alpha, nsim_p, methodpp);
            ind_PP[j] = interim(prob_PP, tau);
            res_PP(j, i) += ind_PP[j];
            res_PP(j, nc+1) += ind_PP[j]*n[i];
          }
          
          if(ind_CP1[j] == 0){
            prob_CP1 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 1, MU1, MU2, sigma);
            
            ind_CP1[j] = interim(prob_CP1, tau);
            res_CP1(j, i) += ind_CP1[j];
            res_CP1(j, nc+1) += ind_CP1[j]*n[i];
          }
          
          if(ind_CP2[j] == 0){
            prob_CP2 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 2, MU1, MU2, sigma);
            
            ind_CP2[j] = interim(prob_CP2, tau);
            res_CP2(j, i) += ind_CP2[j];
            res_CP2(j, nc+1) += ind_CP2[j]*n[i];
          }
          
          if(ind_CP3[j] == 0){
            prob_CP3 = Pcond_Normal_cpp(n[i], alpha, N, x1bar[i], x2bar[i], psd[i], 3, MU1, MU2, sigma);
            
            ind_CP3[j] = interim(prob_CP3, tau);
            res_CP3(j, i) += ind_CP3[j];
            res_CP3(j, nc+1) += ind_CP3[j]*n[i];
          }
        }
        
        else {
          if(ind_BP[j] == 1) res_BP(j, i) += 1;
          else {
            res_BP(j, nc+1) += n[i];
            prob_BP = Post_Normal_cpp(N, x1bar[i], x2bar[i], psd[i],
                                         m1star[j], m2star[j], sstar[j], delta);
            if (prob_BP < neta) res_BP(j, i) += 1;
          }
          
          if(ind_PP[j] == 1) res_PP(j, i) += 1;
          else {
            res_PP(j, nc+1) += n[i];
            if (ZN < z_final) res_PP(j, i) += 1;
          }
          
          if(ind_CP1[j] == 1) res_CP1(j, i) += 1;
          else {
            res_CP1(j, nc+1) += n[i];
            if (ZN < z_final) res_CP1(j, i) += 1;
          }
          
          if(ind_CP2[j] == 1) res_CP2(j, i) += 1;
          else {
            res_CP2(j, nc+1) += n[i];
            if (ZN < z_final) res_CP2(j, i) += 1;
          }
          
          if(ind_CP3[j] == 1) res_CP3(j, i) += 1;
          else {
            res_CP3(j, nc+1) += n[i];
            if (ZN < z_final) res_CP3(j, i) += 1;
          }
        }

        //Rcout << "Ind_CP after this look: " << ind_CP[j] << "\n";
        //Rcout <<  "PredProb/PredPow/CondPow"<< prob_BP << " " << prob_PP << " " << prob_CP << "\n"; 
        // Rcout << "res_CP: " << res_CP <<"\n";
      }
    }
  }
  
  res_BP = res_BP/nsim;
  res_PP = res_PP/nsim;
  res_CP1 = res_CP1/nsim;
  res_CP2 = res_CP2/nsim;
  res_CP3 = res_CP3/nsim;
  res_BP(_, nc) = 1 - res_BP(_, nc-1);
  res_PP(_, nc) = 1 - res_PP(_, nc-1);
  res_CP1(_, nc) = 1 - res_CP1(_, nc-1);
  res_CP2(_, nc) = 1 - res_CP2(_, nc-1);
  res_CP3(_, nc) = 1 - res_CP3(_, nc-1);
  
  List output = List::create(Named("PredProb") = res_BP, _["PredPower"] = res_PP, 
                             _["CondPower1"] = res_CP1(1,_), _["CondPower2"] = res_CP2(1,_),
                             _["CondPower3"] = res_CP3(1,_));
  return(output);
}










