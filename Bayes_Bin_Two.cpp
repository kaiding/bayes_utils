#define RCPPDIST_DONT_USE_ARMA
#include <Rcpp.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Rcpp;

//Binary Endpoints
//Bayesian predictive power (simulation-based), 
//Conditional power, and predictive power.

// predictive probability afer interim given beta prior
// which is a beta-binom
// [[Rcpp::export]]
double beta_binom(int s, int t, double a, double b, int N, int n){
  double num = R::beta(s + t + a, N - s - t + b);
  double den = (N - n - t)*R::beta(t + 1,N - n- t)*R::beta(s + a, n - s + b);
  return(num/den);
}

// posterior probability of detecting a treatment difference
// Simulation-based
// [[Rcpp::export]]
double post_Bin_MC_cpp(double a1, double b1, double a2, double b2, double delta,
                       double N1, double N2,  double x1, double x2, double nsim){
  double a1s = a1 + x1;
  double b1s = N1 - x1 + b1;
  double a2s = a2 + x2;
  double b2s = N2 - x2 +b2;
  
  NumericVector p1 = rbeta(nsim, a1s, b1s);
  NumericVector p2 = rbeta(nsim, a2s, b2s);
  NumericVector diff = p1 - p2;
  return(mean(diff > delta));
}

// Bayesian predictive probability
// Simulation-based
// [[Rcpp::export]]
double Pbayes_Bin_MC_cpp(double a1, double b1, double a2, double b2, 
                         double n1, double s1, double n2, double s2,
                         double N, double delta, double neta, double nsim){
  double res=0;
  
  for (int t1 = 0; t1 < N - n1; ++t1){
    double p1;
    double x1 = s1 + t1;
    p1 = beta_binom(s1, t1, a1, b1, N, n1);
    //Rcout << "t1: "<<t1<<" p1: "<<p1<<"\n";
    for(int t2 = 0; t2 < N - n2; ++t2){
      double p2;
      double post;
      double x2 = s2 + t2;
      post = post_Bin_MC_cpp(a1, b1, a2, b2, delta, N, N, x1, x2, nsim);
      //Rcout << "t2: "<<t2<<" p2: "<<p2<<" post: "<<post<<"\n";
      if (post > neta) {
        p2 = beta_binom(s2, t2, a2, b2, N, n2);
        res += p1*p2;
      }
    }
  }
  return(res);
}

//Binary Endpoints
//Bayesian predictive power (analytical-based), 
// Integrand used in Bayes Pred Prob
class bayes_int: public Numer::Func
{
private:
  double a1;
  double b1;
  double a2;
  double b2;
  double delta;
  
public:
  bayes_int(double a1_, double b1_, double a2_, double b2_, double delta_) : 
  a1(a1_), b1(b1_), a2(a2_), b2(b2_), delta(delta_) {}
  
  double operator()(const double& x) const
  {
    return std::pow(x, a1-1)*std::pow(1-x, b1-1)*R::pbeta(x-delta, a2, b2, 1, 0);
  }
};

double integrate_bayes(double a1, double b1, double a2, double b2, double delta){
  const double upper = 1;
  bayes_int f(a1, b1, a2, b2, delta);
  double err_est;
  int err_code;
  const double res = integrate(f, delta, upper, err_est, err_code);
  return res;
}

// posterior probability of detecting a treatment difference
// Analytical-based
// [[Rcpp::export]]
double post_Bin_Analytical_cpp(double a1, double b1, double a2, double b2, double delta,
                       double N1, double N2, double x1, double x2){
  double a1s = a1 + x1;
  double b1s = N1 - x1 + b1;
  double a2s = a2 + x2;
  double b2s = N2 - x2 +b2;
  
  double res = integrate_bayes(a1s, b1s, a2s, b2s, delta);
  double pgds = (1/(R::beta(a1s, b1s))); 
  return(pgds*res);
}

// Bayesian predictive probability
// Analytical-based
// [[Rcpp::export]]
double Pbayes_Bin_Analytical_cpp(double n1, double n2, double neta, double delta, 
                                 double N, double s1, double s2, 
                                 double a1, double a2, double b1, double b2){
  double ans = 0;
  
  for (int t1 = 0; t1 < N - n1; ++t1){
    double ptgs1 = beta_binom(s1, t1, a1, b1, N, n1);
    double pgds = (1/(R::beta(s1+t1+a1, N - s1 - t1 + b1)));
    double a1s = s1 + t1 + a1;
    double b1s = N - s1 - t1 + b1;
    //Rcout << "t1: "<<t1<<" p1: "<<ptgs1<<"\n";
    for (int t2 = 0; t2 < N - n2; ++t2) {
      double a2s = s2 + t2 + a2;
      double b2s = N - s2 - t2 + b2;
      double integral = integrate_bayes(a1s, b1s, a2s, b2s, delta);
      double pgd = pgds*integral;
      //Rcout << "t2: "<<t2<<" p2: "<<ptgs2<<" post: "<<pgd<<"\n";
      if (pgd > neta) {
        double ptgs2 = beta_binom(s2, t2, a2, b2, N, n2);
        ans += ptgs1*ptgs2;
      }
    }
  }
  return ans;
} 


//Wrapper for predictive prob using either analytical or simulation methods
// [[Rcpp::export]]
// Predictive power
double Pbayes_Bin_cpp(double a1, double b1, double a2, double b2, 
                         double n1, double s1, double n2, double s2,
                         double N, double delta, double neta, int method, double nsim){
  double tmp;
  //method = 1 for analytical
  if (method == 1) tmp = Pbayes_Bin_Analytical_cpp(n1, n2, neta, delta, N, s1, s2, a1, a2, b1, b2);

  //else (by default) for simulation
  else tmp = Pbayes_Bin_MC_cpp(a1, b1, a2, b2, n1, s1, n2, s2, N, delta, neta, nsim);
  
  return(tmp);
}

//Conditional Power
// [[Rcpp::export]]
double Pcond_Bin_cpp(double n1, double n2, double s1, double s2,
                     double N, double alpha, double es, double P1, double P2){
  double p1 = s1/n1;
  double p2 = s2/n2;
  double n = (n1 + n2)/2;
  double pbar = (s1 + s2)/(n1 + n2);
  double theta = 0;
  double PBAR = (P1+P2)/2;
  double sigma = std::sqrt(PBAR*(1-PBAR));
  if (es == 1) {
    theta = P1 - P2;
  }
  if (es == 2) {
    theta = p1 - p2;
  }
  if (es == 3) {
    theta = 0;
  }
  double Zn = (p1 - p2)/std::sqrt(pbar*(1 - pbar)*(1/n1 + 1/n2));
  return(R::pnorm((sqrt(n)*Zn+(N-n)*theta/(sigma*sqrt(2))-sqrt(N)*R::qnorm(1-alpha, 0, 1, 1, 0))/sqrt(N-n), 0, 1, 1, 0));
}

//Predictive Power
// [[Rcpp::export]]
double Ppred_Bin_cpp(double n1, double n2, double s1, double s2, double N,
                     double alpha, double a1, double b1, double a2, double b2){
  double res=0;
  double Za = R::qnorm(1-alpha, 0, 1, 1, 0);
  
  for (int t1 = 0; t1 < N - n1; ++t1){
    double x1 = s1 + t1;
    double p1 = beta_binom(s1, t1, a1, b1, N, n1);
    double p1hat = x1/N;
    //Rcout << "t1: "<<t1<<" p1: "<<p1<<"\n";
    for(int t2 = 0; t2 < N - n2; ++t2){
      double x2 = s2 + t2;
      double p2 = beta_binom(s2, t2, a2, b2, N, n2);
      double p2hat = x2/N;
      double pbar = (p1hat + p2hat)/2;
      double ZN = (p1hat - p2hat)/sqrt(2*pbar*(1-pbar)/N);
      if (ZN > Za) res += p1*p2;
      //Rcout << "t2: "<<t2<<" p2: "<<p2<<" post: "<<post<<"\n";
    }
  }
  return(res);
}


// [[Rcpp::export]]
// Interim monitoring
List IA_Bin_cpp(double n1, double n2, double s1, double s2, double N, int methodpb, double nsim_p, 
                   double delta, double neta, double es, double alpha, 
                   NumericVector a1, NumericVector b1, NumericVector a2, NumericVector b2,
                   double P1, double P2){
  
  int nr = a1.length();
  NumericVector prob_BP(nr);
  NumericVector prob_PP(nr);
  
  double prob_CP = Pcond_Bin_cpp(n1, n2, s1, s2, N, alpha, es, P1, P2);
  
  for (int i=0; i < nr; i++) {
    prob_BP[i] = Pbayes_Bin_cpp(a1[i], b1[i], a2[i], b2[i], n1, s1, n2, s2,
                                   N, delta, neta, methodpb, nsim_p);
    
    prob_PP[i] = Ppred_Bin_cpp(n1, n2, s1, s2, N,
                               alpha, a1[i], b1[i], a2[i], b2[i]);
  }
  
  List output = List::create(Named("PredProb") = prob_BP, _["PredPower"] = prob_PP, 
                             _["CondPower"] = prob_CP);
  return(output);
  
}

//Modified cumsum for selected subsets
NumericVector cumsum_bin_cpp(NumericVector sim_data, NumericVector n){
  double ns = sim_data.length();
  double tp = 0;
  double j = 0;
  NumericVector res(n.length());
  
  for (int i=0; i < ns; i++){
    tp += sim_data[i];
    if (i == n[j]-1) {
      res[j] = tp;
      j += 1;
    }
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
// Operating Characteristics
List OC_Bin_cpp(double nsim, double a1s, double b1s, double a2s, double b2s, 
                NumericVector n, int methodpb, double nsim_p, 
                double delta, double neta, double es, double alpha, double tau,
                NumericVector a1, NumericVector b1, NumericVector a2, NumericVector b2,
                double P1, double P2){
  
  int nr = a1.length();
  int nc = n.length();
  double N = n[nc-1];
  
  NumericMatrix res_BP(nr, nc+2);
  NumericMatrix res_PP(nr, nc+2);
  NumericMatrix res_CP(nr, nc+2);
  
  //Rcout << "Interim+final: " << nc << "\n";
  
  double z_final;
  z_final = R::qnorm(1-alpha, 0, 1, 1, 0);
  
  for (int k = 0; k < nsim; ++k){
    double prop1 = rbeta(1, a1s, b1s)[0];
    double prop2 = rbeta(1, a2s, b2s)[0];
    NumericVector x1 = rbinom(N, 1, prop1);
    NumericVector x2 = rbinom(N, 1, prop2);
    NumericVector r1 = cumsum_bin_cpp(x1, n);
    NumericVector r2 = cumsum_bin_cpp(x2, n);
    //Rcout << "Current n: " << n << "\n";  
    //Rcout << "Current x1: " << x1 << "\n"; 
    //Rcout << "Current x2: " << x2 << "\n";
    //Rcout << "Current s1: " << r1 << "\n"; 
    //Rcout << "Current s2: " << r2 << "\n"; 
    
    double p1 = r1[nc-1]/N;
    double p2 = r2[nc-1]/N;
    double pbar = (p1 + p2)/2;
    double ZN = (p1 - p2)/std::sqrt(2*pbar*(1 - pbar)/N);
    //Rcout << "sim: "<< k+1 << "\n";
    
    for (int j = 0; j < nr; ++j){
      NumericVector ind_BP(nr);
      NumericVector ind_PP(nr);
      NumericVector ind_CP(nr);
     // Rcout << "..."<< "scenario: "<< j+1 << "\n";

      for (int i = 0; i < nc; ++i){
        double prob_BP;
        double prob_PP;
        double prob_CP;
        
       // Rcout << "......"<<"look: "<< i+1 << "\n";
       // Rcout << "......"<< "Current Ind_BP: " << ind_BP[j] << "\n";  
       // Rcout << "......"<< "Current Ind_PP: " << ind_PP[j] << "\n";  
       // Rcout << "......"<< "Current Ind_CP: " << ind_CP[j] << "\n";  
        
        double n1 = n[i];
        double n2 = n[i];
        double s1 = r1[i];
        double s2 = r2[i];
       // Rcout << "......"<< "Current n1: " << n1 << "\n";  
       // Rcout << "......"<< "Current s1: " << s1 << "\n"; 
       // Rcout << "......"<< "Current n2: " << n2 << "\n";  
       // Rcout << "......"<< "Current s2: " << s2 << "\n"; 
        if (i == 0) {
          
          prob_BP = Pbayes_Bin_cpp(a1[j], b1[j], a2[j], b2[j], n1, s1, n2, s2,
                                   N, delta, neta, methodpb, nsim_p);
          
          prob_PP = Ppred_Bin_cpp(n1, n2, s1, s2, N,
                                  alpha, a1[j], b1[j], a2[j], b2[j]);
          
          prob_CP = Pcond_Bin_cpp(n1, n2, s1, s2, N, alpha, es, P1, P2);
          
          ind_BP[j] = interim(prob_BP, tau);
          ind_PP[j] = interim(prob_PP, tau);
          ind_CP[j] = interim(prob_CP, tau);
          
          res_BP(j, i) += ind_BP[j];
          res_PP(j, i) += ind_PP[j];
          res_CP(j, i) += ind_CP[j];
          
          res_BP(j, nc+1) += ind_BP[j]*n[i];
          res_PP(j, nc+1) += ind_PP[j]*n[i];
          res_CP(j, nc+1) += ind_CP[j]*n[i];
        }
        
        else if (i < nc - 1) {
          if(ind_BP[j] == 0)  {
            prob_BP = Pbayes_Bin_cpp(a1[j], b1[j], a2[j], b2[j], n1, s1, n2, s2,
                                     N, delta, neta, methodpb, nsim_p);
            ind_BP[j] = interim(prob_BP, tau);
            res_BP(j, i) += ind_BP[j];
            res_BP(j, nc+1) += ind_BP[j]*n[i];
          }
          
          if(ind_PP[j] == 0) {
            prob_PP = Ppred_Bin_cpp(n1, n2, s1, s2, N,
                                    alpha, a1[j], b1[j], a2[j], b2[j]);
            ind_PP[j] = interim(prob_PP, tau);
            res_PP(j, i) += ind_PP[j];
            res_PP(j, nc+1) += ind_PP[j]*n[i];
          }
          
          if(ind_CP[j] == 0){
            prob_CP = Pcond_Bin_cpp(n1, n2, s1, s2, N, alpha, es, P1, P2);
            
            ind_CP[j] = interim(prob_CP, tau);
            res_CP(j, i) += ind_CP[j];
            res_CP(j, nc+1) += ind_CP[j]*n[i];
          }
        }
        
        else {
          if(ind_BP[j] == 1) res_BP(j, i) += 1;
          else {
            res_BP(j, nc+1) += n[i];
            prob_BP = post_Bin_Analytical_cpp(a1[j], b1[j], a2[j], b2[j], delta,
                                              N, N, s1, s2);
            if (prob_BP < neta) res_BP(j, i) += 1;
          }
          
          if(ind_PP[j] == 1) res_PP(j, i) += 1;
          else {
            res_PP(j, nc+1) += n[i];
            if (ZN < z_final) res_PP(j, i) += 1;
          }
          
          if(ind_CP[j] == 1) res_CP(j, i) += 1;
          else {
            res_CP(j, nc+1) += n[i];
            if (ZN < z_final) res_CP(j, i) += 1;
          }
        }
       // Rcout << "......"<<  "PredProb/PredPow/CondPow"<< prob_BP << " " << prob_PP << " " << prob_CP << "\n"; 
       // Rcout << "......"<< "Ind_BP after this look: " << ind_BP[j] << "\n";
       // Rcout << "......"<< "Ind_PP after this look: " << ind_PP[j] << "\n";
       // Rcout << "......"<< "Ind_CP after this look: " << ind_CP[j] << "\n";
       // Rcout << "......"<< "res_BP: " << res_BP <<"\n";
       // Rcout << "......"<< "res_PP: " << res_PP <<"\n";
       // Rcout << "......"<< "res_CP: " << res_CP <<"\n";
      }
    }
  }
  res_BP = res_BP/nsim;
  res_PP = res_PP/nsim;
  res_CP = res_CP/nsim;
  res_BP(_, nc) = 1 - res_BP(_, nc-1);
  res_PP(_, nc) = 1 - res_PP(_, nc-1);
  res_CP(_, nc) = 1 - res_CP(_, nc-1);
  
  List output = List::create(Named("PredProb") = res_BP, _["PredPower"] = res_PP, 
                             _["CondPower"] = res_CP(1,_));
  return(output);
}




  

