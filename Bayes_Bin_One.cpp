#include <Rcpp.h>
using namespace Rcpp;

///////////////////////////////////////////////////////////
// One-sample

// predictive probability afer interim given beta prior
// which is a beta-binom
double beta_binom(int s, int t, double a, double b, int N, int n){
  double num = R::beta(s + t + a, N - s - t + b);
  double den = (N - n - t)*R::beta(t + 1,N - n- t)*R::beta(s + a, n - s + b);
  return(num/den);
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

// Conditional power
// [[Rcpp::export]]
double Pcond_Bin_One_cpp(double s, double n, double N, double alpha,
                         double p1, double p0){
  double pbar = (p1 + p0)/2;
  double s2 = pbar*(1 - pbar);
  double phat = s/n;
  double theta = p1 - p0;
  double In = n/s2;
  double IN = N/s2;
  double Zn = (phat - p0)/std::sqrt(phat*(1 - phat)/n);
  double Za = R::qnorm(1-alpha, 0, 1, 1, 0);
  double num = Zn*std::sqrt(In) + theta*(IN - In) - Za*std::sqrt(IN);
  double den = std::sqrt(IN - In);
  return(R::pnorm(num/den, 0, 1, 1, 0));
}

// Predictive power
// [[Rcpp::export]]
double Ppred_Bin_One_cpp(double a, double b, double s, double n, 
                         double N, double alpha, double p1, double p0){
  double res=0;
  double Za = R::qnorm(1-alpha, 0, 1, 1, 0);
  
  for (int t = 0; t < N - n; ++t){
    double x = s + t;
    double p = beta_binom(s, t, a, b, N, n);
    double phat = x/N;
    double ZN = (phat - p0)/std::sqrt(phat*(1 - phat)/N);
    if (ZN > Za) res += p;
  }
  return(res);
}

// Predictive probability
// [[Rcpp::export]]
double Pbayes_Bin_One_cpp(double a, double b, double s, double n, 
                          double N, double neta, double p0, double delta){
  double res=0;
  
  for (int t = 0; t < N - n; ++t){
    double x = s + t;
    double p = beta_binom(s, t, a, b, N, n);
    double post = 1 - R::pbeta(p0 + delta, a+x, b+N-x, 1, 0);
    if (post > neta) res += p;
  }
  return(res);
}

// [[Rcpp::export]]
// Interim monitoring
List IA_Bin_One_cpp(double n, double s, double N, double neta, 
                    double p1, double p0, double alpha, double delta, 
                    NumericVector a, NumericVector b){
  
  int nr = a.length();
  NumericVector prob_BP(nr);
  NumericVector prob_PP(nr);
  
  double prob_CP = Pcond_Bin_One_cpp(s, n, N, alpha, p1, p0);
  
  for (int i=0; i < nr; i++) {
    prob_BP[i] = Pbayes_Bin_One_cpp(a[i], b[i], s, n, N, neta, p0, delta);
    
    prob_PP[i] = Ppred_Bin_One_cpp(a[i], b[i], s, n, N, alpha, p1, p0);
  }
  
  List output = List::create(Named("PredProb") = prob_BP, _["PredPower"] = prob_PP, 
                             _["CondPower"] = prob_CP);
  return(output);
  
}


// [[Rcpp::export]]
// Operating Characteristics
List OC_Bin_cpp_One(double nsim, double as, double bs, double p1, double p0, 
                    NumericVector n, double delta, double neta, double alpha, double tau,
                    NumericVector a, NumericVector b){
  
  int nr = a.length();
  int nc = n.length();
  double N = n[nc-1];
  
  NumericMatrix res_BP(nr, nc+2);
  NumericMatrix res_PP(nr, nc+2);
  NumericMatrix res_CP(nr, nc+2);
  
  //Rcout << "Interim+final: " << nc << "\n";
  
  double z_final;
  z_final = R::qnorm(1-alpha, 0, 1, 1, 0);
  
  for (int k = 0; k < nsim; ++k){
    double prop1 = rbeta(1, as, bs)[0];
    NumericVector x = rbinom(N, 1, prop1);
    NumericVector r = cumsum_bin_cpp(x, n);
    //Rcout << "Current n: " << n << "\n";  
    //Rcout << "Current x1: " << x1 << "\n"; 
    //Rcout << "Current x2: " << x2 << "\n";
    //Rcout << "Current s1: " << r1 << "\n"; 
    //Rcout << "Current s2: " << r2 << "\n"; 
    
    double phat = r[nc-1]/N;
    double ZN = (phat - p0)/std::sqrt(phat*(1 - phat)/N);
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
        
        double ni = n[i];
        double s = r[i];
        // Rcout << "......"<< "Current n1: " << n1 << "\n";  
        // Rcout << "......"<< "Current s1: " << s1 << "\n"; 
        // Rcout << "......"<< "Current n2: " << n2 << "\n";  
        // Rcout << "......"<< "Current s2: " << s2 << "\n"; 
        if (i == 0) {
          
          prob_BP = Pbayes_Bin_One_cpp(a[j], b[j], s, ni, N, neta, p0, delta);
          
          prob_PP = Ppred_Bin_One_cpp(a[j], b[j], s, ni, N, alpha, p1, p0);
          
          prob_CP =Pcond_Bin_One_cpp(s, ni, N, alpha, p1, p0);
          
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
            prob_BP = Pbayes_Bin_One_cpp(a[j], b[j], s, ni, N, neta, p0, delta);
            
            ind_BP[j] = interim(prob_BP, tau);
            res_BP(j, i) += ind_BP[j];
            res_BP(j, nc+1) += ind_BP[j]*n[i];
          }
          
          if(ind_PP[j] == 0) {
            prob_PP = Ppred_Bin_One_cpp(a[j], b[j], s, ni, N, alpha, p1, p0);
            
            ind_PP[j] = interim(prob_PP, tau);
            res_PP(j, i) += ind_PP[j];
            res_PP(j, nc+1) += ind_PP[j]*n[i];
          }
          
          if(ind_CP[j] == 0){
            prob_CP =Pcond_Bin_One_cpp(s, ni, N, alpha, p1, p0);
            
            ind_CP[j] = interim(prob_CP, tau);
            res_CP(j, i) += ind_CP[j];
            res_CP(j, nc+1) += ind_CP[j]*n[i];
          }
        }
        
        else {
          if(ind_BP[j] == 1) res_BP(j, i) += 1;
          else {
            res_BP(j, nc+1) += n[i];
            prob_BP = 1 - R::pbeta(p0 + delta, a[j] + s, b[j] + N - s, 1, 0);
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
