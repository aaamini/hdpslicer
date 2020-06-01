#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include<random>


using namespace Rcpp;
//using namespace arma;
using arma::mat;
using arma::vec;
using arma::cube;
using arma::rowvec;
using arma::zeros;
using arma::ones;
using arma::span;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
int sample_int(int N) {
  Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[0];
} 

// [[Rcpp::export]]
arma::vec stick_break(arma::vec x) {
  int n = x.n_elem;
  arma::vec out = arma::ones(n+1);
  out(arma::span(1,n)) = arma::cumprod( arma::ones<arma::vec>(n)-x );
  return( out(arma::span(0,n-1)) %  x);
  
  //arma::vec F = arma::zeros(n);
  //arma::vec P = arma::zeros(n);
  
  //P(0) = 1;
  //F(0) = x(0);
  // for (int j = 1; j < x.n_elem; j++) {
  //   //F(j) = x(j)*P(j-1);
  //   P(j) = (1-x(j))*P(j-1);
  //   Rcout << x(j) << std::endl;
  // }
  // return( P );
}

// [[Rcpp::export]]
int find_tunc(arma::vec beta, double threshold) {
  
  int n = beta.n_elem;
  arma::vec temp = arma::zeros(n+1);
  temp(arma::span(1,n)) = cumsum(beta);
  
  int idx = 0;
  for (int j = 0; j < n+1; j++) {
    idx = n-j;
    if (temp(idx) < 1-threshold) break;
  }
  if (idx > n-1) idx--; 
  return ( idx );
   
}

// [[Rcpp::export]]
int my_sampler(arma::rowvec prob_vec) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::vector<double> weights = arma::conv_to< std::vector<double>>::from(prob_vec);
  
  std::discrete_distribution<> d(weights.begin(), weights.end());
  return( d(gen) );
}


// [[Rcpp::export]]
List hdp_infer_C_base(List y, double beta0=3, double gam0=1, 
               int ITRmax=50, int Kcap=20, int Tcap=20, int W=10) {
  
  
  // double beta0  = (opts.containsElementNamed("beta0") ?  opts["beta0"] : 3);
  // int    gam0   = (opts.containsElementNamed("gam0") ?  opts["gam0"] : 1);
  // int    ITRmax = (opts.containsElementNamed("ITRmax") ?  opts["ITRmax"] : 50);
  // int    Kcap = (opts.containsElementNamed("Kcap") ?  opts["Kcap"] : 20);
  // int    Tcap = (opts.containsElementNamed("Tcap") ?  opts["Tcap"] : 20);
  // int    W = (opts.containsElementNamed("W") ?  opts["W"] : 10);
  
  int J = y.length();
  NumericVector n(J);
  
  
  for (int j=0; j < J; j++) {
    NumericVector temp = y[j];
    n(j) = temp.length();
    
  }
  
  int nmax = max(n);
  
  // convert y from list to a matrix yb
  IntegerMatrix yb(J,nmax);
  for (int j=0; j < J; j++) {
    NumericVector temp = y[j];
    for (int i=0; i < n[j]; i++)  {
      yb(j,i) = temp(i)-1;  // fix the 0-based index issue!!!, y[j,i] will be in 0,1,...,W-1
    }
  }
  
  W = max(yb) + 1;
  
  //cube<int> Zb(J,nmax,ITRmax);
  std::vector<IntegerMatrix> zb_list;
  
  
  NumericMatrix u(J,nmax); //u = clone(y), zb = clone(y);
  IntegerMatrix tb(J,nmax), zb(J,nmax);
  IntegerMatrix T_all(J,nmax);
  IntegerMatrix K_all(J,Kcap);
  // 
  mat v(J,Tcap), gamp(J,Tcap), gam(J,Tcap);
  IntegerMatrix kb(J,Tcap);
  vec betap(Kcap), beta(Kcap);
  mat phi(W,Kcap);
  
  NumericVector my_sample(1);
  NumericVector my_probvec(Kcap);
  
  
  for (int j=0; j < J; j++) {
    for (int i=0; i < n[j]; i++){
      tb(j,i) = 2;// sample_int(10);
      
      zb(j,i) = 1;
      u(j,i) = R::runif(0,1);
    }
    for (int t=0; t < Tcap; t++){
      kb(j,t) = 1;
      v(j,t) = R::runif(0,1);
    }
    //u(j, _) = runif(nmax);
  }
  
  int itr = 0;
  bool CONVERGED = false;
  
  while (!CONVERGED && itr < ITRmax) {
    
    for (int j=0; j < J; j++) {
      for (int t = 0; t < Tcap; t++){
        
        int count1=0, count2=0;
        for (int i = 0; i < nmax; i++){
          if (tb(j,i) == t) {
            count1++;
          } else if (tb(j,i) > t) {
            count2++;
          }
        } 
        gamp(j,t) = R::rbeta(count1 + 1, count2 + gam0);
        
      } // t
     gam(j,span::all) = stick_break( gamp(j,span::all).t() ).t();
      
     for (int i=0; i < n[j]; i++){
        T_all(j,i) = find_tunc( gam(j,span::all).t(), u(j,i) );
     }
    } // j
    
    for (int k=0; k < Kcap; k++) {
      int count1=0, count2=0;
      for (int j=0; j < J; j++) {
        for (int t = 0; t < Tcap; t++){
          if (kb(j,t) == k) {
            count1++;
          } else if (kb(j,t) > k) {
            count2++;
          }
        } // t
      } // j
      betap(k) = R::rbeta(count1 + 1, count2 + beta0);
    } // k
    
    beta = stick_break( betap );
    
    for (int j=0; j < J; j++) {
      for (int t = 0; t < Tcap; t++){
        K_all(j,t) = find_tunc( beta, v(j,t) );
      }
    }
    
    
   // udpdate phi
   mat alphap(Kcap,W);
   for (int k = 0; k < Kcap; k++){
     vec dir_rand(W);
     for (int w=0; w < W; w++) {
       int count = 0;
       for (int j=0; j < J; j++) {
         for (int i=0; i < n[j]; i++) {
           if ( (zb(j,i) == k) && (yb(j,i) == w) ) {
             count++; 
           }
         } // i
       } // j
       alphap(k,w) = count + 1./W; //alpha_prior(w);
       dir_rand(w) = R::rgamma(1, alphap(k,w));
     } // w
     
     // phi is W x Kcap
     phi(span::all,k) = dir_rand / sum(dir_rand);
   } // k
   
   
   // update kb
   // w index is 1 based
   for (int j=0; j < J; j++) {
     for (int t = 0; t < Tcap; t++){
       vec nup(W);
       for (int w=0; w < W; w++) {
         int count = 0;
         for (int i = 0; i < n[j]; i++){
           if ( (tb(j,i) == t) && ( yb(j,i) == w ) )  {
             count++;  
           } 
         } //i
         nup(w) = count;   
       } //w
       int K = K_all(j,t);
       //if (K+1 > Kcap) Rcout << K;
       rowvec temp = nup.t() * log(phi.cols(0,K) + 1e-11); // (W x 1)^T (W x K) = 1 x K
       temp = temp - max(temp)*ones<rowvec>(K+1);
       rowvec prob_vec = exp(temp) + 1e-11;
       kb(j,t) = my_sampler( prob_vec );
       
       //IntegerVector my_range = seq_len(Kcap);
       //my_sample =  RcppArmadillo::sample(my_range, 1, TRUE, wrap(prob_vec));
       
     } //t
   } //j
   
   
   // update tb
   for (int j=0; j < J; j++) {
     for (int i = 0; i < n[j]; i++){
       int T = T_all(j,i);
       // if (T+1  >= Tcap) Rcout << T;
       
       rowvec prob_vec(T+1);
       for (int t = 0; t <= T; t++){
         //if (kb(j,t) >= Kcap) Rcout << kb(j,t);
         //if (yb(j,i) >= W) Rcout << yb(j,i) <<" ";
         prob_vec(t) = phi( yb(j,i), kb(j,t) );
       } //t
       tb(j,i) = my_sampler( prob_vec );
     }//i
   } //j
   
   // update u
   for (int j=0; j < J; j++) {
     for (int i = 0; i < n[j]; i++){
     u(j,i) = R::runif(0, gam(j, tb(j,i)) );
     }
   }
   
   //update v
   for (int j=0; j < J; j++) {
     for (int t = 0; t < Tcap; t++){
       v(j,t) = R::runif(0, beta( kb(j,t) ) );
     }
   }
   
   //udpdate zb
   for (int j=0; j < J; j++) {
     for (int i = 0; i < n[j]; i++){
       zb(j,i) = kb(j, tb(j,i));
     }
   }
   
   //Zb.slice(itr) = zb;
   zb_list.push_back(clone(zb));
   
   //Rcout <<  (zb_list[itr]) << "\n";
   //Rcout <<  (zb) << "\n\n";
     
   if (itr % 10 == 0) Rcout << '*';
   
    itr++;
  }
  
  // List out_list = wrap(zb_list);

  return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list );
  // return List::create(
  //      _["tb"]= tb,
  //      _["kb"] = kb,
  //      _["zb"] = zb,
  //      _["Zb"] = zb_list
  // );
  
  
  // int    n = As.n_rows;
  // arma::vec delta = arma::zeros(T);
  // 
  // arma::mat As_rescaled = (1./rho)*As, 
  //   U = arma::zeros(n,n),
  //   V = arma::zeros(n,n),
  //   X = arma::zeros(n,n),
  //   Xold = arma::zeros(n,n),
  //   Y = arma::zeros(n,n),
  //   Z = arma::zeros(n,n);
  // 
  // double alpha = (n*1.)/K;
  // 
  // 
  // int t = 0;
  // bool CONVERGED = false;
  // while (!CONVERGED && t<T) {
  //   Xold = X;
  //   X = projAXB( 0.5*(Z-U+Y-V+As_rescaled), alpha, n);
  //   Z = max(X+U, arma::zeros(n,n));
  //   Y = projToSDC(X+V);
  //   U = U+X-Z;
  //   V = V+X-Y;
  //   
  //   delta(t) = norm(X-Xold,2);
  //   CONVERGED = delta(t) < tol;
  //   
  //   if ((t+1) % report_interval == 0) {
  //     std::printf("%4d | %15e\n", t+1, delta(t));  
  //   }
  //   
  //   t++;
  // }
  // 
  // return List::create(
  //   _["X"]=X,
  //   _["delta"]=delta,
  //   _["T_term"]=t
  // );
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#res <- hdp_infer_C_base(curr_y)

#res
# temp <- seq(0.2,1,by = .1)
# stick_break_func(temp)
# stick_break(temp)

 # th_vec <- seq(1e-4,0.95,length.out = 10)
 # out = matrix(0,nrow=10,ncol=1000)
 # for (t in 1:1000) {
 #   temp <- stick_break(rbeta(10,1,1))
 #   for (tt in 1:10){
 #    th <- th_vec[tt]
 #    out[tt,t] <- find_tunc_idx(temp,th)-find_tunc(temp,th); 
 #        
 #   }
 # }
 # 
 # out
 # temp <- stick_break(rbeta(10,1,1))
 # temp
 # th <- 0.1
 # find_tunc_idx(temp,th)-  find_tunc(temp,th)
 # sum(out)
  
*/
