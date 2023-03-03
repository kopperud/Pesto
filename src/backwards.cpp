#include <Rcpp.h>
using namespace Rcpp;

// NumericVector rcpp_extinction_ode(NumericVector y,
//                                   NumericVector lambda,
//                                   NumericVector mu,
//                                   double eta,
//                                   int k){
//   NumericVector dy(k);
//
//   for (int i = 0; i < k; i++){
//     dy[i] = mu[i] - (lambda[i] + mu[i] + eta) * y[i] + lambda[i] * pow(y[i], 2.0);
//
//     for (int j = 0; j < k; j++){
//       double tmp = 0;
//       if (j != i){
//         // (η/(K-1)) * sum(E[j] for j in i_not_j)
//         tmp += y[j];
//       }
//       dy[i] += (eta / (k-1)) * tmp;
//     }
//   }
//   return dy;
// }



