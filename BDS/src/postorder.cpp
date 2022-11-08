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

NumericVector ED_ode(NumericVector y,
                     NumericVector lambda,
                     NumericVector mu,
                     double eta,
                     int k){
  NumericVector dy(2*k);

  for (int i = 0; i < k; i++){
    // Extinction prob
    dy[i] = mu[i] - (lambda[i] + mu[i] + eta) * y[i] + lambda[i] * pow(y[i], 2.0);

    // Branch prob
    //dD[i] = - (λ[i] + μ[i] + η) * D[i] + 2*λ[i]*D[i]*Et[i] + (η/(K-1)) * sum(D[j] for j in i_not_j)
    dy[k+i] = mu[i] - (lambda[i] + mu[i] + eta) * y[k+i] + 2 * lambda[i] * y[k+i] * y[i];

    for (int j = 0; j < k; j++){
      double tmpE = 0;
      double tmpD = 0;
      if (j != i){
        tmpE += y[j];
        tmpD += y[k+j];
      }
      // E
      dy[i] += (eta / (k-1)) * tmpE;
      // D
      dy[k+i] += (eta / (k-1)) * tmpD;
    }
  }
  return dy;
}

// Solving the ODE using Euler's method, fixed time step
// [[Rcpp::export]]
NumericMatrix rcpp_backwards(NumericVector lambda,
                              NumericVector mu,
                              double eta,
                              NumericVector u0,
                              double bl,
                              int nsteps){
  double t = 0.0;
  int k = lambda.length();
  double dt = bl / nsteps;

  NumericMatrix y(nsteps, k*2);

  // E(t=0) and D(t=0)
  for (int i = 0; i < k*2; i++){
    y(0, i) = u0[i];
  }

  NumericVector dy;

  for (int i = 1; i < nsteps; i++){
    // Euler's formula:
    // y_{n+1} = y_n + dt * f(y_n, t_n)
    NumericVector y_current = y(i-1, _);
    // dy =
    y(i, _) = y_current + dt * ED_ode(y_current, lambda, mu, eta, k);
  }
  return y;
}

// [[Rcpp::export]]
IntegerVector rcpp_get_descendants(IntegerMatrix edge, int node_idx){
  int nrow = edge.nrow();

  IntegerVector res;

  for (int i = 0; i < nrow; i++){
    int current_anc = edge(i, 0);

    if (current_anc == node_idx){
      res.push_front(i);
    }
  }
  return res;
}

// [[Rcpp::export]]
int rcpp_postorder(NumericVector lambda, NumericVector mu, double eta, IntegerVector po,
                   IntegerMatrix edge, NumericVector branch_lengths, int rootnode) {
  int nrow = edge.nrow();
  int nsols = 100; // solutions per branch

  int k = lambda.length();
  //int rootnode =

    // list
  //
  // forward_results <- list()
  // backward_results <- list()

  NumericMatrix D_ends(nrow, k);
  NumericMatrix D_ends_unnormalized(nrow, k);
  NumericMatrix E_ends(nrow, k);

  NumericVector sf(nrow);




  for (int edge_idx = 0; edge_idx < nrow; edge_idx++){
    int anc = edge(edge_idx, 0);
    int dec = edge(edge_idx, 1);

    double bl = branch_lengths[edge_idx];

    if (dec > rootnode){
      // This is an internal branch
      IntegerVector my_children = rcpp_get_descendants(edge, dec);

      NumericVector D_start = lambda;
      for (int j = 0; j < my_children.length(); j++){
        D_start = D_start * D_ends(j, _);
      }

      NumericVector E_start = E_ends(my_children[0], _ );
    }else{
      // This is a tip node
      NumericVector D_start;
      NumericVector E_start;
      for (int i = 0; i < k; i++){
        D_start.push_front(1.0);
        E_start.push_front(0.0);
      }
      //Rprintf("%f %f \n", D_start[0], D_start[1]);
    }

  // solve the ODE here
  // tmp_res <- backwards(lambda, mu, eta, tstart = 0, tend = bl, E_start, D_start, ntimes = ntimes)



    Rprintf("Anc: %i \t dec: %i \t bl: %f \n", anc, dec, bl);
  }

  return 0;
}

