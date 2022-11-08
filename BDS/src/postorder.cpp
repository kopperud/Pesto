#include <Rcpp.h>
//#include "./backwards.cpp"
using namespace Rcpp;

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
    dy[k+i] = - (lambda[i] + mu[i] + eta) * y[k+i] + 2 * lambda[i] * y[k+i] * y[i];

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
  // double t = 0.0;
  int k = lambda.length();
  double dt = bl / nsteps;

  NumericMatrix y(nsteps, k*2);

  // E(t=0) and D(t=0)
  for (int i = 0; i < k*2; i++){
    y(0, i) = u0[i];
  }

  // NumericVector dy;

  for (int m = 1; m < nsteps; m++){
    // Euler's formula:
    // y_{n+1} = y_n + dt * f(y_n, t_n)
    NumericVector y_current = y(m-1, _);
    // dy =
    y(m, _) = y_current + dt * ED_ode(y_current, lambda, mu, eta, k);
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
List rcpp_postorder(NumericVector lambda, NumericVector mu, double eta, IntegerVector po,
                   IntegerMatrix edge, NumericVector branch_lengths, int rootnode, int nsteps) {
  int nrow = edge.nrow();
  int k = lambda.length();

  // forward_results <- list()
  // backward_results <- list()

  NumericMatrix D_ends(nrow, k);
  NumericMatrix D_ends_unnormalized(nrow, k);
  NumericMatrix E_ends(nrow, k);

  NumericVector sf(nrow);

  for (int idx = 0; idx < nrow; idx++){
    int edge_idx = po[idx] - 1;
    int anc = edge(edge_idx, 0);
    int dec = edge(edge_idx, 1);

    double bl = branch_lengths[edge_idx];
    NumericVector E_start(k,0.0);
    NumericVector D_start(k,1.0);

    if (dec > rootnode){
      // This is an internal branch
      IntegerVector my_children = rcpp_get_descendants(edge, dec);

      //D_start = 1.0;
      NumericVector c1 = D_ends(my_children[0], _);
      NumericVector c2 = D_ends(my_children[1], _);

      for (int i = 0; i < k; i++){
        D_start[i] = c1[i] * c2[i] * lambda[i];
      }


      E_start = E_ends(my_children[0], _);
    }else{
      // This is a tip node
      NumericVector D_start(k,1.0);

      for (int i = 0; i < k; i++){
        E_start.push_front(0.0);
      }
    }

    NumericVector u0(2*k);
    for (int i = 0; i < k; i++){
      u0[i] = E_start[i];
      u0[i+k] = D_start[i];
    }

    // solve the ODE here
    NumericMatrix y = rcpp_backwards(lambda, mu, eta, u0, bl, nsteps);
    NumericVector yfinal = y(nsteps-1, _);



    E_ends(edge_idx, _) = yfinal[Range(0, k)];
    NumericVector D_end(k);
    for (int i = 0; i < k; i++){
      D_end[i] = yfinal[i+k];
    }

    double sfi = 0;
    for (int i = 0; i < k; i++){
      sfi += D_end[i];
    }
    sf[edge_idx] = log(sfi);

    for (int i = 0; i < k; i++){
      D_end[i] = D_end[i] / sfi;
    }

    D_ends(edge_idx, _) = D_end;

    //Rprintf("Anc: %i \t dec: %i \t bl: %f \n", anc, dec, bl);
  }

  // root states
  IntegerVector root_children = rcpp_get_descendants(edge, rootnode+1);

  NumericVector root_probs(k);
  for (int i = 0; i < k; i++){
    root_probs[i] = 1.0;
    for (int j = 0; j < root_children.length(); j++){
      int child = root_children[j];
      root_probs[i] *= D_ends(child, i);
    }
  }
  for (int i = 0; i < k; i++){
    root_probs[i] = root_probs[i] * lambda[i];
  }

  IntegerVector root_descendants = rcpp_get_descendants(edge, rootnode+1);
  NumericVector nonextinct(k);
  for (int i = 0; i < k; i++){
    double nx = 1.0 - E_ends(root_descendants[0], i);
    nonextinct[i] = pow(nx, 2.0);
  }

  double logL = 0;
  double tmp = 0;
  for (int i = 0; i < k; i++){
    tmp += (1.0 / k) * root_probs[i] / (nonextinct[i] * lambda[i]);
  }
  logL += log(tmp);
  logL += sum(sf);

  //Rprintf("logL: %f \n", logL);

  List res;
  res["logL"] = logL;
  res["D_ends"] = D_ends;

  return res;
}

