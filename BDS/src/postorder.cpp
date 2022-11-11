#include <Rcpp.h>
#include "ODE.h"
using namespace Rcpp;

// [[Rcpp::export]]
void ED_ode(NumericVector dy,
            NumericVector y,
            NumericVector lambda,
            NumericVector mu,
            double eta,
            int k){
  //NumericVector dy(2*k);

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
  //return dy;
}

// [[Rcpp::export]]
void EDF_ode(NumericVector dy,
             NumericVector y,
             NumericVector lambda,
             NumericVector mu,
             double eta,
             int k){
  //NumericVector dy(2*k);

  for (int i = 0; i < k; i++){
    // E
    // Extinction prob
    dy[i] = - (mu[i] - (lambda[i] + mu[i] + eta) * y[i] + lambda[i] * pow(y[i], 2.0));

    // D
    // Branch prob
    dy[k+i] = (-1) * (- (lambda[i] + mu[i] + eta) * y[k+i] + 2 * lambda[i] * y[k+i] * y[i]);

    // F
    // Branch prob conditional on ancestral state
    dy[2*k + i] = (-1) * (- (lambda[i] + mu[i] + eta) * y[2*k+i] + 2 * lambda[i] * y[2*k+i] * y[i]);

    for (int j = 0; j < k; j++){
      double tmpE = 0;
      double tmpD = 0;
      double tmpF = 0;
      if (j != i){
        tmpE += y[j];
        tmpD += y[k+j];
        tmpF += y[2*k+j];
      }
      // E
      dy[i] += (eta / (k-1)) * tmpE;
      // D
      dy[k+i] += (eta / (k-1)) * tmpD;
      // F
      dy[2*k+i] -= (eta / (k-1)) * tmpF;
    }
  }
  //return dy;
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

  NumericVector dy(2*k);

  for (int m = 1; m < nsteps; m++){
    // Euler's formula:
    // y_{n+1} = y_n + dt * f(y_n, t_n)
    NumericVector y_current = y(m-1, _);
    // dy =
    ED_ode(dy, y_current, lambda, mu, eta, k); // updates dy
    y(m, _) = y_current + dt * dy;
  }
  return y;
}

// Solving the ODE forwards using Euler's method, fixed time step
// [[Rcpp::export]]
NumericMatrix rcpp_forwards(NumericVector lambda,
                             NumericVector mu,
                             double eta,
                             NumericVector u0,
                             double bl,
                             int nsteps){
  // double t = 0.0;
  int k = lambda.length();
  double dt = - bl / nsteps; // dt is negative, since time is reversed

  NumericMatrix y(nsteps, k*2);

  // E(t=0) and D(t=0)
  for (int i = 0; i < k*2; i++){
    y(0, i) = u0[i];
  }

  NumericVector dy(3*k);

  for (int m = 1; m < nsteps; m++){
    // Euler's formula:
    // y_{n+1} = y_n + dt * f(y_n, t_n)
    NumericVector y_current = y(m-1, _);
    ED_ode(dy, y_current, lambda, mu, eta, k); // updates dy
    y(m, _) = y_current + dt * dy;
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
int rcpp_get_ancestor(IntegerMatrix edge, int node_idx){
  int nrow = edge.nrow();

  int res;

  for (int i = 0; i < nrow; i++){
    int current_anc = edge(i, 1);

    if (current_anc == node_idx){
      res = i;
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

    D_ends_unnormalized(edge_idx, _) = D_end;

    // Normalize D
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
  res["root_probs"] = root_probs;
  res["D_ends"] = D_ends;
  res["D_ends_unnormalized"] = D_ends_unnormalized;
  res["E_ends"] = E_ends;

  return res;
}

// [[Rcpp::export]]
List rcpp_preorder(NumericVector lambda,
                   NumericVector mu,
                   double eta,
                   IntegerVector po,
                   IntegerMatrix edge,
                   NumericVector branch_lengths,
                   NumericVector root_probs,
                   NumericMatrix E_ends,
                   NumericMatrix D_ends,
                   NumericMatrix D_ends_unnormalized,
                   int rootnode,
                   int nsteps) {
  int nrow = edge.nrow();
  int k = lambda.length();

  int left_root_edge = rcpp_get_descendants(edge, rootnode+1)[0];
  NumericVector E_start(k);
  for (int i = 0; i < k; i++){
    E_start[i] = E_ends(left_root_edge, i);
  }
  // ############################
  //
  //        PREORDER
  //
  // ############################
  NumericMatrix F_ends(nrow, k);

  for (int idx = edge.nrow(); idx > 0; idx--){
    int edge_idx = po[idx-1]-1;

    //Rprintf("edge idx: %i \n", edge_idx);
    int parent = edge(edge_idx, 0);
    int child = edge(edge_idx, 1);
    double bl = branch_lengths[edge_idx];
    // Rprintf("parent: %i \t child: %i \n", parent, child);

    NumericVector F_start(k);

    if (parent == rootnode+1){
      IntegerVector root_children = rcpp_get_descendants(edge, rootnode+1);

      int other_child;
      for (int i = 0; i < 2; i++){
        int root_child = root_children[i];

        if (root_child != child){
          other_child = root_children[i];
        }
        F_start = D_ends(other_child, _) * lambda;
      }
    }else{
      int parent_edge = rcpp_get_ancestor(edge, parent);

      IntegerVector children = rcpp_get_descendants(edge, parent)-1;

      int other_child;
      for (int i = 0; i < 2; i++){
        if (children[i] != child){
          other_child = children[i];
        }
      }
      F_start = F_ends(parent_edge, _) * lambda * D_ends(other_child, _);
    }


    NumericVector E_start(k);
    E_start = E_ends(edge_idx, _);
    NumericVector D_start(k);
    D_start = D_ends_unnormalized(edge_idx, _);

    NumericVector u0(2*k);
    for (int i = 0; i < k; i++){
      u0[i] = E_start[i];
      u0[i+k] = D_start[i];
      u0[i+2*k] = F_start[i];
    }
    // Solve the ODE
    NumericMatrix y = rcpp_forwards(lambda, mu, eta, u0, bl, nsteps);
    NumericVector yfinal = y(nsteps-1, _);
    NumericVector F_end(k);
    for (int i = 0; i < k; i++){
      F_end[i] = yfinal[2*k + i];
    }

    // Normalize F
    double fsum = 0.0;
    for (int i = 0; i < k; i++){
      fsum += F_end[i];
    }
    for (int i = 0; i < k; i++){
      F_end[i] = F_end[i] / fsum;
    }
    // Store result
    F_ends(edge_idx, _) = F_end;
  }

  List res;
  res["F_ends"] = F_ends;

  return res;
}


