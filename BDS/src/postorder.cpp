#include <Rcpp.h>
using namespace Rcpp;

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
int rcpp_postorder(NumericVector lambda, NumericVector mu, double eta, IntegerVector po, IntegerMatrix edge, NumericVector branch_lengths, int rootnode) {
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
      // Rprintf("%f %f %f \n", D_start[0], D_start[1], D_start[2]);
    }

  // solve the ODE here
  // tmp_res <- backwards(lambda, mu, eta, tstart = 0, tend = bl, E_start, D_start, ntimes = ntimes)


    Rprintf("Anc: %i, \t dec: %i \t bl: %f \n", anc, dec, bl);
  }

  return 0;
}

