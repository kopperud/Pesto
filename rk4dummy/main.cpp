#include <stdio.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

// #define k 4

typedef struct SSE{
  int k;
  double eta;
  std::vector<double> lambda;
  std::vector<double> mu;
} SSE;

int
func (double t, const double y[], double dy[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
//  double *mu = (double *) params;
  SSE x = *(struct SSE*)params;
  int k                      = x.k;
  double eta                 = x.eta;
  std::vector<double> lambda = x.lambda;
  std::vector<double> mu     = x.mu;

  // dE[:] = μ .- (λ .+ μ .+ η) .* E .+ λ .* E.^2 .+ (η/(K-1)) .* (sum(E) .- E) 
  // lambda = [0.1, 0.2]
  // mu = [0.05, 0.15]
  // eta = 0.05
      
  for (int i = 0; i < k; i++){
      dy[i] = mu[i] - (lambda[i] + mu[i] + eta) * y[i] + lambda[i] * y[i] * y[i];

      double tmpE = 0.0;
      for (int j = 0; j < k; j++){
          if (j != i){
              tmpE += y[j];
          }
      }
      dy[i] += eta * tmpE;
  }
  return GSL_SUCCESS;
}


int
main (void)
{
  //double mu = 10.0;
//#  gsl_odeiv2_system sys = {func, 2, &mu};
  int k = 2;
  double eta = 0.05;
  std::vector<double> lambda{0.1, 0.2};
  std::vector<double> mu{0.05, 0.15};

  SSE parms = {k, eta, lambda, mu};

  // Declare the GSL system
  gsl_odeiv2_system sys; 
  sys.function = func;
  sys.dimension = k;
  sys.params = &parms;

  gsl_odeiv2_driver * d =
//    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
                                  1e-12, 1e-12, 0.0);



  int n_reps = 1000;
  std::vector<double> times(n_reps);
  int ni = 100;
  std::vector<std::vector<double>> res(ni+1, std::vector<double>(2));

  double sum_times = 0;

  for (int m = 0; m < n_reps; m++){
    std::time_t time1 = std::clock();
    double y[2] = { 0.0, 0.0 };
    double t = 0.0, t1 = 10.0;


    for (int i = 1; i <= ni; i++)
     {
       //double ti = i * t1 / 100.0;
       double ti = i * t1 / ni;
  //     printf("t: %f \t ti: %f \n", t, ti);
       int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
 
       if (status != GSL_SUCCESS)
         {
           printf ("error, return value=%d\n", status);
           break;
         }

       for (int j = 0; j < k; j++){
         res[i][j] = y[j];
       }
    }

    std::time_t time2 = std::clock();
    double time_diff = (time2 - time1)*(1000000.0) / ((double)CLOCKS_PER_SEC);
    sum_times += time_diff;
//    times.push_back(time_diff);
  }
 //  sleep(2);
  double avg_time = sum_times / n_reps;
  printf("Clocks per second: %ld \n ", CLOCKS_PER_SEC);  
  printf("average time per %i replicates: %f µs \n", n_reps, avg_time);

  printf("final solution: ");
  for (int j = 0; j < k; j++){
    printf("y[%i] = %.15f \t", j, res[100][j]);
  }
  printf("\n");

  gsl_odeiv2_driver_free (d);
  return 0;
}
