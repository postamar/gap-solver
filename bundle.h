#ifndef BUNDLE_H
#define BUNDLE_H

#include <stdint.h>


typedef struct
{
	double epsilon;    // numerical error tolerance
	int m;             // number of subgradients in bundle
	int max_m;         // maximum bundle size
	int n;             // subgradients dimension
	double *a;         // n-by-m matrix of subgradients in bundle
	double *b;         // rhs for penalized bundle QP centered in x1
					   // for all i, if subg_i is a z0-subgradient in x0,
					   // b_i = <subg_i, x1 - x0> + z0
	double *next_b;    // rhs for subsequent bundle iteration
	double *aat;       // m-by-m matrix
	double scale;      // penalty parameter for QP

	int n_iterations;  // bundle iterations counter 
	double *x;         // x found by bundle method
	double best_z;     // best z used by bundle method
	double *best_x;    // best x used by bundle method
	double max_z;      // maximum z found by bundle method
	double *max_x;     // x yielding maximal z found by bundle method
	double guessed_z;  //
	double actual_z;   // 
	double *best_subg; // subgradient at last best_x
	int most_recent_i; // index of most recent subgradient (active in next QP)

	double *agg_subg;  // aggregated subgradient at best_x
	double agg_b;      // aggregated rhs at best_x
	double agg_next_b; // aggregated rhs at x = best_x + kkt_x

	int kkt_m;         // number of nonzero multipliers
	int *kkt_i;        // active constraint indices
	double *kkt_x;     // primal solution vector for penalized bundle QP
	double *kkt_mul;   // dual solution vector for penalized bundle QP
	int *ipiv;         // pivot array used by LAPACK
	double *kkt_a;     // matrix which is LU-factorized to find the multipliers
	double kkt_z;


	double time_bundle;    // time spent in last call to bundle_solve
	double time_qp;        // cumulative time spent solving QP in last call to bundle_solve
	double time_callback;  // cumulative time spent in callback in last call to bundle_solve

} bundle_t;


bundle_t* bundle_create (int max_m, int n, double epsilon); 
	/*
	 * Initializes a bundle struct:
	 * - max_m: bundle size (~10 recommended),
	 * - n: vector dimension
	 * - epsilon: numerical error tolerance
	 */


void bundle_destroy (bundle_t *bundle);
	/*
	 * Frees memory allocated in a bundle struct 
	 */

double bundle_solve (bundle_t *bundle, int max_iterations, int max_qp_iterations, double subgnorm_opt_tol, double linerr_opt_tol, double z_cutoff, 
		void *data, double (*bundle_callback) (void*, double*, double*), 
		double init_scale, double acceptable_model_exactness, double *init_x, double init_z, double *init_subg);
	/* 
	 * Maximizes a convex non-differentiable function using the bundle method, until reaching a value or iteration threshold (z_cutoff and max_iterations, respectively), or until fulfilling optimality criteria (see subgnorm_opt_tol and linerr_opt_tol):
	 * - max_qp_iterations: maximum enumerations of mask to solve QP
	 * - subgnorm_opt_tol: threshold below which a subgradient norm are deemed to be zero,
	 * - linerr_opt_tol: threshold below which a linearization error is deemed to be zero, 
	 * - z_cutoff: threshold above which z no longer needs to be maximized any further,
	 * - data: a pointer passed on to bundle_callback (see bundle_step),
	 * - init_scale: initial QP penalty,
	 * - acceptable_model_exactness: if the ratio between actual and predicted improvement is above this value, then the bundle method performs a major step,
	 * - init_x: initial values for best_x, may be NULL in which case 0 is used,
	 * - init_z: initial value for best_z, ignored if init_x or init_subg are NULL,
	 * - init_subg: subgradient at init_x, may be NULL.
	 */



#endif

