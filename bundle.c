
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bundle.h"

// Resource usage
double getutime () 
{
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);

	return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1e-6;
}


// BLAS and LAPACK functions

extern double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void dscal_(const int *n, const double *alpha, const double *x, const int *incx);
extern void dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *OUT_ipiv, double *b, const int *ldb, int *OUT_info);

double ddot (int n, const double *x, const double *y)
{
	int one = 1;
	return ddot_(&n, x, &one, y, &one);
}

void daxpy (int n, double alpha, const double *x, double *y)
{
	int one = 1;
	daxpy_(&n, &alpha, x, &one, y, &one);
}

void dscal (int n, double alpha, double *x)
{
	int one = 1;
	dscal_(&n, &alpha, x, &one);
}

int dgesv (int n, double *a, int *piv, double *b)
{
	int info, one = 1;
	dgesv_(&n, &one, a, &n, piv, b, &n, &info);
	assert(info >= 0);
	return info;
}



int bundle_qp_solve_mask (bundle_t *bundle, uint64_t mask)
	/*
	* Attempts to solve the penalized bundle QP assuming the zero/nonzero constraint multiplier combination specified in 'mask' is correct:
	* max z - 1/2t <x^T,x>
	* s.c
	*     z - <A_i,x> <= b_i,   i = 1,...,m : mul
	*     z scalar, x n-dimensional vector
	*  
	*  This QP is convex, therefore solving it is equivalent to solving the KKT conditions system.
	*
	*  For all possible values of mask,
	*    if the i-th bit of mask is 1:
	*      => z* - <A_i,x*> = b_i.
	*    if the i-th bit of mask is 0:
	*      => mul*_i = 0.
	*   
	*  Also, supposing the KKT conditions are verified, we have, with k iterating over all 1-bits of mask,
	*    grad(z* - 1/2t <x^T,x*>) = sum over k of (mul*_k . grad(z* - <A_k,x*>))
	*      <=> [1 | -1/t x*] = sum over i of (mul*_k . [1 | -A_k])
	*      <=> t . sum over k of (mul*_k . A_k) = x*
	*          and sum over k of mul*_k         = 1.
	*  
	*  Therefore by substituting x* we obtain the system:
	*    z* - <A_i, t . sum over k of (mul*_k . A_k)> = b_i,  for all i with the i-th bit of 'mask' is 1
	*      <=>  | -t(A'.A'^T) e |  . | mul | = | b |,         with e = (1, ..., 1) and dim(e) = number of bits of mask set to 1
	*           |      e^T    0 |    |  z  |   | 1 |,         and with A' = A with the rows corresponding to the bits set to 0 in mask removed.
	*
	*  If the system is not linearly independent, then the optimal solution cannot be found with this mask.
	*  Otherwise, having mul* and z*, we verify dual feasibility, for all i set to 1 in mask:
	*     <=> mul_i* >= 0 
	*  Also, we verify dual optimality:
	*     <=> z* is minimal for all dual-feasible z* found with other masks
	*  Finally, we verify primal feasibility, for all i set to 0 in mask:
	*     <=> z* - <A_i,x*> <= b_i, and substituting x*,
	*     <=> z* - t . sum over k of mul*_k <A_i, A_k^T> <= b_i
	*
	*  If these hold then all KKT conditions are verified and z* is optimal.
	*/
{
	double z;
	int info, isize;
	uint64_t bit;
	int i, j, k;


	// get system size (in bundle->kkt_m) and active constraint indices (in bundle->kkt_i[])
	// set system rhs (in bundle->kkt_mul[])
	for (bundle->kkt_m = i = 0, bit = 1; i < bundle->m; i++, bit <<= 1) 
		if (mask & bit) {
			bundle->kkt_i[bundle->kkt_m] = i;
			bundle->kkt_mul[bundle->kkt_m] = bundle->b[i];
			++bundle->kkt_m;
		}
	bundle->kkt_mul[bundle->kkt_m] = 1.0;
	isize = 1 + bundle->kkt_m; // set matrix dimension size for LAPACK

	// set system lhs (in bundle->kkt_a[][])
	for (i = 0; i < bundle->kkt_m; i++) {
		for (j = 0; j < bundle->kkt_m; j++) 
			bundle->kkt_a[i * isize + j] = bundle->kkt_a[j * isize + i] = - (bundle->aat[bundle->kkt_i[i] * bundle->max_m + bundle->kkt_i[j]] / bundle->scale);
		bundle->kkt_a[i * isize + bundle->kkt_m] = 1.0;
		bundle->kkt_a[bundle->kkt_m * isize + i] = 1.0;
	}
	bundle->kkt_a[bundle->kkt_m * isize + bundle->kkt_m] = 0.0;

	// solve system using LAPACK, to retrieve mul* and z* (in bundle->kkt_mul[])
	info = dgesv(isize, bundle->kkt_a, bundle->ipiv, bundle->kkt_mul);

	// check if system has full rank
	if (info > 0)
		return 0;

	// system is linearly independent and has a unique solution
	z = bundle->kkt_mul[bundle->kkt_m];
	// check dual feasibility
	for (i = 0; i < bundle->kkt_m; i++) 
		if (bundle->kkt_mul[i] < 0.0)
			return 0;

	// check dual optimality
	if (z > bundle->kkt_z * 1.01)
		return 0;
	else if (z < bundle->kkt_z)
		bundle->kkt_z = z;

	// check primal feasibility
	for (k = i = 0, bit = 1; i < bundle->m; i++, bit <<= 1) {
		if (mask & bit) {
			bundle->next_b[i] = z;
		}
		else {
			bundle->next_b[i] = 0.0;
			for (j = 0; j < bundle->kkt_m; j++)
				bundle->next_b[i] += bundle->kkt_mul[j] * bundle->aat[i * bundle->max_m + bundle->kkt_i[j]];
			bundle->next_b[i] /= bundle->scale;
			bundle->next_b[i] += bundle->b[i];
			if (z > bundle->next_b[i] + bundle->epsilon)
				return 0;
		}
	}

	// the KKT conditions are all verified
	return 1;
}


int bundle_qp_solve (bundle_t* bundle, uint64_t mask, int rank, int *n_iter)
{
	/*
	 *	Enumerates all 'mask' parameters for bundle_qp_solve_mask(),
	 *	  in partial order of decreasing bit weight, 
	 *	  favoring masks with highest rank, i.e. where the most recent subgradients are active.
	 *
	 *	If the preceding iteration of the bundle was a minor step (=> bundle->most_recent_i != -1), then we know that the most recent subgradient is active:
	 *	  in this case, masks where it is set to inactive are skipped. 
	 */
	assert(rank >= 0);

	if (rank == 0) {
		--*n_iter;
		return bundle_qp_solve_mask(bundle, mask);
	}
	if (*n_iter > 0 && bundle_qp_solve(bundle, mask | (1ULL << (rank - 1)), rank - 1, n_iter))
		return 1;
	if (*n_iter > 0 && rank != bundle->most_recent_i && bundle_qp_solve(bundle, mask, rank - 1, n_iter))
		return 1;

	return 0;
}


double bundle_guess (bundle_t* bundle, int max_qp_iterations)
	/*
	 * Guesses where the x giving the estimated best z lies by solving the penalized bundle QP.
	 *
	 * If the resolution of the QP takes more than max_qp_iterations mask guesses, it is aborted,
	 *   and the bundle is repopulated as if it were aggregated in the previous bundle iteration,
	 *   in other words reduced to only 2 subgradients, the previous aggregate and the previous evaluated subgradient.
	 *   The QP is then solved anew, and completely.
	 * The procedure then computes x, the aggregate subgradient and the linearization error for this QP.
	 *
	 */
{
	int solved;
	int i, n_iter = max_qp_iterations;

	bundle->time_qp -= getutime(1);
	bundle->kkt_z = INFINITY;
	solved = bundle_qp_solve(bundle, 0ULL, bundle->m, &n_iter);

	if (!solved) {
		bundle->time_qp += getutime(1);
		++bundle->n_iterations;
		//printf("trim\nit: %i\tpen: %f\t", bundle->n_iterations, bundle->scale);
		// trim a
		if (bundle->m != 2)
			memcpy(&bundle->a[bundle->n], &bundle->a[(bundle->m - 1) * bundle->n], bundle->n * sizeof(double));
		memcpy(&bundle->a[0], bundle->agg_subg, bundle->n * sizeof(double));
		// recompute aat
		bundle->aat[0] = ddot(bundle->n, bundle->agg_subg, bundle->agg_subg);
		bundle->aat[bundle->max_m + 1] = bundle->aat[(bundle->m - 1) * bundle->max_m + bundle->m - 1];
		bundle->aat[1] = bundle->aat[bundle->max_m] = ddot(bundle->n, bundle->agg_subg, &bundle->a[bundle->n]);
		// trim b
		bundle->b[0] = (bundle->most_recent_i == -1) ? bundle->agg_next_b : bundle->agg_b;
		bundle->b[1] = bundle->b[bundle->m - 1];
		// update bundle information and solve again
		bundle->most_recent_i = -1;
		bundle->m = 2;
		bundle->time_qp -= getutime(1);
		n_iter = (max_qp_iterations < 3) ? 3 : max_qp_iterations;
		bundle->kkt_z = INFINITY;
		solved = bundle_qp_solve(bundle, 0ULL, bundle->m, &n_iter);
		assert(solved);
	}
	
	// compute the aggregate subgradient using BLAS (in bundle->agg_subg[])
	bundle->agg_b = bundle->agg_next_b = 0.0;
	bzero(bundle->agg_subg, bundle->n * sizeof(double));
	for (i = 0; i < bundle->kkt_m; i++) {
		daxpy(bundle->n, bundle->kkt_mul[i], &bundle->a[bundle->kkt_i[i] * bundle->n], bundle->agg_subg);
		bundle->agg_b += bundle->kkt_mul[i] * bundle->b[bundle->kkt_i[i]];
		bundle->agg_next_b += bundle->kkt_mul[i] * bundle->next_b[bundle->kkt_i[i]];
	}
	// compute local x 
	memcpy(bundle->kkt_x, bundle->agg_subg, bundle->n * sizeof(double));
	dscal(bundle->n, 1.0 / bundle->scale, bundle->kkt_x);
	// compute global x
	memcpy(bundle->x, bundle->best_x, bundle->n * sizeof(double));
	daxpy(bundle->n, 1.0, bundle->kkt_x, bundle->x);

	bundle->time_qp += getutime(1);
	// return guessed z
	return bundle->kkt_mul[bundle->kkt_m];
}

int bundle_update (bundle_t* bundle, double agg_subg_square)
	/*
	 * Updates the bundle after a guess, guaranteeing space for a new subgradient.
	 *
	 * Specifically, this involves updating:
	 *   - the dot product matrix bundle->aat,
	 *   - the right-hand sides bundle->b and bundle->next_b,
	 *   - the bundle size bundle->m.
	 *
	 * If under capacity:
	 *   do nothing.
	 * If saturated:
	 *   if not all subgradients are active:
	 *     remove those which aren't
	 *   otherwise:
	 *     remove all subgradients and add aggregate
	 */
{
	int i, j;

	assert(bundle->m <= bundle->max_m);

	// check if bundle is under capacity
	if (bundle->m < bundle->max_m) 
		return bundle->m++;

	// bundle is saturated, either trim it or aggregate it

	if (bundle->kkt_m < bundle->m) {
		// trim the bundle
		for (i = 0; i < bundle->kkt_m; i++) {
			for (j = 0; j < bundle->kkt_m; j++)
				bundle->aat[i * bundle->max_m + j] = bundle->aat[bundle->kkt_i[i] * bundle->max_m + bundle->kkt_i[j]];
			if (i < bundle->kkt_i[i]) {
				memcpy(&bundle->a[i * bundle->n], &bundle->a[bundle->kkt_i[i] * bundle->n], bundle->n * sizeof(double));
				bundle->b[i] = bundle->b[bundle->kkt_i[i]];
				bundle->next_b[i] = bundle->next_b[bundle->kkt_i[i]];
			}
		}
		bundle->m = bundle->kkt_m + 1;
		return bundle->kkt_m;
	}
	else {
		// replace first subgradient by aggregate 
		memcpy(&bundle->a[0], bundle->agg_subg, bundle->n * sizeof(double));
		bundle->b[0] = bundle->agg_b;
		bundle->next_b[0] = bundle->agg_next_b;
		bundle->aat[0] = agg_subg_square;
		// throw away all other subgradients
		bundle->m = 2;
		return 1;
	}
}
	



typedef enum {
	bundle_status_minor_step,
	bundle_status_major_step,
	bundle_status_cutoff,
	bundle_status_tolerably_optimal,
	bundle_status_optimal,
} bundle_status_t;


bundle_status_t bundle_step (bundle_t *bundle, int max_qp_iterations, double subgnorm_opt_tol, double linerr_opt_tol, double z_cutoff, void* data, double (*bundle_callback) (void*, double*, double*), double init_scale, double acceptable_model_exactness)
	/*
	 * Performs one iteration of the bundle method:
	 * - max_qp_iterations: maximum enumerations of mask to solve QP
	 * - subgnorm_opt_tol: threshold below which a subgradient norm are deemed to be zero,
	 * - linerr_opt_tol: threshold below which a linearization error is deemed to be zero, 
	 * - z_cutoff: threshold above which z no longer needs to be maximized any further,
	 * - data: a pointer passed on to bundle_callback,
	 * - double bundle_callback(void *data, double *x, double *subg): must evalute and return z at x, writing the subgradient at x in subg.
	 * - init_scale: initial penalty parameter
	 * - acceptable_model_exactness: if the ratio between actual and predicted improvement is above this value, then the bundle method performs a major step,
	 */
{
	int new_subg_i = -1;
	double roh;
	double previous_best_z = bundle->best_z;
	double *new_subg = NULL;
	double subg_delta, linerr_opt_delta, agg_subg_square;
	int i;


	//printf("it: %i\tpen: %f\t", bundle->n_iterations, bundle->scale);
	
	// guess where the x yielding the optimal z lies
	bundle->time_qp -= getutime(1);
	bundle->guessed_z = bundle_guess(bundle, max_qp_iterations);
	bundle->time_qp += getutime(1);
	linerr_opt_delta = bundle->agg_b - bundle->best_z;
	agg_subg_square = ddot(bundle->n, bundle->agg_subg, bundle->agg_subg);

	//printf("z-est: %f ", bundle->guessed_z);

	// update the bundle
	new_subg_i = bundle_update(bundle, agg_subg_square);

	// evaluate guessed x and subgradient at guessed x
	new_subg = &bundle->a[new_subg_i * bundle->n];
	bundle->time_callback -= getutime(1);
	bundle->actual_z = bundle_callback(data, bundle->x, new_subg);
	bundle->time_callback += getutime(1);

	// update maximums
	if (bundle->actual_z > bundle->max_z) {
		bundle->max_z = bundle->actual_z;
		memcpy(bundle->max_x, bundle->x, bundle->n * sizeof(double));
		if (bundle->max_z >= z_cutoff) {
			//printf("(%f)*\tcutoff\n", bundle->actual_z);
			return bundle_status_cutoff;
		}
	}

	//printf("(%f)%c\t", bundle->actual_z, (bundle->actual_z == bundle->max_z) ? '*' : ' ');

	// update A.A^T with new subgradient
	for (i = 0; i < bundle->m; i++) 
		bundle->aat[i * bundle->max_m + new_subg_i] = bundle->aat[new_subg_i * bundle->max_m + i] = ddot(bundle->n, &bundle->a[i * bundle->n], new_subg);

	// check for optimality, i.e. if new subgradient is a null vector
	if (bundle->aat[new_subg_i * bundle->max_m + new_subg_i] < bundle->epsilon) {
		//printf("optimal\n");
		return bundle_status_optimal;
	}
	//printf("agg.err: %f\tagg.norm: %f\t", linerr_opt_delta, sqrt(agg_subg_square));
	if (linerr_opt_delta <= linerr_opt_tol && agg_subg_square <= subgnorm_opt_tol * subgnorm_opt_tol) {
		//printf("optimal withing tolerances\n");
		return bundle_status_tolerably_optimal;
	}

	// check improvement
	roh = (bundle->actual_z - previous_best_z) / (bundle->guessed_z - previous_best_z + bundle->epsilon); 
	//printf("roh: %.2f\t", roh);
	if (roh >= acceptable_model_exactness) {
		// Major Step
		// update penalization parameter bundle->scale
		daxpy(bundle->n, -1.0, new_subg, bundle->best_subg);
		subg_delta = ddot(bundle->n, bundle->best_subg, bundle->best_subg);
		bundle->scale = (subg_delta == 0.0) ? init_scale : (1.0 / (1.0 / bundle->scale  +  ddot(bundle->n, bundle->kkt_x, bundle->best_subg) / subg_delta));
		// center the penalized QP to new maximum
		bundle->best_z = bundle->actual_z;
		memcpy(bundle->best_x, bundle->x, bundle->n * sizeof(double));
		memcpy(bundle->best_subg, new_subg, bundle->n * sizeof(double));
		// update QP rhs for new center
		memcpy(bundle->b, bundle->next_b, bundle->max_m * sizeof(double));
		bundle->b[new_subg_i] = bundle->actual_z;
		// the most recent subgradient may not be active
		bundle->most_recent_i = -1;
		//printf("major step\n");
		return bundle_status_major_step;
	}
	else {
		// Minor Step
		// update rhs for new subgradient for existing center
		bundle->b[new_subg_i] = bundle->actual_z - ddot(bundle->n, new_subg, bundle->kkt_x);
		//printf("\n");
		// the most recent subgradient is certainly active
		bundle->most_recent_i = new_subg_i;
		return bundle_status_minor_step;
	}
}



double bundle_solve (bundle_t *bundle, int max_iterations, int max_qp_iterations, double subgnorm_opt_tol, double linerr_opt_tol, double z_cutoff, 
		void *data, double (*bundle_callback) (void*, double*, double*), 
		double init_scale, double acceptable_model_exactness, double *init_x, double init_z, double *init_subg)
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
{
	int i;
	bundle_status_t status;

	// initialize the bundle 
	if (NULL == init_x) {
		for (i = 0; i < bundle->n; i++)
			bundle->best_x[i] = 0.0;
		bundle->actual_z = bundle_callback(data, bundle->best_x, &bundle->a[0]);
	}
	else {
		memcpy(bundle->best_x, init_x, bundle->n * sizeof(double));
		if (NULL == init_subg) {
			bundle->actual_z = bundle_callback(data, bundle->best_x, &bundle->a[0]);
		}
		else {
			memcpy(&bundle->a[0], init_subg, bundle->n * sizeof(double));
			bundle->actual_z = init_z;
		}
	}
	bundle->max_z = bundle->best_z = bundle->actual_z;
	memcpy(bundle->max_x, bundle->best_x, bundle->n * sizeof(double));
	memcpy(bundle->best_subg, &bundle->a[0], bundle->n * sizeof(double));
	bundle->m = 1;
	bundle->most_recent_i = 0;
	bundle->b[0] = bundle->actual_z;
	bundle->aat[0] = ddot(bundle->n, &bundle->a[0], &bundle->a[0]);
	bundle->scale = init_scale;
	bundle->time_bundle = bundle->time_qp = bundle->time_callback = 0.0;
	//printf("it: 0\tpen: %f\tmax: %f\n", bundle->scale, bundle->actual_z);

	// find a maximum within max_iterations or until greater than z_cutoff
	bundle->time_bundle -= getutime(1);
	for (bundle->n_iterations = 1; bundle->n_iterations < max_iterations; bundle->n_iterations++) {
		// perform bundle step
		status = bundle_step(bundle, max_qp_iterations, subgnorm_opt_tol, linerr_opt_tol, z_cutoff, data, bundle_callback, init_scale, acceptable_model_exactness); 

		// check for termination criteria
		if (bundle_status_cutoff == status || bundle_status_optimal == status || bundle_status_tolerably_optimal == status)
			break;
	}
	bundle->time_bundle += getutime(1);

	// terminate the bundle search
	return bundle->max_z;
}



bundle_t* bundle_create (int max_m, int n, double epsilon)
{
	bundle_t* bundle = (bundle_t*) malloc(sizeof(bundle_t));

	bundle->time_qp = bundle->time_callback = bundle->time_bundle = 0.0;

	bundle->epsilon = epsilon;
	bundle->max_m = max_m;
	bundle->m = 0;
	bundle->n = n;
	bundle->scale = 1.0;
	bundle->b = (double*) calloc(max_m, sizeof(double));
	bundle->next_b = (double*) calloc(max_m, sizeof(double));
	bundle->a = (double*) calloc(max_m * n, sizeof(double));
	bundle->aat = (double*) calloc(max_m * max_m, sizeof(double));

	bundle->most_recent_i = -1;
	bundle->max_z = bundle->best_z = -INFINITY;
	bundle->max_x = (double*) calloc(n, sizeof(double));
	bundle->best_x = (double*) calloc(n, sizeof(double));
	bundle->best_subg = (double*) calloc(n, sizeof(double));
	bundle->x = (double*) calloc(n, sizeof(double));

	bundle->agg_subg = (double*) calloc(n, sizeof(double));
	bundle->agg_b = bundle->agg_next_b = 0.0;

	bundle->kkt_m = 0;
	bundle->kkt_i = (int*) calloc(max_m, sizeof(int));
	bundle->kkt_x = (double*) calloc(n, sizeof(double));
	bundle->kkt_mul = (double*) calloc(max_m + 1, sizeof(double));
	bundle->kkt_a = (double*) calloc((max_m + 1) * (max_m + 1), sizeof(double));
	bundle->ipiv = (int*) calloc(max_m + 1, sizeof(int));

	bundle->n_iterations = 0;
	bundle->actual_z = bundle->guessed_z = -INFINITY;


	return bundle;
}


void bundle_destroy (bundle_t *bundle)
{
	free(bundle->b);
	free(bundle->next_b);
	free(bundle->a);
	free(bundle->aat);
	free(bundle->x);
	free(bundle->kkt_x);
	free(bundle->kkt_i);
	free(bundle->kkt_a);
	free(bundle->ipiv);
	free(bundle->kkt_mul);
	free(bundle->max_x);
	free(bundle->best_x);
	free(bundle->best_subg);
	free(bundle->agg_subg);
	free(bundle);
}


