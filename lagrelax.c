#include "lagrelax.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <float.h>



double my_cb (void *data, double *mul, double *subg)
{
	lagrelax_t *lr = (lagrelax_t*) data;

	lagrelax_eval(lr, lr->bundle_lrs->x0, lr->bundle_lrs->x1, mul, subg);

	return lr->z;
}


int lagrangian_dual_solve (lagrelax_t* lr, lagrelax_state_t* lrs, int bundle_iterations)
{
	double old_z_ub = *(lr->z_ub_ptr);

	lr->bundle_lrs = lrs;
	lr->z = lrs->z_lb = bundle_solve(lr->bundle, 
			bundle_iterations,  // max bundle iterations
			1,
			1e-5, // aggregate subg norm opt tolerance
			1e-5, // aggregate linearization error opt tolerance
			old_z_ub, // z cutoff value 
			lr, // lagrangian relaxation data, passed to my_cb
			my_cb, // my callback function for evaluating the lagrangian relaxation
			1.0, // initial penalty value
			0.8, // major step threshold
			lrs->mul, // initial point 
			0.0, // initial value
			NULL); // initial subgradient

	memcpy(lrs->mul, lr->bundle->max_x, lr->inst->n * sizeof(double));

	return (old_z_ub != *(lr->z_ub_ptr));
}


lagrelax_state_t lagrelax_state_create (gap_inst_t* inst)
{
	int i;
	lagrelax_state_t lrs;

	lrs.inst = inst;
	lrs.z_lb = -INFINITY;
	lrs.mul = (double*) calloc(inst->n, sizeof(double));
	lrs.x0 = (uint8_t**) calloc(inst->m, sizeof(uint8_t*));
	lrs.x1 = (uint8_t**) calloc(inst->m, sizeof(uint8_t*));
	for (i = 0; i < inst->m; i++) {
		lrs.x0[i] = (uint8_t*) calloc(inst->n, sizeof(uint8_t));
		lrs.x1[i] = (uint8_t*) calloc(inst->n, sizeof(uint8_t));
	}

	return lrs;
}



void lagrelax_state_copy (lagrelax_state_t* dst, const lagrelax_state_t* src)
{
	int i;

	dst->z_lb = src->z_lb;
	memcpy(dst->mul, src->mul, src->inst->n * sizeof(double));
	for (i = 0; i < src->inst->m; i++) {
		memcpy(dst->x0[i], src->x0[i], src->inst->n * sizeof(uint8_t));
		memcpy(dst->x1[i], src->x1[i], src->inst->n * sizeof(uint8_t));
	}
}



void lagrelax_state_destroy (lagrelax_state_t* lrs)
{
	int i;

	for (i = 0; i < lrs->inst->m; i++) {
		free(lrs->x0[i]);
		free(lrs->x1[i]);
	}
	free(lrs->x0);
	free(lrs->x1);
	free(lrs->mul);
}


lagrelax_t lagrelax_create (gap_inst_t* inst, double* z_ub_ptr)
{
	lagrelax_t lr;
	int i, j, maxb = -1;

	lr.inst = inst;
	lr.z_ub_ptr = z_ub_ptr;
	lr.n_eval = 0;
	lr.ks = (knapsack_t*) calloc(inst->m, sizeof(knapsack_t));
	lr.sum_x = (int*) calloc(inst->n, sizeof(int));
	lr.z = -INFINITY;
	lr.feasible_z = lr.feasible = 0;
	lr.slack = (int*) calloc(inst->m, sizeof(int));
	lr.rc = (double**) calloc(inst->m, sizeof(double*));
	lr.minknap_data = minknap_prepare(inst->n);
	lr.bundle_lrs = NULL;
	for (i = 0; i < inst->m; i++) {
		lr.ks[i] = knapsack_create(inst, i, lr.minknap_data);
		lr.rc[i] = (double*) calloc(inst->n, sizeof(double));
		if (inst->capacity[i] > maxb)
			maxb = inst->capacity[i];
	}
	++maxb;
	lr.index_list = (int*) calloc(inst->n, sizeof(int));
	lr.f_val_cache = (float**) calloc(inst->n, sizeof(float*));
	lr.g_val_cache = (float**) calloc(inst->n, sizeof(float*));
	for (j = 0; j < inst->n; j++) {
		lr.f_val_cache[j] = (float*) calloc(maxb, sizeof(float));
		lr.g_val_cache[j] = (float*) calloc(maxb, sizeof(float));
	}

	lr.bundle = bundle_create(20, inst->n, 1e-6);
	return lr;
} 

void lagrelax_destroy (lagrelax_t* lr)
{
	int i, j;

	for (i = 0; i < lr->inst->m; i++) {
		knapsack_destroy(&lr->ks[i]);
		free(lr->rc[i]);
	}
	for (j = 0; j < lr->inst->n; j++) {
		free(lr->f_val_cache[j]);
		free(lr->g_val_cache[j]);
	}
	free(lr->f_val_cache);
	free(lr->g_val_cache);
	free(lr->index_list);
	free(lr->ks);
	free(lr->sum_x);
	free(lr->rc);
	free(lr->slack);
	minknap_free(lr->minknap_data);
	bundle_destroy(lr->bundle);
}



int lagrelax_eval (lagrelax_t* lr, uint8_t** x0, uint8_t** x1, double *mul, double *subg)
{
	int i, j, z_ub_update = 0;
	int best_i;
	double r, best_r;

	++lr->n_eval;
	lr->z = 0.0;
	for (j = 0; j < lr->inst->n; j++) {
		lr->z += mul[j];
		lr->sum_x[j] = 0;
	}

	lr->feasible = 1;
	for (i = 0; i < lr->inst->m; i++) {
		knapsack_reset(&lr->ks[i], lr->inst);
		for (j = 0; j < lr->inst->n; j++) 
			lr->ks[i].p[j] -= mul[j];
		knapsack_solve(&lr->ks[i], x0[i], x1[i]);
		assert(lr->ks[i].opt != -INFINITY);
		lr->z += lr->ks[i].opt;
		for (j = 0; j < lr->inst->n; j++) {
			if (lr->ks[i].x[j]) 
				++lr->sum_x[j];
			assert(!lr->ks[i].x[j] || !x0[i][j]);
			assert(lr->ks[i].x[j] || !x1[i][j]);
		}
		if (lr->ks[i].opt == INFINITY) {
			lr->feasible = 0;
			for (j = 0; j < lr->inst->n; j++) 
				if (lr->ks[i].x[j])
					assert(x1[i][j]);
		}
		lr->slack[i] = lr->inst->capacity[i];
	}

	for (j = 0; j < lr->inst->n; j++) {
		if (lr->sum_x[j] != 1) 
			lr->feasible = 0;
	}

	if (NULL != subg) {
		for (j = 0; j < lr->inst->n; j++) 
			subg[j] = (double) (1 - lr->sum_x[j]);
	}

	lr->feasible_z = 0;
	for (j = 0; j < lr->inst->n; j++)
		if (lr->sum_x[j] == 1)
			for (i = 0; i < lr->inst->m; i++)
				if (lr->ks[i].x[j]) {
					lr->feasible_z += lr->inst->cost[i][j];
					lr->slack[i] -= lr->inst->size[i][j];
				}

	if (!lr->feasible) {
		lr->feasible = 1;
		for (i = 0; i < lr->inst->m; i++)
			if (lr->slack[i] < 0)
				lr->feasible = 0;

		for (j = 0; j < lr->inst->n && lr->feasible; j++) 
			if (lr->sum_x[j] != 1) {
				best_i = -1;
				best_r = INFINITY;
				for (i = 0; i < lr->inst->m; i++)
					if (lr->inst->size[i][j] < lr->slack[i] && (lr->sum_x[j] == 0 || lr->ks[i].x[j])) {
						r = lr->ks[i].p[j] / ((double) lr->inst->size[i][j]);
						if (r < best_r) {
							best_i = i;
							best_r = r;
						}
					}
				if (best_i == -1)
					lr->feasible = 0;
				else {
					lr->feasible_z += lr->inst->cost[best_i][j];
					lr->slack[best_i] -= lr->inst->size[best_i][j];
				}
			}
	}

	if (lr->feasible && (double) (lr->feasible_z - 1) <= *(lr->z_ub_ptr)) {
		*(lr->z_ub_ptr) = (double) (lr->feasible_z - 1);
		z_ub_update = 1;
	}

	return z_ub_update;
}


int lagrangian_dual_fix_variables (lagrelax_t* lr, lagrelax_state_t* lrs)
{
	int i, j, k = -1, n, w, rval = 0;
	int nnz, cap;
	knapsack_t* ks;
	uint8_t* x0;
	uint8_t* x1;
	double* rc;
	double base_z, z;
	float f;

	n = lr->inst->n - 1;
	for (i = 0; i < lr->inst->m; i++) {
		ks = &lr->ks[i];
		x0 = lrs->x0[i];
		x1 = lrs->x1[i];
		rc = lr->rc[i];

		base_z = 0.0;
		cap = ks->cap;
		for (j = 0; j < lr->inst->n; j++) {
			rc[j] = INFINITY;
			if (x1[j]) {
				cap -= ks->w[j];
				base_z += ks->p[j];
			}
		}
		nnz = 0;
		for (j = 0; j < lr->inst->n; j++) {
			if (!x0[j] && !x1[j]) {
				if (ks->w[j] > cap) {
					rc[j] = -INFINITY;
				}
				else 
				lr->index_list[nnz++] = j;
			}
		}
		
		if (nnz > 0) {

			j = lr->index_list[0];
			for (w = 0; w < ks->w[j]; w++)
				lr->f_val_cache[0][w] = 0.0f;
			f = (ks->p[j] < 0.0f) ? ks->p[j] : 0.0f;
			for (; w <= cap; w++) 
				lr->f_val_cache[0][w] = f;
			for (k = 1; k < nnz; k++) {
				j = lr->index_list[k];
				memcpy(lr->f_val_cache[k], lr->f_val_cache[k-1], (cap + 1) * sizeof(float));
				if (ks->p[j] < 0.0f) {
					for (w = ks->w[j]; w <= cap; w++) {
						f = lr->f_val_cache[k - 1][w - ks->w[j]] + ks->p[j];
						if (f < lr->f_val_cache[k][w]) 
							lr->f_val_cache[k][w] = f;
					}
				}
			}
	
			j = lr->index_list[nnz-1];
			for (w = cap - ks->w[j] + 1; w <= cap; w++) 
				lr->g_val_cache[nnz-1][w] = 0.0f;
			f = (ks->p[j] < 0.0f) ? ks->p[j] : 0.0f;
			for (w = cap - ks->w[j]; w >= 0; w--) 
				lr->g_val_cache[nnz-1][w] = f;
			for (k = nnz - 2; k >= 0; k--) {
				j = lr->index_list[k];
				memcpy(lr->g_val_cache[k], lr->g_val_cache[k+1], (cap + 1) * sizeof(float));
				if (ks->p[j] < 0.0f) {
					for (w = cap - ks->w[j]; w >= 0; w--) {
						f = lr->g_val_cache[k + 1][w + ks->w[j]] + ks->p[j];
						if (f < lr->g_val_cache[k][w]) 
							lr->g_val_cache[k][w] = f;
					}
				}
			}
	
			// compute rc
			for (k = 0; k < nnz; k++) {
				j = lr->index_list[k];
				if (ks->x[j]) {
					if (k == 0) 
						rc[j] = ks->opt - base_z - lr->g_val_cache[1][0];
					else if (k == nnz - 1) {
						rc[j] = ks->opt - base_z - lr->f_val_cache[nnz-2][cap];
					}
					else
						for (w = 0; w <= cap; w++) {
							f = ks->opt - base_z - lr->f_val_cache[k - 1][w] - lr->g_val_cache[k + 1][w];
							if (rc[j] == INFINITY || f > rc[j])
								rc[j] = f;
						}
				}
				else {
					rc[j] = -INFINITY;
					if (k == 0) 
						rc[j] = ks->opt - base_z - lr->g_val_cache[1][ks->w[j]] - ks->p[j];
					else if (k == nnz - 1) {
						if (cap >= ks->w[j]) {
							rc[j] = ks->opt - base_z - lr->f_val_cache[nnz-2][cap - ks->w[j]] - ks->p[j];
						}
					}
					else 
						for (w = ks->w[j]; w <= cap; w++) {
							f = ks->opt - base_z - lr->f_val_cache[k - 1][w - ks->w[j]] - ks->p[j] - lr->g_val_cache[k + 1][w];
							if (f > rc[j])
								rc[j] = f;
						}
				}
			}
		}
	}


	for (i = 0; i < lrs->inst->m; i++) 
		for (j = 0; j < lrs->inst->n; j++) 
			if (!lrs->x0[i][j] && !lrs->x1[i][j]) {
				z = 0.0;
				if (lr->ks[i].x[j]) {
					if (lr->sum_x[j] == 1) {
						z = -INFINITY;
						for (k = 0; k < lrs->inst->m; k++)
							if (k != i && !lrs->x0[k][j] && lr->rc[k][j] > z)
								z = lr->rc[k][j];
					}
					if (lr->z - lr->rc[i][j] - z > *(lr->z_ub_ptr)) {
						rval = lrs->x1[i][j] = 1;
						for (k = 0; k < lrs->inst->m; k++) 
							if (k != i)
								lrs->x0[k][j] = 1;
					}
				}
				else {
					for (k = 0; k < lrs->inst->m; k++) 
						if (k != i && lr->ks[i].x[j])
							z += lr->rc[i][j];
					if (lr->z - lr->rc[i][j] - z > *(lr->z_ub_ptr))
						rval = lrs->x0[i][j] = 1;
				}
			}

	for (j = 0; j < lrs->inst->n; j++) 
		if (lr->sum_x[j] <= 1) {
			for (n = i = 0; i < lrs->inst->m; i++) {
				if (lrs->x0[i][j]) {
					assert(!lrs->x1[i][j]);
					++n;
				}
				else
					k = i;
			}
			if (n == lrs->inst->m) 
				return rval;
			if (lrs->x1[k][j] == 0 && n + 1 == lrs->inst->m) {
				assert(!lrs->x0[k][j]);
				rval = lrs->x1[k][j] = 1;
			}
		}

	return rval;
}


