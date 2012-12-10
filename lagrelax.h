#ifndef LAGRELAX_H
#define LAGRELAX_H

#include "knapsack.h"
#include "bundle.h"

#include <stdint.h>


typedef struct
{
	gap_inst_t* inst;
	uint8_t** x0;
	uint8_t** x1;
	double z_lb;
	double* mul;
} lagrelax_state_t;


typedef struct
{
	gap_inst_t* inst;
	double* z_ub_ptr;
	int n_eval;

	knapsack_t* ks;

	double z;
	int feasible;
	int feasible_z;
	int* sum_x;
	double** rc;

	bundle_t *bundle;
	lagrelax_state_t* bundle_lrs;

	int* index_list;
	float** f_val_cache;
	float** g_val_cache;

	int* slack;
	minknap_t* minknap_data;

} lagrelax_t;



lagrelax_t lagrelax_create (gap_inst_t* inst, double* z_ub_ptr);
void lagrelax_destroy (lagrelax_t* lr);
int lagrelax_eval (lagrelax_t* lr, uint8_t** x0, uint8_t** x1, double *mul, double *subg);

lagrelax_state_t lagrelax_state_create (gap_inst_t* inst);
void lagrelax_state_copy (lagrelax_state_t* dst, const lagrelax_state_t* src);
void lagrelax_state_destroy (lagrelax_state_t* lrs);
int lagrelax_state_check (lagrelax_state_t* lrs); 
void lagrelax_state_print_fixed (lagrelax_state_t* lrs); 

int lagrangian_dual_fix_variables (lagrelax_t* lr, lagrelax_state_t* lrs);
int lagrangian_dual_solve (lagrelax_t* lr, lagrelax_state_t* lrs, int bundle_iterations);


#endif

