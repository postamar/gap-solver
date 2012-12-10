#include "knapsack.h"


#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>



knapsack_t knapsack_create (gap_inst_t* inst, int i, minknap_t* minknap_data)
{
	int j;

	knapsack_t ks;
	ks.i = i;
	ks.n = inst->n;
	ks.slack = ks.cap = inst->capacity[i];
	ks.p = (double*) calloc(inst->n, sizeof(double));
	ks.w = (int*) calloc(inst->n, sizeof(int));
	ks.x = (int*) calloc(inst->n, sizeof(int));
	for (j = 0; j < inst->n; j++) 
		ks.w[j] = inst->size[i][j];
	ks.item_list = (knapsack_enum_item_t*) calloc(inst->n, sizeof(knapsack_enum_item_t));
	ks.enum_x = (int*) calloc(inst->n, sizeof(int));
	knapsack_reset(&ks, inst);
	ks.minknap_n = ks.minknap_z = ks.minknap_c = 0;
	ks.minknap_p = (int*) calloc(inst->n, sizeof(int));
	ks.minknap_w = (int*) calloc(inst->n, sizeof(int));
	ks.minknap_x = (int*) calloc(inst->n, sizeof(int));
	ks.minknap_j = (int*) calloc(inst->n, sizeof(int));
	ks.minknap_data = minknap_data;
	return ks;
}


void knapsack_destroy(knapsack_t* ks)
{
	free(ks->p);
	free(ks->w);
	free(ks->x);
	free(ks->item_list);
	free(ks->enum_x);
	free(ks->minknap_p);
	free(ks->minknap_w);
	free(ks->minknap_x);
	free(ks->minknap_j);
}


void knapsack_reset (knapsack_t* ks, gap_inst_t* inst)
{
	int j;

	ks->item_list_size = 0;
	ks->opt = INFINITY;
	ks->slack = ks->cap;
	for (j = 0; j < ks->n; j++) {
		ks->p[j] = (double) inst->cost[ks->i][j];
	}
}






void knapsack_solve (knapsack_t* ks, uint8_t* x0, uint8_t* x1)
{
	int j, k, wsum = 0;
	double z = 0.0;
	double minp = 0.0;

	for (j = 0; j < ks->n; j++) {
		ks->x[j] = 0;
		if (x1[j]) {
			ks->x[j] = 1;
			ks->slack -= ks->w[j];
			z += ks->p[j];
		}
		if (ks->p[j] < 0.0 && !x0[j]) {
			wsum += ks->w[j];
			if (ks->p[j] < minp)
				minp = ks->p[j];
		}
	}
	if (ks->slack < 0)
		return;

	ks->minknap_n = 0;
	for (j = 0; j < ks->n; j++) 
		if (ks->p[j] < 0.0 && !x0[j] && !x1[j]) {
			ks->minknap_p[ks->minknap_n] = (int) lround(ks->p[j] * 1e5 / minp);
			assert(ks->minknap_p[ks->minknap_n] >= 0);
			ks->minknap_w[ks->minknap_n] = ks->w[j];
			ks->minknap_j[ks->minknap_n] = j;
			++ks->minknap_n;
		}


	ks->minknap_c = ks->slack;
	ks->opt = z;

	if (wsum <= ks->slack) {
		ks->slack -= wsum;
		for (k = 0; k < ks->minknap_n; k++) {
			j = ks->minknap_j[k];
			ks->x[j] = 1;
			ks->opt += ks->p[j];
		}
	}
	else if (ks->minknap_n > 1) {
		ks->minknap_z = minknap(ks->minknap_data, ks->minknap_n, ks->minknap_p, ks->minknap_w, ks->minknap_x, ks->minknap_c);
		for (k = 0; k < ks->minknap_n; k++) 
			if (ks->minknap_x[k]) {
				j = ks->minknap_j[k];
				ks->x[j] = 1;
				ks->opt += ks->p[j];
				ks->slack -= ks->w[j];
			}
	}
}


