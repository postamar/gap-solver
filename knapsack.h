#ifndef KNAPSACK_H
#define KNAPSACK_H

#include "instance.h"
#include "minknap.h"

#include <stdint.h>

typedef struct
{
	double p;
	long w, j;
} knapsack_enum_item_t;


typedef struct
{
	double* p;
	int* w;
	int cap;
	int i, n;

	int* x;
	int slack;
	double opt;


	knapsack_enum_item_t* item_list;
	int* enum_x;
	int item_list_size;

	int minknap_n;
	int minknap_z;
	int minknap_c;
	int* minknap_p;
	int* minknap_w;
	int* minknap_x;
	int* minknap_j;

	minknap_t* minknap_data;
} knapsack_t;


knapsack_t knapsack_create (gap_inst_t* inst, int i, minknap_t* minknap_data);
void knapsack_destroy (knapsack_t* ks);
void knapsack_reset (knapsack_t* ks, gap_inst_t* inst);
void knapsack_solve (knapsack_t* ks, uint8_t* x0, uint8_t* x1);


#endif


