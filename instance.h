#ifndef INSTANCE_H
#define INSTANCE_H

#include <stdio.h>
#include <math.h>

#ifndef INFINITY
	#define INFINITY DBL_MAX
#endif

typedef struct
{
	char name[128];
	int m, n;
	int* capacity;
	int** cost;
	int** size;

	double* init_mul;
	double init_z_lb;
	int optimal_value;

	int max_subgradient_iterations;
	int max_search_seconds;
	double search_progress_output_period;
	FILE* progress_output;
	FILE* log;

} gap_inst_t;


gap_inst_t gap_instance_load (int argc, char** argv);
void gap_instance_free (gap_inst_t* inst);


#endif
