#include "instance.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>



double var_weight (int cost, int size)
{
	return ((double) abs(cost)) / ((double) size);
}


gap_inst_t gap_instance_load (int argc, char** argv) 
{
	gap_inst_t inst;
	int i, j;
	FILE* f;
	float flo;
	char buf[256];

	for (i = strlen(argv[1]); i > 0; i--) 
		if (argv[1][i - 1] == '/')
			break;
	strcpy(inst.name, &argv[1][i]);

	inst.log = stdout;
	inst.progress_output = stderr;

	if ((argc != 6)
		|| (NULL == (f = fopen(argv[1], "r")))
		|| (1 != sscanf(argv[2], "%i", &inst.optimal_value))
		|| (1 != sscanf(argv[3], "%i", &inst.max_subgradient_iterations))
		|| (1 != sscanf(argv[4], "%i", &inst.max_search_seconds))
		|| (1 != sscanf(argv[5], "%f", &flo)))
		goto load_err;
	inst.search_progress_output_period = flo;

	if (2 != fscanf(f, "%i %i\n", &inst.m, &inst.n)) goto load_err;

	inst.capacity = (int*) calloc(inst.m, sizeof(int));
	inst.cost = (int**) calloc(inst.m, sizeof(int*));
	inst.size = (int**) calloc(inst.m, sizeof(int*));
	if (NULL == inst.capacity || NULL == inst.cost || NULL == inst.size) goto load_err;
		
	for (i = 0; i < inst.m; i++) {
		inst.cost[i] = (int*) calloc(inst.n, sizeof(int));
		inst.size[i] = (int*) calloc(inst.n, sizeof(int));
		if (NULL == inst.cost[i] || NULL == inst.size[i]) goto load_err;
	}

	for (i = 0; i < inst.m; i++) 
		for (j = 0; j < inst.n; j++) 
			if (1 != fscanf(f, "%i", &inst.cost[i][j])) goto load_err;

	for (i = 0; i < inst.m; i++) 
		for (j = 0; j < inst.n; j++) 
			if (1 != fscanf(f, "%i", &inst.size[i][j])) goto load_err;

	for (i = 0; i < inst.m; i++) 
		if (1 != fscanf(f, "%i", &inst.capacity[i])) goto load_err;
	fclose(f);

	inst.init_mul = (double*) calloc(inst.n, sizeof(double));
	inst.init_z_lb = -INFINITY;

	strcpy(buf, argv[1]);
	strcat(buf, ".mul");
	f = fopen(buf, "r");
	if (f != NULL) {
		fscanf(f, "%f\n", &flo);
		fscanf(f, "%f\n", &flo);
		inst.init_z_lb = flo;
		for (j = 0; j < inst.n; j++) {
			fscanf(f, "%f\n", &flo);
			inst.init_mul[j] = flo;
		}
	}
	fclose(f);


	return inst;

load_err:
	fprintf(stderr, "\ninstance file \t%s\n", argv[1]);
	fprintf(stderr, "ERROR\t%i\t%i\t%i\t%f\n", inst.optimal_value, inst.max_subgradient_iterations, inst.max_search_seconds, inst.search_progress_output_period);
	inst.n = -1;
	return inst;
}



void gap_instance_free (gap_inst_t* inst)
{
	int i;

	for (i = 0; i < inst->m; i++) {
		free(inst->cost[i]);
		free(inst->size[i]);
	}
	free(inst->cost);
	free(inst->init_mul);
	free(inst->size);
	free(inst->capacity);
}


