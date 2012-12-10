#include "instance.h"
#include "bb.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/resource.h>


int main (int argc, char** argv)
{
	bb_t bb;
	gap_inst_t inst = gap_instance_load(argc, argv);
	lagrelax_t root_lagrelax;
	lagrelax_state_t root_lagrelax_state;
	double global_z_ub = INFINITY;
	int lower_z_ub, z_ub, rval, done_flag = 0;
	int n_lag_eval;
	int n_node_eval = 1;
	struct rusage resource_usage;

	if (inst.n == -1)
		return -1;

	fprintf(inst.log, "\ninstance\t%s\n", inst.name);
	fflush(inst.log);

	root_lagrelax = lagrelax_create(&inst, &global_z_ub);
	root_lagrelax_state = lagrelax_state_create(&inst);

	memcpy(root_lagrelax_state.mul, inst.init_mul, inst.n * sizeof(double));
	root_lagrelax_state.z_lb = inst.init_z_lb;
	lagrangian_dual_solve(&root_lagrelax, &root_lagrelax_state, 10000);
	n_lag_eval = root_lagrelax.n_eval;
	lower_z_ub = (inst.optimal_value == -1) ? (int) round(ceil(root_lagrelax_state.z_lb)) : inst.optimal_value - 1;

	fprintf(inst.log, "%i", lower_z_ub - 1);
	if (global_z_ub != INFINITY && round(global_z_ub) == lower_z_ub) {
		done_flag = 1;
		fprintf(inst.log, "*");
	}
	getrusage(RUSAGE_SELF, &resource_usage);
	fprintf(inst.log, "\t%.2f", (double) resource_usage.ru_utime.tv_sec + (double)resource_usage.ru_utime.tv_usec / 1000000.0);
	fprintf(inst.log, "\t%i\t%i\n", n_lag_eval, n_node_eval);
	fflush(inst.log);

	for (z_ub = lower_z_ub; !done_flag && (inst.optimal_value == -1 || z_ub < inst.optimal_value); z_ub++) {
		bb = bb_create(&inst, z_ub, &root_lagrelax_state);
		bb.lagr = lagrelax_create(bb.inst, &bb.z_ub);
		rval = bb_search(&bb, 0);

		fprintf(inst.log, "%i", z_ub);
		if (rval == 1) {
			fprintf(inst.log, "*");
			done_flag = 1;
		}
		else if (rval == -1) {
			fprintf(inst.log, "t");
			done_flag = 1;
		}
		getrusage(RUSAGE_SELF, &resource_usage);
		fprintf(inst.log, "\t%.2f", (double) resource_usage.ru_utime.tv_sec + (double)resource_usage.ru_utime.tv_usec / 1000000.0);
		n_lag_eval += bb.lagr.n_eval;
		n_node_eval += bb.n_node_eval;
		fprintf(inst.log, "\t%i\t%i\n", n_lag_eval, n_node_eval);
		fflush(inst.log);

		bb_free(&bb);
	}

	lagrelax_state_destroy(&root_lagrelax_state);
	lagrelax_destroy(&root_lagrelax);
	gap_instance_free(&inst);

	return 0;
}


