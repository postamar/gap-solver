#include "bb.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <sys/resource.h>

const int alloc_step = 10;



void bb_branch (bb_t* bb, int parent_node_index)
{
	int i, j, k, branch_j;
	int new_node_index;
	bb_node_t* node;
	bb_node_t* new_node;

	node = &bb->queue[parent_node_index];
	for (branch_j = -1, j = 0; j < bb->inst->n; j++) {
		if (bb->lagr.sum_x[j] != 1 && (branch_j == -1 || node->lrs.mul[j] > node->lrs.mul[branch_j])) 
			branch_j = j;
	}
	assert(branch_j != -1);

	for (i = 0; i < bb->inst->m; i++) {
		if (!bb->queue[parent_node_index].lrs.x0[i][branch_j]) {
			new_node_index = bb_node_create_copy(bb, parent_node_index);
			new_node = &bb->queue[new_node_index];
			new_node->lrs.x1[i][branch_j] = 1;
			for (k = 0; k < bb->inst->m; k++)
				if (k != i) 
					new_node->lrs.x0[k][branch_j] = 1;
			if (!bb_node_pre_eval(bb, new_node_index)) 
				bb->queue[new_node_index].id = -1;
		}
	}
}


void bb_node_move_to_garbage (bb_t* bb, int node_index) 
{
	if (bb->garbage_queue_size == bb->garbage_queue_current_alloc) {
		bb->garbage_queue_current_alloc += alloc_step;
		bb->garbage_queue = (bb_node_t*) realloc(bb->garbage_queue, bb->garbage_queue_current_alloc * sizeof(bb_node_t));
	}

	bb->garbage_queue[bb->garbage_queue_size++] = bb->queue[node_index];
	bb->queue[node_index] = bb->queue[--bb->queue_size];
}



int bb_search (bb_t* bb, int best_first_search)
{
	int k, n, node_index, old_queue_size, time_flag = 0;
	struct rusage resource_usage;
	struct rusage prev_rusage;
	double init_z_lb = bb->queue[0].lrs.z_lb;

	getrusage(RUSAGE_SELF, &prev_rusage);
	fprintf(bb->inst->progress_output, "\niter\teval\tgap\t\tlb\n----\t----\t---\t\t--\n");

	bb_node_pre_eval(bb, 0);
	for (n = 0; bb->queue_size > 0 && !time_flag && !bb->found_feasible; n++)
	{
		if (bb->garbage_queue_size > 1000) {
			for (k = 0; k < bb->garbage_queue_size; k++) {
				lagrelax_state_destroy(&bb->garbage_queue[k].lrs);
			}
			bb->garbage_queue_size = 0;
		}

		old_queue_size = bb->queue_size;
		bb->queue_size = 0;
		for (node_index = -1, bb->z_lb = bb->z_ub, k = 0; k < old_queue_size; k++) 
			if (bb->queue[k].id == -1 || bb->queue[k].lrs.z_lb > bb->z_ub || bb->queue[k].lrs.z_lb == INFINITY) {
				if (bb->garbage_queue_size == bb->garbage_queue_current_alloc) {
					bb->garbage_queue_current_alloc += alloc_step;
					bb->garbage_queue = (bb_node_t*) realloc(bb->garbage_queue, bb->garbage_queue_current_alloc * sizeof(bb_node_t));
				}
				bb->garbage_queue[bb->garbage_queue_size++] = bb->queue[k];
			}
			else {
				if (bb->queue[k].lrs.z_lb <= bb->z_lb) {
					bb->z_lb = bb->queue[k].lrs.z_lb;
					node_index = bb->queue_size;
				}
				bb->queue[bb->queue_size++] = bb->queue[k];
			}

		if (node_index == -1) 
			break;
		else if (!best_first_search)
			node_index = bb->queue_size - 1;

		getrusage(RUSAGE_SELF, &resource_usage);
		if (bb->inst->max_search_seconds >= 0 && resource_usage.ru_utime.tv_sec >= bb->inst->max_search_seconds) 
			time_flag = 1;

		if (time_flag || bb->inst->search_progress_output_period + ((double) prev_rusage.ru_utime.tv_sec) + ((double) prev_rusage.ru_utime.tv_usec) / 1000000.0
							<= ((double) resource_usage.ru_utime.tv_sec) + ((double) resource_usage.ru_utime.tv_usec) / 1000000.0) {
			if (n > 0) 
				fprintf(bb->inst->progress_output, "%i\t%i\t%f%%\t%f\n", n, bb->lagr.n_eval, 100.0 * (1.0 - (bb->z_lb - init_z_lb) / (bb->z_ub - init_z_lb)), bb->z_lb);
			prev_rusage = resource_usage;
		}

		if (!time_flag) {
			if (bb_node_eval(bb, node_index)) {
				bb_branch(bb, node_index);
			}
			bb->queue[node_index].id = -1;
		}
	}


	if (!time_flag)
		fprintf(bb->inst->progress_output, "%i\t%i\t0.000000%%\tdone\n\n", n, bb->lagr.n_eval);


	return (bb->found_feasible) ? 1 : ((time_flag) ? -1 : 0);
}



int bb_node_create_copy (bb_t* bb, int src_node_index)
{
	bb_node_t* src;
	bb_node_t* dst;

	if (bb->queue_size == bb->queue_current_alloc) {
		bb->queue_current_alloc += alloc_step;
		bb->queue = (bb_node_t*) realloc(bb->queue, bb->queue_current_alloc * sizeof(bb_node_t));
	}

	if (bb->garbage_queue_size > 0)		
		bb->queue[bb->queue_size++] = bb->garbage_queue[--bb->garbage_queue_size];
	else 
		bb->queue[bb->queue_size++].lrs = lagrelax_state_create(bb->inst);

	src = &bb->queue[src_node_index];
	dst = &bb->queue[bb->queue_size - 1];
	lagrelax_state_copy(&dst->lrs, &src->lrs);
	dst->id = bb->n_nodes++;
	dst->depth = src->depth + 1;

	return bb->queue_size - 1;
}


int bb_node_pre_eval (bb_t* bb, int node_index)
{
	assert(node_index >= 0);
	assert(node_index < bb->queue_size);

	bb_node_t* node = &bb->queue[node_index];

	if (lagrelax_eval(&bb->lagr, node->lrs.x0, node->lrs.x1, node->lrs.mul, NULL)) {
		bb->found_feasible = 1;
		return 0;
	}
	node->lrs.z_lb = bb->lagr.z;

	return (node->lrs.z_lb != INFINITY && node->lrs.z_lb <= bb->z_ub);
}

int bb_node_eval (bb_t* bb, int node_index)
{
	assert(node_index >= 0);
	assert(node_index < bb->queue_size);

	int i, j, n0, n1;
	bb_node_t* node = &bb->queue[node_index];

	assert(node->lrs.z_lb != INFINITY && node->lrs.z_lb <= bb->z_ub);

	++bb->n_node_eval;

	if (lagrangian_dual_solve(&bb->lagr, &node->lrs, bb->inst->max_subgradient_iterations)) {
		bb->found_feasible = 1;
		return 0;
	}


	if (node->lrs.z_lb == INFINITY || node->lrs.z_lb > bb->z_ub) 
		return 0;

	if (lagrangian_dual_fix_variables(&bb->lagr, &node->lrs)) {
		for (j = 0; j < bb->inst->n; j++) {
			n0 = n1 = 0;
			for (i = 0; i < bb->inst->m; i++) 
				if (node->lrs.x0[i][j]) {
					assert(!node->lrs.x1[i][j]);
					++n0;
				}
				else if (node->lrs.x1[i][j]) {
					assert(n1 == 0);
					assert(!node->lrs.x0[i][j]);
					++n1;
				}
			if (n0 == bb->inst->m) {
				return 0;
			}
		}
	}

	return 1;
}


bb_t bb_create (gap_inst_t* inst, int init_z_ub, lagrelax_state_t* root_lrs)
{
	bb_t bb;
	bb_node_t* node;

	bb.inst = inst;

	bb.garbage_queue_size = 0;
	bb.garbage_queue_current_alloc = inst->n * inst->m;
	bb.garbage_queue = (bb_node_t*) calloc(bb.garbage_queue_current_alloc, sizeof(bb_node_t));

	bb.queue_size = 1;
	bb.queue_current_alloc = inst->n * inst->m;
	bb.queue = (bb_node_t*) calloc(bb.queue_current_alloc, sizeof(bb_node_t));

	node = &bb.queue[0];
	node->id = 0;
	node->depth = 0;
	node->lrs = lagrelax_state_create(inst);
	lagrelax_state_copy(&node->lrs, root_lrs);

	bb.n_nodes = 1;
	bb.n_node_eval = 0;
	bb.old_z_ub = bb.z_ub = (double) init_z_ub;
	bb.found_feasible = 0;

	return bb;
}



void bb_free (bb_t* bb)
{
	int k;

	for (k = 0; k < bb->garbage_queue_size; k++) {
		lagrelax_state_destroy(&bb->garbage_queue[k].lrs);
	}
	free(bb->garbage_queue);

	for (k = 0; k < bb->queue_size; k++) {
		lagrelax_state_destroy(&bb->queue[k].lrs);
	}
	free(bb->queue);

	lagrelax_destroy(&bb->lagr);
}

