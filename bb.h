#ifndef BB_H
#define BB_H

#include "lagrelax.h"

#include <stdint.h>


extern const int alloc_step;

typedef struct
{
	int id;
	int depth;
	lagrelax_state_t lrs;
} bb_node_t;



typedef struct
{
	gap_inst_t* inst;

	bb_node_t* garbage_queue;
	int garbage_queue_size, garbage_queue_current_alloc;

	bb_node_t* queue;
	int queue_size, queue_current_alloc;

	lagrelax_t lagr;

	int n_nodes;
	int n_node_eval;
	double z_lb;
	double z_ub;
	double old_z_ub;
	int found_feasible;

} bb_t;


bb_t bb_create (gap_inst_t* inst, int init_z_ub, lagrelax_state_t* root_lrs);
int bb_search (bb_t* bb, int best_first_search);
void bb_free (bb_t* bb);

void bb_node_alloc (bb_t* bb);
int bb_node_create_copy (bb_t* bb, int src_node_index);
int bb_node_eval (bb_t* bb, int node_index);
int bb_node_pre_eval (bb_t* bb, int node_index);
int bb_node_fix_x (bb_t* bb, bb_node_t* node);
void bb_node_remove (bb_t* bb, int node_index);
void bb_node_move_to_garbage (bb_t* bb, int node_index);
void bb_node_free(bb_t* bb, bb_node_t* node);

#endif

