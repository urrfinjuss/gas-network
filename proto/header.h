#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define	DEBUG_MODE	1
#define	FTYPE	double
#define	pi	acosl(-1.0L)

typedef struct gas_pipe				gpipe, *gpipe_ptr;
typedef struct gas_node				gnode, *gnode_ptr;
typedef struct gas_compressor	gcomp, *gcomp_ptr;

typedef struct parameters {
	char					name[80];
	char					intmethod[80];
	const gsl_interp_type *iType;
	unsigned int	skip;
	unsigned int	nsteps;
	unsigned int	curr;
	FTYPE					dt;
	FTYPE					dx;
	FTYPE					tmax;
	FTYPE					sound;
	FTYPE					friction;
} params;

typedef struct gas_network {
	unsigned int	nodes;
	unsigned int	pipes;
	unsigned int	comps;
	unsigned int	consistent;
	gcomp_ptr	comp;
	gpipe_ptr	pipe;
	gnode_ptr	node;
} network;

struct gas_node {
	FTYPE	*var;
	FTYPE p[2];
	unsigned int	type;
	unsigned int	nleft;
	unsigned int	nright;
	unsigned int	ncomp;
	gcomp_ptr	*comp;
	gpipe_ptr	*left;
	gpipe_ptr	*right;
};

struct gas_pipe {
	long int N;
	FTYPE	dx;
	FTYPE	*y[2];
	FTYPE	fs[2];
	FTYPE	fd[2];
	FTYPE	len;
	FTYPE	wid;
	gcomp_ptr	cmp;
	gnode_ptr	left, right;
};

struct gas_compressor {
	unsigned int	idx;
	gnode_ptr			loc;
	gpipe_ptr			dest;
	FTYPE					*boost;
};

// main.c

// input.c
extern FILE* fopen_set(char *fname, char* flag, const unsigned skip_lines);
extern long int count_data(char *fname);
extern void read_input_file(char* fname);
extern void read_network_list(char *fname);
extern void read_nodes_list(char *fname);
extern void read_pipes_list(char *fname);
extern void read_comps_list(char *fname);
extern void update_nodes_comps();
extern void update_nodes_pipes();

// output.c
extern void network_snapshot();
extern void dump_pipe(gpipe_ptr in, char *line);

// message.c
extern void dmesg(char *line, unsigned int flag);

// init.c
extern void call_init_network(char *filename);
extern void call_init_nodes();
extern void call_init_comps();
extern void call_init_pipes();

// network.c

// pipes.c
extern void init_forward_interior();
extern void free_forward_interior();
extern void forward_interior();
extern void forward_pipe(int j);

// nodes.c
extern void sync_nodes();
extern void forward_nodes();
extern void forward_pressure_node(gnode_ptr in);
extern void forward_transport_node(gnode_ptr in);
extern void calculate_pressure(gnode_ptr in);
extern void forward_incoming_pipes(gnode_ptr in);
extern void forward_outgoing_pipes(gnode_ptr in);

// comps.c

// verify.c
extern void verify_input_conf();
extern void verify_network_conf();
extern void verify_node_conf();
extern void verify_pipe_conf();
extern void verify_comp_conf();
extern void verify_interp_nodes(long int j, long int nbase, FTYPE *Tbase, FTYPE *Vbase);
extern void verify_interp_comps(long int j, long int nbase, FTYPE *Tbase, FTYPE *Vbase);
extern void verify_set_pipes(long int j);

// List global variables
extern network	net;
extern params		par;
