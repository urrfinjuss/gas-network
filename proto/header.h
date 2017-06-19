#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define	DEBUG_MODE	1
#define	FTYPE	double

typedef struct gas_pipe				gpipe, *gpipe_ptr;
typedef struct gas_node				gnode, *gnode_ptr;
typedef struct gas_compressor	gcomp, *gcomp_ptr;

typedef struct parameters {
	char					name[80];
	char					intmethod[80];
	unsigned int	skip;
	unsigned int	nsteps;
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
	FTYPE	len;
	FTYPE	wid;
	gnode_ptr	left, right;
};

struct gas_compressor {
	unsigned int	idx;
	gnode_ptr			loc;
	gpipe_ptr			dest;
	FTYPE					*boost;
};

// message.c
extern void dmesg(char *line, unsigned int flag);

// input.c
extern void read_input_file(char *fname);
extern void read_network_list(char *fname);
extern void read_nodes_list(char *fname);
extern void read_pipes_list(char *fname);
extern void read_comps_list(char *fname);
extern void update_nodes_comps();
extern void update_nodes_pipes();

extern int read_pipe(char *fname, gpipe_ptr pipe);
extern int install_compressors(char *fname);
extern int load_initial_data();
extern FILE* fopen_set(char *fname, char* flag, const unsigned skip_lines);
extern long int count_data(char *fname);

// network.c
extern int call_init_network(char *filename);
extern int build_network();

// temporal.c
extern int call_init_temporal();

// verify.c
extern void verify_consistency();
extern void verify_input_conf();
extern void verify_network_conf();
extern void verify_node_conf();
extern void verify_pipe_conf();
extern void verify_comp_conf();

// output.c
extern void network_snapshot();
extern void dump_pipe(gpipe_ptr in, char *line);

// nodes.c
extern void call_init_nodes();

// comps.c
extern void call_init_comps();

// List global variables
extern network	net;
extern params		par;
