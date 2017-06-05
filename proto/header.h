#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define	DEBUG_MODE	1
#define	FTYPE	double

typedef struct gas_pipe				gpipe, *gpipe_ptr;
typedef struct gas_node				gnode, *gnode_ptr;
typedef struct gas_compressor	gcomp, *gcomp_ptr;

typedef struct parameters {
	char					name[80];
	unsigned int	skip;
	FTYPE					tmax;
	FTYPE					sound;
	FTYPE					friction;
} params;

typedef struct gas_network {
	unsigned int	nodes;
	unsigned int	pipes;
	unsigned int	compressors;
	unsigned int	consistent;
	gcomp_ptr	comp;
	gpipe_ptr	pipe;
	gnode_ptr	node;
} network;

struct gas_node {
	char node_type[80];
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
extern int read_input_file(char *fname);
extern int read_network_list(char *fname);
extern int read_pipe(char *fname, gpipe_ptr pipe);
extern int install_compressors(char *fname);
extern int load_initial_data();

// init.c
extern int call_init_network(char *filename);
extern int verify_consistency();
extern int build_network();

// output.c
extern void network_snapshot();
extern void dump_pipe(gpipe_ptr in, char *line);

// List global variables
extern network	net;
extern params		par;
