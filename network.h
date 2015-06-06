#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define DEBUGMODE 1

typedef struct network network;
typedef struct node node, *node_ptr;
typedef struct noise_p noise_p;
typedef struct pipe {
  double L;
  int N;
  int BCflag;
  node_ptr *adjacent;
  double *flx, *prss;
} pipe, *pipe_ptr;

struct node {
  int adj_n;
  pipe_ptr *adj_p;
  node_ptr *adj_k;
  double prss;
  double compress;
};

struct noise_p {
  double tc;	// correlation time 
  double lc;    // correlation distance
  double A;	// noise amplitude
};

struct network {
  char fname[80], nname[80], matname[80];
  int nnodes;
  int nlinks;
  noise_p *noise;
  node_ptr *knot;
  pipe_ptr *link; 
};

// aux.c
extern void err_msg(char* msg);
extern void debug_msg(char* msg);

//
extern int noise_status(noise_p *noise);

// init.c
extern void init_links(network *net);
extern void allocate_memory(network *net, int *lm);

// gsl_noise.c
extern int init_network(char *filename, network *net);
