#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct network network;
typedef struct node node, *node_ptr;
typedef struct noise_p noise_p;
typedef struct pipe {
  double L;
  int N;
  int BCflag;
  node *left, *right;
  double *flx, *prss;
} pipe, *pipe_ptr;

struct node {
  pipe_ptr left, right;
  double prss;
  double compress;
};

struct noise_p {
  double tc;	// correlation time 
  double lc;    // correlation distance
  double A;	// noise amplitude
};

struct network {
  char fname[80], nname[80];
  int nnodes;
  int nlinks;
  noise_p *noise;
  node_ptr node;
  pipe_ptr link; 
};

extern int err_msg(char* msg);
extern int noise_status(noise_p *noise);
