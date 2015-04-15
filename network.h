#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct node node, *node_ptr;
typedef struct noise params;
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

struct noise {
  double tc;
  double A; // noise ampl
};
