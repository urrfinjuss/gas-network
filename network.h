#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define DEBUGMODE 1

typedef struct network network;
typedef struct node node, *node_ptr;
typedef struct noise_p noise_p;
typedef struct gpipe{
  double L;
  double d;
  double W_l, Wb_r; // characteristics required to update pipe
  double Wb_l, W_r; // characteristics required to update node
  int N;
  int BCflag;
  node_ptr left;
  node_ptr right;
  double *f, *p, *q;
} gpipe, *gpipe_ptr;

struct node {
  int adj_n;
  int nl, nr;
  gpipe_ptr *outg;
  gpipe_ptr *incm;
  gpipe_ptr *adj_p;
  node_ptr *adj_k;
  node_ptr *left;
  node_ptr *right;
  double P, F;
  double cratio;
};

struct noise_p {
  double tc;	// correlation time 
  double lc;    // correlation distance
  double A;	// noise amplitude
};

struct network {
  char fname[80], nname[80], dname[80], current_dir[1024];
  char incname[256];
  int nnodes;
  int nlinks;
  int npcent;
  double c, diss, tmax;
  noise_p *noise;
  node_ptr *knot;
  gpipe_ptr *link; 
};

// aux.c
extern void err_msg(char* msg);
extern void debug_msg(char* msg);

//
extern int noise_status(noise_p *noise);

// init.c
extern void init_data(network *net);
extern void init_links(network *net);
extern void init_arrays(network *net);
extern void allocate_memory(network *net, double *lm);
extern void rescale_data(network *net);
extern int load_data(FILE *fh, gpipe_ptr lnk);
extern void save_data(FILE* fh, gpipe_ptr lnk, network* net);

// gsl_noise.c
extern int init_network(char *filename, network *net);

// drawing.c
extern void mgl_init_draw(network *net);
extern void mgl_draw_pipe(gpipe_ptr in, network *net, char* fname, char* title);
extern void mgl_draw_pressure(network *net, char* fname, char* title);
extern void mgl_draw_flux(network *net, char* fname, char* title);
extern void mgl_draw_network(network *net, char* fname, char* title);

//evolve.c
extern void init_nodes(network *net);
extern void init_evolve(network *net);
extern void hyperbolic_step(network *net, double time);
extern void nonlinear_hstep(network* net, double dt);
extern void evolve_network(network *net);
extern void split_step2(network *net, double dt, double time);




