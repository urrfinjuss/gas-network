#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define DEBUGMODE 0

typedef struct network network;
typedef struct node node, *node_ptr;
typedef struct noise_p noise_p;
typedef struct gpipe{
  double L;
  double d;
  int N;
  int BCflag;
  node_ptr left;
  node_ptr right;
  double *f, *p, *q;
} gpipe, *gpipe_ptr;

struct node {
  int adj_n;
  gpipe_ptr *adj_p;
  node_ptr *adj_k;
  double P, F;
  double cratio;
};

struct noise_p {
  double tc;	// correlation time 
  double lc;    // correlation distance
  double A;	// noise amplitude
};

struct network {
  char fname[80], nname[80], matname[80], dname[80], current_dir[1024];
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
extern void init_links(network *net);
extern void init_arrays(network *net);
extern void allocate_memory(network *net, int *lm);

// gsl_noise.c
extern int init_network(char *filename, network *net);

// ic.c
extern void init_data(network *net);
extern void load_data(FILE *fh, gpipe_ptr lnk);
extern void rescale_data(network *net);
// drawing.c
extern void mgl_init_draw(network *net);
extern void mgl_draw_pipe(gpipe_ptr in, network *net, char* fname, char* title);
extern void mgl_draw_pressure(network *net, char* fname, char* title);
extern void mgl_draw_flux(network *net, char* fname, char* title);

//evolve.c
extern void init_evolve(network *net);
extern void hyperbolic_step(network *net);
extern void evolve_network(network *net);



