#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>


#define LINEAR_BC 1
#define DEBUGMODE 0
#define BALANCECHECK 0
#define ADIABATIC 0
#define ADIABATIC_SAVE 0
#define ZERO_GAMMA 1
#define NOISE 0
#define POOR_MANS_NOISE 1
#define PIPE_SAVE 1
#define pi 2*acos(0)
#define RFLAG_P_FIXED 0

typedef struct network network, *network_ptr;
typedef struct node node, *node_ptr;
typedef struct noise_p noise_p;
typedef struct compressor compressor, *compressor_ptr;
typedef struct params params, *params_ptr;
typedef struct gpipe{
  double L, Fl, Fr, Fm;
  double Pl, Pm, Pr;
  double d, *gamma;
  double W_l, Wb_r; // characteristics required to update pipe
  double Wb_l, W_r; // characteristics required to update node
  int N, key;
  int BCflag;
  node_ptr left;
  node_ptr right;
  fftw_complex *fx;
  double *f, *p, *q;
  compressor_ptr c_id;
} gpipe, *gpipe_ptr;

struct node {
  int key;
  int adj_n;
  int nl, nr;
  gpipe_ptr *outg;
  gpipe_ptr *incm;
  gpipe_ptr *adj_p;
  node_ptr *adj_k;
  node_ptr *left;
  node_ptr *right;
  double D, D_prev;
  double P, F, mult, Q, G; 
  double cratio;
};

struct compressor {
  node_ptr position;
  gpipe_ptr tp;
  double cratio;
  double cratio_prev;
};

struct params {
  double time; 
};

typedef struct noise_p {
  double tau;	// correlation time 
  double delta; // noise generator time step 
  double lc;    // correlation distance
  double A;	// noise amplitude
} *noise_ptr;

struct network {
  char fname[80], nname[80], dname[80], current_dir[1024];
  char incname[256];
  int nnodes;
  int nlinks;
  int ncomps;
  int npcent;
  int mglf;
  int nskip;
  int LFLAG;			// LFLAG = 1 -- fixed pressure on inlet, 0 -- fixed flux on inlet
  double c, diss, tmax;
  double curr_T;
  double P0, F0;
  noise_p *noise;
  node_ptr *knot;
  gpipe_ptr *link;
  compressor_ptr cssr;
  noise_ptr nse;
};

// aux.c
extern void err_msg(char* msg);
extern void debug_msg(char* msg);
extern void view_network(network_ptr ns);
extern void copy_network(network_ptr ns, network_ptr nd);

//
extern int noise_status(noise_p *noise);

// init.c
extern void init_data(network *net);
extern void init_links(network *net);
extern void init_arrays(network *net);
extern void allocate_memory(network *net, double *lm);
extern void rescale_data(network *net);
extern int load_data(FILE *fh, gpipe_ptr lnk);
extern void save_data(FILE* fh, gpipe_ptr lnk, network* net, double time);
extern void init_save_balance(network* net);
extern void save_balance(network* net, network* netb, double dt, double dx);

// gsl_noise.c
extern int init_network(char *filename, network *net);

// drawing.c 
/*
extern void mgl_init_draw(network *net);
extern void mgl_draw_pipe(gpipe_ptr in, network *net, char* fname, char* title);
extern void mgl_draw_pressure(network *net, char* fname, char* title);
extern void mgl_draw_flux(network *net, char* fname, char* title);
extern void mgl_draw_network(network *net, char* fname, char* title);
*/
//evolve.c
extern void init_nodes(network *net);
extern void init_evolve(network *net);
extern void hyperbolic_step(network *net, double time);
extern void nonlinear_hstep(network* net, double dt);
extern void evolve_network(network *net);
extern void evolve_network_balance(network *net, network *netb);
extern void split_step2(network *net, double dt, double time);
extern void anatoly_bc(network *net, double time);

//bc.c
extern void init_compressors(network *net);
extern void init_demands(network *net);
extern void update_compressors(network* net, double ctime);
extern void update_demands(network *net, double ctime);

//noise.c
extern void play_with_noise(network_ptr net);
extern void init_noise(network_ptr net);
extern void generate_noise_step();
extern void generate_poor_mans_noise(double time);

//adiabatic.c
extern void prepare_arrays_adiabatic(network_ptr net);
extern void adiabatic_rk4(fftw_complex *in, double dt);
extern void save_adiabatic_temporal();
extern void save_adiabatic();

