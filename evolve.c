#include "network.h"

static double ovsq2, pl, pr;
static double **w, **wb;
static int *Np, Nl, tl;
static gpipe_ptr p_k;
static node_ptr left, right;

void init_evolve(network *net) {
  w  = malloc(sizeof(double *)*(net->nlinks));
  wb = malloc(sizeof(double *)*(net->nlinks));
  Np = malloc((net->nlinks)*sizeof(int));
  for (int k = 0; k < net->nlinks; k++){
    p_k = net->link[k];
    Np[k] = p_k->N;
    w[k]  = malloc(sizeof(double)*(p_k->N));
    wb[k] = malloc(sizeof(double)*(p_k->N));
  }
  Nl = net->nlinks;
  ovsq2 = 1./sqrt(2.);	
}

void hyperbolic_step(network *net) {
  for (int k = 0; k < Nl; k++) {
    tl = Np[k]-1;
    p_k = net->link[k];
    left = (p_k->left);
    right = (p_k->right);
    left->P  = 	p_k->p[0];	// something simple first
    right->P = 	p_k->p[tl];	// something simple first
    for (int j = 0; j < Np[k]; j++) {
      w[k][j]  = (p_k->p[j] + p_k->f[j])*ovsq2;
      wb[k][j] = (p_k->p[j] - p_k->f[j])*ovsq2;
    }
    //printf("At t = 0:\t%e\t%e\t%e\t%e\n", w[k][0], w[k][1], wb[k][0], wb[k][1]);
    memmove( w[k]+1,  w[k], sizeof(double)*(Np[k]-1));
    memmove(wb[k], wb[k]+1, sizeof(double)*(Np[k]-1));
    /*printf("At t = dt:\t%e\t%e\t%e\t%e\n", w[k][0], w[k][1], wb[k][0], wb[k][1]);
    err_msg("Check Complete");*/
    /*w[k][0] =  
    wb[k][Np[k]-1] = */

    p_k->f[0] = left->P - p_k->p[1] + p_k->f[1];
    p_k->f[tl] = p_k->p[tl-1] + p_k->f[tl-1] - right->P;
    for (int j = 1; j < Np[k]-1; j++) {
      p_k->p[j] = (w[k][j] + wb[k][j])*ovsq2;
      p_k->f[j] = (w[k][j] - wb[k][j])*ovsq2;
    }
    p_k->p[Np[k]-1] = right->P;
    p_k->p[0] = left->P; 	
  }
}

void evolve_network(network *net) {
  char msg[1024], outdir[1024], msg2[80];
  double dx = 1000./(net->npcent);
  double dt = dx/(net->c);  			// physical time-step (seconds)
  int n_steps = round((net->tmax)/dt);
  int n_curr = 0;
  printf("Simulating network for %f hours\nTime step is %f (in secs)\nTotal steps %d\n", (net->tmax)/3600., dt, n_steps);
  while (1) {
    hyperbolic_step(net);
    n_curr++;
    for (int n = 0; n < net->nlinks; n++){ 
      sprintf(outdir, "%s/figures_%03d", net->current_dir, n);
      sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr);
      sprintf(msg2, "t = %.3f sec", n_curr*dt);
      mgl_draw_pipe(net->link[n], net, msg, msg2);
    }
    if (n_curr == n_steps) break; 
  }

}









