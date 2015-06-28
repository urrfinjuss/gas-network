#include "network.h"

static double ovsq2, sq2;
//static double pl, pr, a;
static double a, Q; // consumption at node
static double beta, P_l, P_r, f_l, f_r;
static double **w, **wb;
static int *Np, Nl, tl;
static gpipe_ptr p_k;
static node_ptr n_k;
static double q1;
//static node_ptr n_k, left, right;

void init_evolve(network *net) {
  w  = malloc(sizeof(double *)*(net->nlinks));
  wb = malloc(sizeof(double *)*(net->nlinks));
  Np = malloc((net->nlinks)*sizeof(int));
  ovsq2 = 1./sqrt(2.);	
  sq2 = sqrt(2.);
  //a = (2. + pow(2., 1./3) + pow(2., -1./3))/3.;
  a  = 0.002; 
  Nl = net->nlinks;
  for (int k = 0; k < net->nlinks; k++){
    p_k = net->link[k];
    Np[k] = p_k->N;
    w[k]  = malloc(sizeof(double)*(p_k->N));
    wb[k] = malloc(sizeof(double)*(p_k->N));
    p_k->W_l = (p_k->p[0]+p_k->f[0])*ovsq2;
    p_k->W_r = (p_k->p[Np[k]-1]+p_k->f[Np[k]-1])*ovsq2;
    p_k->Wb_l = (p_k->p[0]-p_k->f[0])*ovsq2;
    p_k->Wb_r = (p_k->p[Np[k]-1]-p_k->f[Np[k]-1])*ovsq2;
  }
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
        n_k->P += p_k->Wb_l;
      } 
      else if ( p_k->right == n_k ) {
	n_k->P += p_k->W_r;
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = sq2*(n_k->P)/(n_k->adj_n);
  }
}

void split_step2(network *net, double dt, double time) {
  nonlinear_hstep(net, 0.5*dt);
  hyperbolic_step(net, time);
  nonlinear_hstep(net, 0.5*dt);
}


void nonlinear_hstep(network* net, double dt) {
  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];
    beta = 0.5*(net->diss)/(p_k->d);
    P_l = (p_k->left)->P;
    P_r = (p_k->right)->P;

    f_l = sq2*(p_k->W_l) - P_l;
    f_r = P_r - sq2*(p_k->Wb_r);
    //printf("Before:\t%e\t%e\t%e\t%e\n", f_l, f_r, P_l, P_r);
    for (int n = 0; n < p_k->N; n++) {
      p_k->f[n] = p_k->f[n]/(1. + beta*dt*fabs(p_k->f[n])/(p_k->p[n]));
    }
    f_l = f_l/(1. + beta*dt*fabs(f_l)/P_l);
    f_r = f_r/(1. + beta*dt*fabs(f_r)/P_r);

    //printf("After:\t%e\t%e\t%e\t%e\n", f_l, f_r, P_l, P_r);

    p_k->W_l = ovsq2*(P_l + f_l);
    p_k->Wb_r = ovsq2*(P_r - f_r);
  }
}


void update_bc(network *net, double time) {
  //  Update of BC below
  for (int j = 1; j < net->nnodes; j++) {
    n_k = net->knot[j];
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
        n_k->P += p_k->Wb_l;
      } 
      else if ( p_k->right == n_k ) {
	n_k->P += p_k->W_r;
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = sq2*(n_k->P)/(n_k->adj_n);
    for (int l = 1; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
	n_k->F = n_k->P - sq2*(p_k->Wb_l);
	p_k->W_l = ovsq2*(n_k->P + n_k->F);
      } 
      else if ( p_k->right == n_k ) {
	n_k->F = sq2*(p_k->W_r) - n_k->P;
	p_k->Wb_r = ovsq2*(n_k->P - n_k->F);
      } else {
        err_msg("Error in network construction");
      }
    }
  }
  //  End Update  
  // Override to have influx at node 0
  n_k = net->knot[0]; 
  n_k->P = 0.; 
  Q = 150.*(1. - 1./cosh(a*time)); 
  for (int l = 0; l < n_k->nr; l++) {
      //printf("%p\n", n_k->outg[l]);
      p_k = n_k->outg[l];
      n_k->P += p_k->Wb_l;
  }
  for (int l = 0; l < n_k->nl; l++) {
      //printf("%p\n", n_k->incm[l]);
      p_k = n_k->incm[l];
      n_k->P += p_k->Wb_l;
  }

  n_k->P = sq2*(n_k->P)/(n_k->adj_n) - Q;
  //printf("%.8e\t%.8e\n", -Q, n_k->P);

  for (int l = 0; l < n_k->nr; l++) {
      p_k = n_k->outg[l];
      n_k->F = n_k->P - sq2*(p_k->Wb_l);
      p_k->W_l = ovsq2*(n_k->P + n_k->F);
  }
  for (int l = 0; l < n_k->nl; l++) {
      p_k = n_k->incm[l];
      n_k->F = sq2*(p_k->W_r) - n_k->P;
      p_k->Wb_r = ovsq2*(n_k->P - n_k->F);
  }


  n_k = net->knot[3]; 
  n_k->P = 0.; 
  Q = -150.*(1. - 1./cosh(a*time)); 
  for (int l = 0; l < n_k->nr; l++) {
      //printf("%p\n", n_k->outg[l]);
      p_k = n_k->outg[l];
      n_k->P += p_k->Wb_l;
  }
  for (int l = 0; l < n_k->nl; l++) {
      //printf("%p\n", n_k->incm[l]);
      p_k = n_k->incm[l];
      n_k->P += p_k->Wb_l;
  }

  n_k->P = sq2*(n_k->P)/(n_k->adj_n) - Q;
  //printf("%.8e\t%.8e\n", -Q, n_k->P);

  for (int l = 0; l < n_k->nr; l++) {
      p_k = n_k->outg[l];
      n_k->F = n_k->P - sq2*(p_k->Wb_l);
      p_k->W_l = ovsq2*(n_k->P + n_k->F);
  }
  for (int l = 0; l < n_k->nl; l++) {
      p_k = n_k->incm[l];
      n_k->F = sq2*(p_k->W_r) - n_k->P;
      p_k->Wb_r = ovsq2*(n_k->P - n_k->F);
  }
  //err_msg("Complete");
}


void hyperbolic_step(network *net, double time) {
  for (int k = 0; k < Nl; k++) {
    tl = Np[k]-1;
    p_k = net->link[k];
    for (int j = 0; j < Np[k]; j++) {
      w[k][j]  = (p_k->p[j] + p_k->f[j])*ovsq2;
      wb[k][j] = (p_k->p[j] - p_k->f[j])*ovsq2;
    }    
    // taking values from characteristics coming from nodes [time (+0)]
    w[k][0] = p_k->W_l;  
    wb[k][tl] = p_k->Wb_r;
    // updating characteristics going into the nodes [time (+1)]
    p_k->W_r = w[k][Np[k]-1];
    p_k->Wb_l = wb[k][0];
    // finished collecting required values from outgoing characteristics
    memmove( w[k]+1,  w[k], sizeof(double)*(Np[k]-1));
    memmove(wb[k], wb[k]+1, sizeof(double)*(Np[k]-1));
    // now we must update p_k->W_l and p_k->Wb_r
    for (int j = 0; j < Np[k]; j++) {
      p_k->p[j] = (w[k][j] + wb[k][j])*ovsq2;
      p_k->f[j] = (w[k][j] - wb[k][j])*ovsq2;
    }
  }
  update_bc(net, time);

  /*
  //  Update of BC below
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
        n_k->P += p_k->Wb_l;
      } 
      else if ( p_k->right == n_k ) {
	n_k->P += p_k->W_r;
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = sq2*(n_k->P)/(n_k->adj_n);
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
	n_k->F = n_k->P - sq2*(p_k->Wb_l);
	p_k->W_l = ovsq2*(n_k->P + n_k->F);
      } 
      else if ( p_k->right == n_k ) {
	n_k->F = sq2*(p_k->W_r) - n_k->P;
	p_k->Wb_r = ovsq2*(n_k->P - n_k->F);
      } else {
        err_msg("Error in network construction");
      }
    }
  }
  //  End Update
  */
}


void evolve_network(network *net) {
  FILE *fh;
  char msg[1024], outdir[1024], msg2[80];
  double dx = 1000./(net->npcent);
  double dt = dx/(net->c);  			// physical time-step (seconds)
  int n_steps = round((net->tmax)/dt);
  //int n_steps = 50;
  int n_skip = 40;
  int n_curr = 0;
  printf("Simulating network for %f hours\nTime step is %f (in secs)\nTotal steps %d\n", (net->tmax)/3600., dt, n_steps);
  while (1) {
    split_step2(net, dx, n_curr*dt);  // note we use dtau = dx here!
    n_curr++;
    if (n_curr % n_skip == 0) {
      printf("t = %.3f sec\n", n_curr*dt);
      for (int n = 0; n < net->nlinks; n++){ 
        sprintf(outdir, "%s/figures_%03d", net->current_dir, n);
        sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/n_skip);
        sprintf(msg2, "t = %.3f sec", n_curr*dt);
        //mgl_draw_pipe(net->link[n], net, msg, msg2);
        sprintf(msg, "%s/pipe_%03d/%s_%03d.txt", net->current_dir, n, net->dname, n_curr/n_skip);
        fh = fopen(msg, "w");
	save_data(fh, net->link[n], net);
        fclose(fh);
      }
      sprintf(outdir, "%s/network", net->current_dir);
      sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/n_skip);
      mgl_draw_network(net, msg, msg2);
    }
    if (n_curr == n_steps) break; 
  }
}









