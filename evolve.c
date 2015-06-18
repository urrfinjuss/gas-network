#include "network.h"

static double ovsq2, sq2, pl, pr;
static double **w, **wb;
static int *Np, Nl, tl;
static gpipe_ptr p_k;
static node_ptr n_k, left, right;

void init_evolve(network *net) {
  w  = malloc(sizeof(double *)*(net->nlinks));
  wb = malloc(sizeof(double *)*(net->nlinks));
  Np = malloc((net->nlinks)*sizeof(int));
  ovsq2 = 1./sqrt(2.);	
  sq2 = sqrt(2.);
  Nl = net->nlinks;
  for (int k = 0; k < net->nlinks; k++){
    p_k = net->link[k];
    Np[k] = p_k->N;
    w[k]  = malloc(sizeof(double)*(p_k->N));
    wb[k] = malloc(sizeof(double)*(p_k->N));
    p_k->W_l = (p_k->p[0]+p_k->f[0])*ovsq2;
    p_k->W_r = (p_k->p[Np[k]-1]+p_k->f[Np[k]-1])*ovsq2;
    //printf("W left %e\n", p_k->W_l);
    p_k->Wb_l = (p_k->p[0]-p_k->f[0])*ovsq2;
    p_k->Wb_r = (p_k->p[Np[k]-1]-p_k->f[Np[k]-1])*ovsq2;
  }
  //init_nodes(net);
}

/*void init_nodes(network *net) {
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      tl = (p_k->N) - 1;
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
	n_k->W[l]  = 0.;
        n_k->Wb[l] = p_k->p[0];
      } 
      else if ( p_k->right == n_k ) {
        n_k->W[l]  = p_k->p[tl];
	n_k->Wb[l] = 0.;
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) n_k->P += (n_k->W[l] + n_k->Wb[l]);
    n_k->P = (n_k->P)/(n_k->adj_n);
  }
}*/

/*void init_nodes(network *net) {
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      tl = (p_k->N) - 1;		 // index of the last point in a pipe
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
	n_k->W[l]  = (p_k->p[0] + p_k->f[0])*ovsq2;
        n_k->Wb[l] = (p_k->p[0] - p_k->f[0])*ovsq2;
        //printf("Outgoing %p Value on W\t%e\n", p_k, p_k->W_l);
      } 
      else if ( p_k->right == n_k ) {
        n_k->W[l]  = (p_k->p[tl] + p_k->f[tl])*ovsq2;
	n_k->Wb[l] = (p_k->p[tl] - p_k->f[tl])*ovsq2;
        //printf("Incoming %p Value on W+\t%e\n", p_k, p_k->Wb_r);
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) n_k->P += (n_k->W[l] + n_k->Wb[l]);
    n_k->P = (n_k->P)/(n_k->adj_n);
  }
}*/



/*void sync_nodes(network *net) {
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    printf("Node %d:\n", j);
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
	p_k->W_l = n_k->W[l];
        //printf("Outgoing %p Value on W\t%.10e\n", p_k, p_k->W_l);
      } 
      else if ( p_k->right == n_k ) {
	p_k->Wb_r = n_k->Wb[l];
        //printf("Incoming %p Value on W+\t%.10e\n", p_k, p_k->Wb_r);
      } else {
        err_msg("Error in network construction");
      }
    }
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      tl = (p_k->N) - 1;
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
        //n_k->W[l]  = (p_k->p[0] + p_k->f[0])*ovsq2;
	//n_k->W[l]  = 0.;
        n_k->Wb[l] = (p_k->p[0] - p_k->f[0])*ovsq2;
        printf("Outgoing %p: W = %e\tW+ = %e\n", p_k, n_k->W[l], n_k->Wb[l]);
        //printf("P = %e\tF = %e\n", p_k->p[0], p_k->f[0]);
      } 
      else if ( p_k->right == n_k ) {
        n_k->W[l]  =  (p_k->p[tl] + p_k->f[tl])*ovsq2;
        //n_k->Wb[l] =  (p_k->p[tl] - p_k->f[tl])*ovsq2;
	//n_k->Wb[l] = 0.;
        printf("Incoming %p: W = %e\tW+ = %e\n", p_k, n_k->W[l], n_k->Wb[l]);
      } else {
        err_msg("Error in network construction");
      }
    }
    //n_k->P = 0.;
    //for (int l = 0; l < n_k->adj_n; l++) n_k->P += (n_k->W[l] + n_k->Wb[l]);
    //n_k->P = (n_k->P)/((n_k->adj_n)*ovsq2);
  }
}*/


/*void sync_nodes(network *net) {
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    printf("Node %d:\n", j);
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
	p_k->W_l = n_k->W[l];
        //printf("Outgoing %p Value on W\t%.10e\n", p_k, p_k->W_l);
      } 
      else if ( p_k->right == n_k ) {
	p_k->Wb_r = n_k->Wb[l];
        //printf("Incoming %p Value on W+\t%.10e\n", p_k, p_k->Wb_r);
      } else {
        err_msg("Error in network construction");
      }
    }
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      tl = (p_k->N) - 1;
      if (p_k->left == n_k) {            // means pipe is connected to node j on the left
        //n_k->W[l]  = (p_k->p[0] + p_k->f[0])*ovsq2;
	//n_k->W[l]  = 0.;
        n_k->Wb[l] = (p_k->p[0] - p_k->f[0])*ovsq2;
        printf("Outgoing %p: W = %e\tW+ = %e\n", p_k, n_k->W[l], n_k->Wb[l]);
        //printf("P = %e\tF = %e\n", p_k->p[0], p_k->f[0]);
      } 
      else if ( p_k->right == n_k ) {
        n_k->W[l]  =  (p_k->p[tl] + p_k->f[tl])*ovsq2;
        //n_k->Wb[l] =  (p_k->p[tl] - p_k->f[tl])*ovsq2;
	//n_k->Wb[l] = 0.;
        printf("Incoming %p: W = %e\tW+ = %e\n", p_k, n_k->W[l], n_k->Wb[l]);
      } else {
        err_msg("Error in network construction");
      }
    }
    //n_k->P = 0.;
    //for (int l = 0; l < n_k->adj_n; l++) n_k->P += (n_k->W[l] + n_k->Wb[l]);
    //n_k->P = (n_k->P)/((n_k->adj_n)*ovsq2);
  }
}*/


void hyperbolic_step(network *net) {
  for (int k = 0; k < Nl; k++) {
    tl = Np[k]-1;
    p_k = net->link[k];
    //left = (p_k->left);
    //right = (p_k->right);
    //left->P  = 	p_k->p[0];	// something simple first
    //right->P = 	p_k->p[tl];	// something simple first
    //printf("Pipe %d:\n",k);
    for (int j = 0; j < Np[k]; j++) {
      w[k][j]  = (p_k->p[j] + p_k->f[j])*ovsq2;
      wb[k][j] = (p_k->p[j] - p_k->f[j])*ovsq2;
    }    

    // taking values from characteristics coming from nodes [time (+0)]
    w[k][0] = p_k->W_l;  // error here
    wb[k][tl] = p_k->Wb_r;
    // updating characteristics going into the nodes [time (+1)]
    p_k->W_r = w[k][Np[k]-1];
    p_k->Wb_l = wb[k][0];
    // finished collecting required values from outgoing characteristics

    //printf("W_left  \t%.10e\n", p_k->W_l);
    //printf("Wb_right\t%.10e\n", p_k->Wb_r);
    memmove( w[k]+1,  w[k], sizeof(double)*(Np[k]-1));
    memmove(wb[k], wb[k]+1, sizeof(double)*(Np[k]-1));

    // now we must update p_k->W_l and p_k->Wb_r


    //p_k->f[0] = left->P - p_k->p[1] + p_k->f[1];
    //p_k->f[tl] = p_k->p[tl-1] + p_k->f[tl-1] - right->P;
    for (int j = 0; j < Np[k]; j++) {
      p_k->p[j] = (w[k][j] + wb[k][j])*ovsq2;
      p_k->f[j] = (w[k][j] - wb[k][j])*ovsq2;
    }
    //p_k->p[Np[k]-1] = right->P;
    //p_k->p[0] = left->P; 	
  }

  

  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    //printf("Node %d:\n", j);  
    n_k->P = 0.;
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
        n_k->P += p_k->Wb_l;
	//p_k->W_l = ;
        //printf("Outgoing %p Value on W\t%.10e\n", p_k, p_k->W_l);
      } 
      else if ( p_k->right == n_k ) {
	n_k->P += p_k->W_r;
	//p_k->Wb_r = ;
        //printf("Incoming %p Value on W+\t%.10e\n", p_k, p_k->Wb_r);
      } else {
        err_msg("Error in network construction");
      }
    }
    n_k->P = sq2*(n_k->P)/(n_k->adj_n);
    //printf("Pressure = %.8e\n", n_k->P);
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
	n_k->F = n_k->P - sq2*(p_k->Wb_l);
	p_k->W_l = ovsq2*(n_k->P + n_k->F);
        //printf("Outgoing %p W\t%.10e\tF = %.10e\n", p_k, p_k->W_l, n_k->F);
      } 
      else if ( p_k->right == n_k ) {
	n_k->F = sq2*(p_k->W_r) - n_k->P;
	p_k->Wb_r = ovsq2*(n_k->P - n_k->F);
	//p_k->Wb_r = ;
        //printf("Incoming %p W+\t%.10e\tF = %.10e\n", p_k, p_k->Wb_r, n_k->F);
      } else {
        err_msg("Error in network construction");
      }
    }
  }

}

void evolve_network(network *net) {
  char msg[1024], outdir[1024], msg2[80];
  double dx = 1000./(net->npcent);
  double dt = dx/(net->c);  			// physical time-step (seconds)
  int n_steps = round((net->tmax)/dt);
  //int n_steps = 4;
  int n_curr = 0;
  printf("Simulating network for %f hours\nTime step is %f (in secs)\nTotal steps %d\n", (net->tmax)/3600., dt, n_steps);
  while (1) {
    //sync_nodes(net);
    hyperbolic_step(net);
    n_curr++;
    if (n_curr % 1 == 0) {
      for (int n = 0; n < net->nlinks; n++){ 
        sprintf(outdir, "%s/figures_%03d", net->current_dir, n);
        sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/1);
        sprintf(msg2, "t = %.3f sec", n_curr*dt);
        //mgl_draw_pipe(net->link[n], net, msg, msg2);


        sprintf(outdir, "%s/network", net->current_dir);
        sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/1);
        mgl_draw_network(net, msg, msg2);
      }
    }
    if (n_curr == n_steps) break; 
  }

}









