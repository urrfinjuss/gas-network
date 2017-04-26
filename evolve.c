#include "network.h"

static double ovsq2, sq2;
//static double pl, pr, a;
static double a, Q; // consumption at node
static double beta, P_l, P_r, f_l, f_r;
static double **w, **wb, *p_prev;
static int *Np, Nl, tl;
static gpipe_ptr p_k;
static node_ptr n_k;
static double q1;
static char msg[80];
static char dmsg[512];
//static node_ptr n_k, left, right;

void init_evolve(network *net) {
  w  = malloc(sizeof(double *)*(net->nlinks));
  wb = malloc(sizeof(double *)*(net->nlinks));
  Np = malloc((net->nlinks)*sizeof(int));
  ovsq2 = 1./sqrt(2.);	
  sq2 = sqrt(2.);
  //printf("My square root 2 and 1/square root 2 %.16e\t%.16e\n", sq2, ovsq2);
  //gam = net->gamma;
  //a = (2. + pow(2., 1./3) + pow(2., -1./3))/3.;
  a  = 0.002; 
  Nl = net->nlinks;
  for (int k = 0; k < net->nlinks; k++){
    p_k = net->link[k];
    Np[k] = p_k->N;
    w[k]  = malloc(sizeof(double)*(p_k->N));
    wb[k] = malloc(sizeof(double)*(p_k->N));
  }
}

void split_step2(network *net, double dt, double time) {
#if NOISE
#if POOR_MANS_NOISE
  generate_poor_mans_noise(time);
#else
  generate_noise_step();
#endif
#endif

  nonlinear_hstep(net, 0.5*dt);
  hyperbolic_step(net, time);
  nonlinear_hstep(net, 0.5*dt);
}

void split_step2nonoise(network *net, double dt, double time) {
  
   // original
  /*
  nonlinear_hstep(net, 0.5*dt);
  hyperbolic_step(net, time);
  nonlinear_hstep(net, 0.5*dt);
  */
  //printf("%.12e\n",net->link[0]->Fl); 
  nonlinear_hstep(net, 0.5*dt);
  //printf("%.12e\n",net->link[0]->Fl); 
  hyperbolic_step(net, time+dt);
  //printf("%.12e\n",net->link[0]->Fl); 
  nonlinear_hstep(net, 0.5*dt);
  //printf("%.12e\n",net->link[0]->Fl); 
}


void nonlinear_hstep(network* net, double dt) {
  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];
    beta = 0.5*(net->diss)/(p_k->d);

    //P_l = (p_k->left)->P;  old dumb
    //P_r = (p_k->right)->P;

    //f_l = sq2*(p_k->W_l) - P_l;
    p_k->Fl = sq2*(p_k->W_l) - p_k->Pl;
    //f_r = P_r - sq2*(p_k->Wb_r);
    p_k->Fr = p_k->Pr - sq2*(p_k->Wb_r);
    //printf("Before:\t%e\t%e\t%e\t%e\n", p_k->Fl/net->c, p_k->Fr/net->c, p_k->Pl, p_k->Pr);
    for (int n = 0; n < p_k->N; n++) {
      p_k->f[n] = p_k->f[n]/(1. + beta*dt*fabs(p_k->f[n])/(p_k->p[n]));
    }
    
    //f_l = f_l/(1. + beta*dt*fabs(f_l)/P_l);
    p_k->Fl = p_k->Fl/(1. + beta*dt*fabs(p_k->Fl)/p_k->Pl);
    //f_r = f_r/(1. + beta*dt*fabs(f_r)/P_r);
    p_k->Fr = p_k->Fr/(1. + beta*dt*fabs(p_k->Fr)/p_k->Pr);

    // end of standard nonlinear part
    // distributed consumption
#if NOISE
    p_k->Fl = p_k->Fl*(1. + beta*dt*fabs(p_k->Fl)/p_k->Pl);
    p_k->Fr = p_k->Fr*(1. + beta*dt*fabs(p_k->Fr)/p_k->Pr);
#if LINEAR_BC
    p_k->Fl = p_k->Fl + dt*beta*((p_k->left)->G)*p_k->Pl;
    p_k->Pl = p_k->Pl + ((p_k->left)->Q)*dt;
#endif
    for (int n = 0; n < p_k->N; n++) {
      p_k->f[n] = p_k->f[n] + dt*beta*p_k->gamma[n]*(p_k->p[n]);
      p_k->p[n] = p_k->p[n] + dt*(p_k->q[n] + p_k->fx[n]);
    }
    //printf("%.12e\t%.12e\t%.12e\n", p_k->p[0], p_k->p[1], p_k->p[2]);
#if LINEAR_BC
    p_k->Fr = p_k->Fr + dt*beta*((p_k->right)->G)*p_k->Pr;
    p_k->Pr = p_k->Pr + ((p_k->right)->Q)*dt;
#endif
#endif
    //printf("After:\t%e\t%e\t%e\t%e\n", p_k->Fl/net->c, p_k->Fr/net->c, p_k->Pl, p_k->Pr);

    //p_k->W_l = ovsq2*(P_l + f_l);
    p_k->W_l = ovsq2*(p_k->Pl + p_k->Fl);
    //p_k->Wb_r = ovsq2*(P_r - f_r); 
    p_k->Wb_r = ovsq2*(p_k->Pr - p_k->Fr);

    // why is there no update for other 2 characteristics?
    //p_k->W_r = ovsq2*(P_r + f_r);
    p_k->W_r = ovsq2*(p_k->Pr + p_k->Fr);
    //p_k->Wb_l = ovsq2*(P_l - f_l); 
    p_k->Wb_l = ovsq2*(p_k->Pl - p_k->Fl); 
  }
}


void update_bc(network *net, double time) {
  double dt = 1000./(net->npcent);
  double Beta; 
  //  Update of BC below
  update_compressors(net, time);
  update_demands(net, time);
  double flux_residual = 0.;
  //double Sinc;
  for (int j = 0; j < net->nnodes; j++) {
    n_k = net->knot[j];
    n_k->P = 0.;//-(net->c)*n_k->D;
    n_k->mult = 0.;
    sprintf(dmsg,"Doing Node %p:\n", n_k);
    debug_msg(dmsg);
    //printf("Node %02d:\t", j);
    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            	      // means left side of pipe p_k is connected to node j
        n_k->P += p_k->Wb_l*0.25*pi*(p_k->d*p_k->d);  // new accounts for distinct crossections of pipes

	//printf("found pipe %p (left)\t", p_k);

	sprintf(dmsg,"Pipe %p: W* = %.15e\n", p_k, p_k->Wb_l);
	debug_msg(dmsg);
#if 1	
	if (p_k->c_id != NULL) {
		n_k->mult += (p_k->c_id)->cratio*0.25*pi*(p_k->d*p_k->d);
		sprintf(msg, "Compressor %p located in pipe %p \n", p_k->c_id, p_k);
		debug_msg(msg);
	}
	else n_k->mult += 0.25*pi*(p_k->d*p_k->d);
#endif
      } 
      else if ( p_k->right == n_k ) { 		      // means right side of pipe p_k is connected to node j
	n_k->P += p_k->W_r*0.25*pi*(p_k->d*p_k->d);   // new

	//printf("found pipe %p (right)\t Diameter = %.12e\n", p_k, p_k->d); exit(1);

	sprintf(dmsg, "Pipe %p: W = %.15e\n", p_k, p_k->W_r);
	debug_msg(dmsg);
	n_k->mult += 0.25*pi*(p_k->d*p_k->d);
      } else {
        err_msg("Error in network construction");
      }
    }
    //printf("Sum Si*Wi = %.16e\nn_k->mult = %.16e\nC = %.16e\n", n_k->P, n_k->mult,(p_k->c_id)->cratio);
    //exit(1);
    sprintf(dmsg,"Current denominator %.15e\n", n_k->mult);
    debug_msg(dmsg);

    
    //unsigned int sim = 3;
    if (j == 0) {
	if (net->LFLAG == 1) {
	  n_k->P = net->P0;
  	  //if (p_k->c_id != NULL) n_k->P = net->P0*(p_k->c_id)->cratio;
	} else {
	  n_k->P = (sq2*n_k->P)/(n_k->mult) + net->F0;  // should derive this again... n_k->P = (sq2*n_k->P + (net->F0))/(n_k->mult);
	}
        /*
        if (sim == 2) {
          double P0 = 1e7;
          double Pval = P0;
          double frq = 0.00174436248289312;

          Pval = P0*exp(144.*1e-9/frq*(1.-cos(1.0*frq*time)));
          n_k->P = Pval;
        } else if (sim == 3) {
	  double A = 10088085.7094053;
	  double B = -88289.1640032395;
	  double frq = 0.00174534439987236;
          double Pval = A;
 	  Pval += B*cos(1.0*frq*time);
 	  n_k->P = Pval;
        }
       */
    } else n_k->P = (sq2*n_k->P - (net->c)*(n_k->D))/(n_k->mult);

    if (j == 1) { 
      if (RFLAG_P_FIXED) {
      /*
 	if (sim == 1) {
          double A1 = 6.40365e+06;
          double B1 = 2.80523e+06;
          double C1 = 324426;
          double D1 = 77302.9;
          double E1 = 22734.1;
          double frq = 0.00174436248289312;

	  double Pval = A1;
          Pval += -B1*cos(1.0*frq*time);
          Pval += -C1*cos(2.0*frq*time);  
          Pval += -D1*cos(3.0*frq*time);  
          Pval += -E1*cos(4.0*frq*time);  
          n_k->P = Pval;
        } else if (sim == 2) {
          double A1 = 3156065.06937776;
          double B1 = 277911.435298177;
          double C1 = -6718.98478372491;
          double D1 = 65.6104975076953;
          double frq = 0.00174436248289312;
          
          double Pval = A1;
          Pval += B1*cos(1.0*frq*time);
          Pval += C1*cos(2.0*frq*time);
          Pval += D1*cos(3.0*frq*time);
          n_k->P = Pval;
        } else if (sim == 3) {
	  double frq = 0.00174534439987236;
          double A1 = 7.79828e+06;
          double B1 = 428934;
          double B2 = -68367.4;
          double C1 = -1846.93;
          double C2 =  16060.2;
          double Pval = A1;
          Pval += B1*sin(1.0*frq*time) + B2*cos(1.0*frq*time);
          Pval += C1*sin(2.0*frq*time) + C2*cos(2.0*frq*time);
          n_k->P = Pval;
        }
      */
      }
    }


    for (int l = 0; l < n_k->adj_n; l++) {
      p_k = n_k->adj_p[l];
      if (p_k->left == n_k) {            // means left side of pipe p_k is connected to node j
#if 1       
	if (p_k->c_id != NULL) {
           sprintf(dmsg,"Left side is connected here\n");
	   debug_msg(dmsg);

	   if (j == 0) {
	      if (net->LFLAG==1) {
		n_k->F = ((p_k->c_id)->cratio)*n_k->P - sq2*(p_k->Wb_l);
        	//printf("Flux = %e\n",n_k->F/net->c);
	      }
	      else {}//n_k->F = net->F0;
	   } else n_k->F = ((p_k->c_id)->cratio)*n_k->P - sq2*(p_k->Wb_l);
           
           p_k->Pl = (p_k->c_id->cratio)*n_k->P;  // experimental 1
	   //p_k->Fl = (n_k->F)/net->c; old dumb
	   p_k->Fl = n_k->F;
	   //p_k->W_l = ovsq2*(((p_k->c_id)->cratio)*n_k->P + n_k->F);
	   p_k->W_l = ovsq2*(p_k->Pl + p_k->Fl);

	   sprintf(dmsg, "To pipe %p: flux is %.15e\tpressure is %.15e\tCompression = %.12e\n", p_k, (p_k->Fl), 0.000001*((p_k->c_id)->cratio)*n_k->P, (p_k->c_id)->cratio);
	   debug_msg(dmsg);
        } else {
	   if (j == 0) {
	      if (net->LFLAG == 1) {
		 n_k->F = n_k->P - sq2*(p_k->Wb_l);
	      } else {
		 n_k->F = net->F0;
		 //Beta = 0.5*(net->diss)/(p_k->d);
		 //n_k->F = p_k->f[1]/(1 + 0.5*Beta*dt*p_k->f[1]/p_k->p[1]) + (n_k->P - p_k->p[1]) - 0.5*dt*(p_k->q[1]+n_k->Q);
	      }
	   } else {
		n_k->F = n_k->P - sq2*(p_k->Wb_l);
           }
	   p_k->Pl = n_k->P;  // experimental 2
	   p_k->Fl = (n_k->F); //net->c;
           p_k->W_l = ovsq2*(p_k->Pl + p_k->Fl);
	   sprintf(dmsg, "To pipe %p: flux is %.15e\tpressure is %.15e\n", p_k, (p_k->Fl), 0.000001*n_k->P);
	   debug_msg(dmsg);
	}
#endif
      } else if ( p_k->right == n_k ) {
	n_k->F = sq2*(p_k->W_r) - n_k->P;
	p_k->Pr = n_k->P;  // experimental 3
	//p_k->Fr = (n_k->F)/net->c; // old dumb
	p_k->Fr = (n_k->F);
	p_k->Wb_r = ovsq2*(p_k->Pr - p_k->Fr);
        sprintf(dmsg, "To pipe %p: flux is %.15e\tpressure is %.15e\n", p_k, (p_k->Fr), 0.000001*n_k->P);
	debug_msg(dmsg);
      } else err_msg("Error in network construction");
    }
    //printf("Node %02d Complete\n", j);   
  }
}


void hyperbolic_step(network *net, double time) {
  for (int k = 0; k < Nl; k++) {
    tl = Np[k]-1;
    p_k = net->link[k];
    //printf("%.16e\t%.16e\n", p_k->p[460], p_k->f[460]/net->c);
    //exit(1);
    for (int j = 0; j < Np[k]; j++) {
      w[k][j]  = (p_k->p[j] + p_k->f[j])*ovsq2;
      wb[k][j] = (p_k->p[j] - p_k->f[j])*ovsq2;
    }
    //printf("%e\n%e\n%e\n",(p_k->W_l - p_k->Wb_l)*ovsq2/net->c, (w[0][0] - wb[0][0])*ovsq2/net->c, (w[0][1] - wb[0][1])*ovsq2/net->c);
    //if (p_k->c_id != NULL)     
    // taking values from characteristics coming from nodes [time (+0)]
    // updating characteristics going into the nodes [time (+1)]
    p_k->W_r = w[k][Np[k]-1];
    p_k->Wb_l = wb[k][0];
    /*if (p_k->c_id != NULL) {
	sprintf(dmsg,"Multiplier at compressor %p is %.4e\n", p_k->c_id, ((p_k->c_id)->cratio));
	debug_msg(dmsg);
    	p_k->Wb_l = (p_k->p[0] - p_k->f[0])*ovsq2;  // *(p_k->c_id->cratio)
    }
    else p_k->Wb_l = wb[k][0];*/
    //printf("step: Wbl = %.12e\nP = %.12e\tF = %.12e\n", p_k->Wb_l, p_k->p[0], p_k->f[0]/net->c);
    // finished collecting required values from outgoing characteristics
    
    memmove( w[k]+1,  w[k], sizeof(double)*(Np[k]-1));
    w[k][0] = p_k->W_l;  
    memmove(wb[k], wb[k]+1, sizeof(double)*(Np[k]-1));
    wb[k][tl] = p_k->Wb_r;
    
    // alternative straight
    /*for (int j = tl-1; j > -1; j--) {
	w[k][j+1] = w[k][j];
    }
    w[k][0] = p_k->W_l;  
    for (int j = 0; j < tl; j++) {
	wb[k][j] = wb[k][j+1];
    }
    wb[k][tl] = p_k->Wb_r;*/
    // end alternative straight


    for (int j = 0; j < Np[k]; j++) {
      p_k->p[j] = (w[k][j] + wb[k][j])*ovsq2;
      p_k->f[j] = (w[k][j] - wb[k][j])*ovsq2;
    }
  }
  update_bc(net, time);
  for (int k = 0; k < Nl; k++) {
    tl = Np[k]-1;
    p_k = net->link[k];
    p_k->Fm = (w[k][Np[k]/2] - wb[k][Np[k]/2])*ovsq2;
    p_k->Pm = (w[k][Np[k]/2] + wb[k][Np[k]/2])*ovsq2;
  }
}

void init_dump(network_ptr net) {
  FILE *fhtime = fopen("temporal_pressure.txt","w");
  int N = 0;
  gpipe_ptr pk;

  p_prev = malloc((net->link[0]->N)*sizeof(double));


  fprintf(fhtime, "# 1. time, sec ");
  for (int n = 0; n < net->nnodes; n++) fprintf(fhtime, "%d. Node %d: Pressure, MPa\t", n+2, n);
  fprintf(fhtime, "\n\n%.12e\t", 0.);
  for (int n = 0; n < net->nlinks; n++) {
 	pk = net->link[n];
	N = pk->N;
	fprintf(fhtime, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t",  
	0.000001*pk->Pl,
        0.000001*pk->p[N/8],  
	0.000001*pk->p[N/4],
        0.000001*pk->p[3*N/8],
	0.000001*pk->p[N/2],
        0.000001*pk->p[5*N/8],
	0.000001*pk->p[3*N/4],  
	0.000001*pk->p[7*N/8],
	0.000001*pk->Pr );
  }
  //for (int n = 0; n < net->nnodes; n++) fprintf(fhtime, "%e\t",  0.000001*((net->knot[n])->P) );
  fprintf(fhtime, "\n"); 
  fclose(fhtime);

  fhtime = fopen("temporal_flux.txt","w");
  fprintf(fhtime, "# 1. time, sec ");
  for (int n = 0; n < net->nlinks; n++) fprintf(fhtime, "%d. Link %d left: Flux, kg/m^2/s\t", n+2, n);
  fprintf(fhtime, "\n\n%.12e\t", 0.);
  for (int n = 0; n < net->nlinks; n++) fprintf(fhtime, "%.12e\t%.12e\t%.12e\t", (net->link[n]->Fl)/net->c, (net->link[n]->Fm)/net->c, (net->link[n]->Fr)/net->c ); 
  //ovsq2*((net->link[n])->W_l - (net->link[n])->Wb_l)/(net->c)
  //for (int n = 0; n < net->nlinks; n++) fprintf(fhtime, "%e\t", (net->link[n]->Fl)/net->c ); //ovsq2*((net->link[n])->W_l - (net->link[n])->Wb_l)/(net->c)
  fprintf(fhtime, "\n");  
  fclose(fhtime);


}

void write_logdiffp(char *fname, network_ptr net) {
   
}

void dump_current(network_ptr net, double time) {
  FILE *fhtime = fopen("temporal_flux.txt", "a");
  int N = 0;
  gpipe_ptr pk;

  fprintf(fhtime, "%.12e\t", time);  
  for (int n = 0; n < net->nlinks; n++) {
    fprintf(fhtime, "%.12e\t%.12e\t%.12e\t", (net->link[n]->Fl)/net->c, (net->link[n]->Fm)/net->c, (net->link[n]->Fr)/net->c );
    //fprintf(fhtime, "%e\t", (net->link[n]->Fl)/net->c  );
  }
  fprintf(fhtime, "\n");  
  fclose(fhtime);

  fhtime = fopen("temporal_pressure.txt", "a");
  fprintf(fhtime, "%.12e\t", time);      
  //for (int n = 0; n < net->nnodes; n++) fprintf(fhtime, "%e\t",  0.000001*((net->knot[n])->P) );
  //for (int n = 0; n < net->nlinks; n++) fprintf(fhtime, "%e\t",  0.000001*((net->link[n])->Pl) );
  for (int n = 0; n < net->nlinks; n++) {
	pk = net->link[n];
	N = pk->N;
	fprintf(fhtime, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t",  
        0.000001*pk->Pl,
        0.000001*pk->p[N/8],
        0.000001*pk->p[N/4],
        0.000001*pk->p[3*N/8],
        0.000001*pk->Pm,
        0.000001*pk->p[5*N/8],
        0.000001*pk->p[3*N/4],
        0.000001*pk->p[7*N/8],
        0.000001*pk->Pr, 
	creal(p_k->fx[N/2]));
  }
  fprintf(fhtime, "\n");  
  fclose(fhtime);
}

void evolve_network(network *net) {
  FILE *fh, *fhtime;
  char msg[1024], outdir[1024], msg2[80];
  double dx = 1000./(net->npcent);
  double dt = dx/(net->c);  			// physical time-step (seconds)
  int n_steps = round((net->tmax)/dt);
  //int n_steps = 16;
  int n_skip = net->nskip;
  //int n_skip = 1;
  int n_curr = 0;
  printf("Simulating network for %f hours\nTime step is %f (in secs)\nTotal steps %d\n", (net->tmax)/3600., dt, n_steps);
  
  init_dump(net);
  
  while (1) {
#if NOISE
    split_step2(net, dx, n_curr*dt);  // note (dx = dtau)
    //printf("Uncomment to run with noise\n"); exit(1);	
#else
    memcpy(p_prev, net->link[0]->p, (net->link[0]->N)*sizeof(double));
    split_step2nonoise(net, dx, n_curr*dt);  // note (dx = dtau)
#endif
    if (ADIABATIC) adiabatic_rk4(net->link[0]->fx, dt);
    n_curr++;
    
    if (n_curr % n_skip == 0) {
      printf("t = %.8e min\n", n_curr*dt/60.);
      dump_current(net, n_curr*dt);
      save_adiabatic();

      if (ADIABATIC_SAVE) save_adiabatic_temporal();
      if (PIPE_SAVE) {	
        for (int n = 0; n < net->nlinks; n++){ 
          sprintf(msg, "%s/pipe_%03d/%s_%03d.txt", net->current_dir, n, net->dname, n_curr/n_skip);
          fh = fopen(msg, "w");
	  save_data(fh, net->link[n], net, n_curr*dt);
          fclose(fh);         
        }
        /*if (net->mglf != 0) {
          sprintf(outdir, "%s/network", net->current_dir);
          sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/n_skip);
          sprintf(msg2, "t = %.3f sec", n_curr*dt);
          mgl_draw_network(net, msg, msg2);
        }*/
      }
    }
    if (n_curr == n_steps) break; 
  }
  printf("t = %.8e min\n", n_curr*dt/60.);
  fhtime = fopen("params.txt", "a");
  fprintf(fhtime, "%e\t%e\t%e\n", n_curr*dt/3600., ((net->knot[1])->P)*0.000001, ovsq2*((net->link[0])->W_l - (net->link[0])->Wb_l)/(net->c) );
  fclose(fhtime);
}

void evolve_network_balance(network *net, network *netb){
  FILE *fh, *fhtime;
  char msg[1024], outdir[1024], msg2[80];
  double dx = 1000./(net->npcent);
  double dt = dx/(net->c);  			// physical time-step (seconds)
  int n_steps = round((net->tmax)/dt);
  int n_skip = net->nskip;
  int n_curr = 0; net->curr_T = n_curr*dt; 
  printf("Simulating network for %f hours\nTime step is %f (in secs)\nTotal steps %d\n", (net->tmax)/3600., dt, n_steps);
  init_save_balance(net);
  while (1) {
    copy_network(net, netb);
    split_step2(net, dx, n_curr*dt);  // note (dx = dtau)
    n_curr++; net->curr_T = n_curr*dt;
    
    if (n_curr % n_skip == 0) {
      printf("t = %.3f min\n", n_curr*dt/60.);
      save_balance(net, netb, dt, dx);
      fhtime = fopen("temporal_flux.txt", "a");
      fprintf(fhtime, "%e\t", n_curr*dt);  
      for (int n = 0; n < net->nlinks; n++) {
	fprintf(fhtime, "%e\t%e\t%e\t", (net->link[n]->Fl), (net->link[n]->Fm), (net->link[n]->Fr)  );
      }
      fprintf(fhtime, "\n");  
      fclose(fhtime);

      fhtime = fopen("temporal_pressure.txt", "a");
      fprintf(fhtime, "%e\t", n_curr*dt);  
      for (int n = 0; n < net->nlinks; n++) {
        fprintf(fhtime, "%e\t%e\t%e\t",  0.000001*((net->link[n])->Pl),  0.000001*((net->link[n])->Pm),  0.000001*((net->link[n])->Pr) );
      } 
      fprintf(fhtime, "\n");  
      fclose(fhtime);

      for (int n = 0; n < net->nlinks; n++){ 
        /*sprintf(outdir, "%s/figures_%03d", net->current_dir, n);
        sprintf(msg, "%s/%s_%03d.png", outdir, net->dname, n_curr/n_skip);
        sprintf(msg2, "t = %.3f sec", n_curr*dt);
        mgl_draw_pipe(net->link[n], net, msg, msg2);*/
        sprintf(msg, "%s/pipe_%03d/%s_%03d.txt", net->current_dir, n, net->dname, n_curr/n_skip);
        fh = fopen(msg, "w");
	save_data(fh, net->link[n], net, n_curr*dt);
        fclose(fh);
      }
    }
    if (n_curr == n_steps) break; 
  }
  printf("t = %.3f min\n", n_curr*dt/60.);
}









