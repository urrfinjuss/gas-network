#include "network.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

static gpipe_ptr p_k;
static node_ptr n_k;
static double ovsq2, sq2, cs;
static double *time;
static double **c2D;
static gsl_interp_accel *acc;
static gsl_spline **spline_arr;


void init_compressors(network *net){
   FILE *fh = fopen("c.txt","r");
   size_t count = 1000;
   char *line = malloc(1000);
   char msg[256];
   int rnum, m = 0;
   int lcount = 0;

   sq2 = sqrt(2.);
   ovsq2 = 1./sq2;
   cs = net->c;
   fgets(msg, 256, fh);
   fgets(msg, 256, fh);
   while(1){
	//printf("%d\n", lcount);
	if(fgets(msg, 256, fh)==NULL) break;
	lcount++;
   }
   sprintf(msg, "\nTemporal data available for %8f hours\n", 1.*(lcount-1)/12.);
   debug_msg(msg);
   rewind(fh);
   fgets(msg, 256, fh);
   fgets(msg, 256, fh);

   double *link_matrix = malloc(sizeof(double)*(net->ncomps + 1)*lcount);
   //err_msg("Complete");

   getline(&line, &count, fh);
   int read = -1, cur = 0, cCount = 0;
   while( sscanf(line+cur, "%lf%n", &link_matrix[cCount], &read) == 1) {
	cur+=read;
	cCount++;
   }
   int rCount = 1;
   while(getline(&line, &count, fh)!=-1) {
	rCount++;
   }
   rewind(fh);
   printf("%d\n%d\n%d\n", rCount, cCount, net->ncomps);
   int i = 0;
   while(getline(&line, &count, fh)!=-1)
   {
     read = -1;
     cur  =  0;
     while(sscanf(line+cur, "%lf%n", &link_matrix[i], &read) == 1) {
	cur+=read;
	i=i+1;
     }
   }
   fclose(fh);
   sprintf(msg, "Compression data has %d entries\n", cCount*rCount);
   debug_msg(msg);

   time = malloc(lcount*sizeof(double));
   c2D = malloc(net->ncomps*sizeof(double *));
   for (int i = 0; i < net->ncomps; i++) c2D[i] = malloc(lcount*sizeof(double));
   for (int j = 0; j < lcount; j++) time[j] = link_matrix[6*j];
   for (int i = 0; i < net->ncomps; i++) {
	for (int j = 0; j < lcount; j++) {
	   c2D[i][j] = link_matrix[6*j+i+1];
	}
   }
   acc = gsl_interp_accel_alloc();
   const gsl_interp_type *t = gsl_interp_cspline_periodic;
   spline_arr = malloc(net->ncomps*sizeof(gsl_spline *));
   for (int i = 0; i < net->ncomps; i++) spline_arr[i] = gsl_spline_alloc(t, lcount);

#if 1
   //double xi, yi;
   //char str[24];
   //FILE *fh2;
   for (int i = 0; i < net->ncomps; i++) {
     //sprintf(str, "spline_%d.txt", i);
     //fh2 = fopen(str, "w");
     gsl_spline_init(spline_arr[i], time, c2D[i], lcount);
     //for (int k = 0; k < 101; k++) {
	//xi = (1-k/100.0)*time[0] + (k/100.0)*time[lcount-1];
	//yi = gsl_spline_eval(spline_arr[i], xi, acc);
	//fprintf(fh2, "%g\t%g\n", xi, yi);
     //}
     //fclose(fh2);
   }
#endif
}

void update_bc2(network *net, double time) {
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
  // Override to have influx at node 0
  /*n_k = net->knot[0]; 
  n_k->P = 0.; 
  Q = 289.*(1. - 1./cosh(a*time))*sin(a*time); 
  //Q = 289.;
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
  Q = -289.*(1. - 1./cosh(a*time))*cos(a*time); 
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
  //err_msg("Complete");*/
}













