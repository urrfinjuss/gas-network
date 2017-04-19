#include "network.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

static gpipe_ptr p_k;
static node_ptr n_k;
static double ovsq2, sq2, cs;
static double *time;
static double **c2D, **d2D;
static char msg[80];
static char dmsg[512];
static gsl_interp_accel *acc;
static gsl_spline **spline_arr, **spline_arr2;


void init_compressors(network *net){
   FILE *fh = fopen("c.txt","r");
   size_t count = 16384;
   char *line = malloc(16384);
   char msg[8192];
   int rnum, m = 0;
   int lcount = 0;

   sq2 = sqrt(2.);
   ovsq2 = 1./sq2;
   cs = net->c;
   fgets(msg, 8192, fh);
   fgets(msg, 8192, fh);
   while(1){
	//printf("%d\n", lcount);
	if(fgets(msg, 8192, fh)==NULL) break;
	lcount++;
   }
   sprintf(msg, "\nTemporal data available for %8f hours\n", 1.*(lcount-1)/1800.);
   debug_msg(msg);
   rewind(fh);
   fgets(msg, 8192, fh);
   fgets(msg, 8192, fh);

   double *link_matrix = malloc(sizeof(double)*(net->ncomps + 1)*lcount);
   //err_msg("Complete");

   fgets(line, (int) count, fh);
   //getline(&line, &count, fh);
   
   int read = -1, cur = 0, cCount = 0;
   while( sscanf(line+cur, "%lf%n", &link_matrix[cCount], &read) == 1) {
	cur+=read;
	cCount++;
   }
   int rCount = 1;
   //while(getline(&line, &count, fh)!=-1) {
   while(fgets(line, (int) count, fh)!=NULL) {
	rCount++;
   }
   rewind(fh);
   sprintf(dmsg,"%d\n%d\n%d\n", rCount, cCount, net->ncomps);
   debug_msg(dmsg);
   int i = 0;
   //while(getline(&line, &count, fh)!=-1) {
   while(fgets(line, (int) count, fh)!=NULL) {
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
   for (int j = 0; j < lcount; j++) time[j] = link_matrix[(cCount)*j];
   printf("time = %f\n", time[1]);
   for (int i = 0; i < net->ncomps; i++) {
	for (int j = 0; j < lcount; j++) {
	   c2D[i][j] = link_matrix[cCount*j+i+1];
	}
   }
   acc = gsl_interp_accel_alloc();
   const gsl_interp_type *t = gsl_interp_cspline_periodic;
   spline_arr = malloc(net->ncomps*sizeof(gsl_spline *));
   for (int i = 0; i < net->ncomps; i++) {
	spline_arr[i] = gsl_spline_alloc(t, lcount);
	gsl_spline_init(spline_arr[i], time, c2D[i], lcount);
   }
#if 0
   double xi, yi;
   char str[24];
   int Nx = 200;
   FILE *fh2;
   for (int i = 0; i < net->ncomps; i++) {
     sprintf(str, "cspline_%d.txt", i);
     fh2 = fopen(str, "w");
     gsl_spline_init(spline_arr[i], time, c2D[i], lcount);
     for (int k = 0; k < Nx+1; k++) {
	xi = (1-1.*k/Nx)*time[0] + (1.*k/Nx)*time[lcount-1];
	yi = gsl_spline_eval(spline_arr[i], xi, acc);
	fprintf(fh2, "%g\t%g\n", xi, yi);
     }
     fclose(fh2);
   }

#endif
   update_compressors(net, 0.);
}

void init_demands(network *net){
   FILE *fh = fopen("d.txt","r");
   size_t count = 16384;
   char *line = malloc(16384);
   char msg[16384];
   int rnum, m = 0;
   int lcount = 0;

   fgets(msg, 16384, fh);
   fgets(msg, 16384, fh);
   while(1){
	//printf("%d\n", lcount);
	if(fgets(msg, 16384, fh)==NULL) break;
	lcount++;
   }
   sprintf(msg, "\nTemporal data available for %8f hours\n", 1.*(lcount-1)/1800.);
   debug_msg(msg);
   rewind(fh);
   fgets(msg, 16384, fh);
   fgets(msg, 16384, fh);

   double *link_matrix = malloc(sizeof(double)*(net->nnodes + 1)*lcount);
   //err_msg("Complete");

   //getline(&line, &count, fh);
   fgets(line, (int) count, fh);
   int read = -1, cur = 0, cCount = 0;
   while( sscanf(line+cur, "%lf%n", &link_matrix[cCount], &read) == 1) {
	cur+=read;
	cCount++;
   }
   int rCount = 1;
   while(fgets(line, (int) count, fh)!=NULL) {
   //while(getline(&line, &count, fh)!=-1) {
	rCount++;
   }
   rewind(fh);
   sprintf(dmsg,"%d\n%d\n%d\n", rCount, cCount, net->nnodes);
   debug_msg(dmsg);
   int i = 0;
   //while(getline(&line, &count, fh)!=-1) {
   while(fgets(line, (int) count, fh)!=NULL) {
     read = -1;
     cur  =  0;
     while(sscanf(line+cur, "%lf%n", &link_matrix[i], &read) == 1) {
	cur+=read;
	i=i+1;
     }
   }
   fclose(fh);
   sprintf(msg, "Demand data has %d entries\n", cCount*rCount);
   debug_msg(msg);

   time = malloc(lcount*sizeof(double));
   d2D = malloc(net->nnodes*sizeof(double *));
   for (int i = 0; i < net->nnodes; i++) d2D[i] = malloc(lcount*sizeof(double));
   for (int j = 0; j < lcount; j++) time[j] = link_matrix[(net->nnodes+1)*j];
   for (int i = 0; i < net->nnodes; i++) {
	for (int j = 0; j < lcount; j++) {
	   d2D[i][j] = link_matrix[(net->nnodes+1)*j+i+1];
	}
   }

   acc = gsl_interp_accel_alloc();
   const gsl_interp_type *t = gsl_interp_cspline_periodic;
   spline_arr2 = malloc(net->nnodes*sizeof(gsl_spline *));
   for (int i = 0; i < net->nnodes; i++) {
	spline_arr2[i] = gsl_spline_alloc(t, lcount);
	gsl_spline_init(spline_arr2[i], time, d2D[i], lcount);
   }
#if 0

   double xi, yi;
   char str[24];
   FILE *fh2;
   for (int i = 0; i < net->nnodes; i++) {
     sprintf(str, "dspline_%d.txt", i);
     fh2 = fopen(str, "w");
     gsl_spline_init(spline_arr2[i], time, d2D[i], lcount);
     for (int k = 0; k < 401; k++) {
	xi = (1-k/400.0)*time[0] + (k/400.0)*time[lcount-1];
	yi = gsl_spline_eval(spline_arr2[i], xi, acc);
	fprintf(fh2, "%g\t%g\n", xi, yi);
     }
     fclose(fh2);
   }
   exit(1);
#endif
   update_demands(net, 0.);
}




void update_compressors(network* net, double ctime){

   for (int j = 0; j < net->ncomps; j++) {
     (net->cssr[j]).cratio_prev = (net->cssr[j]).cratio;
     (net->cssr[j]).cratio = gsl_spline_eval(spline_arr[j], ctime, acc);
     sprintf(msg, "C%d (%p): Compression %.4f\n", j, &(net->cssr[j]), (net->cssr[j]).cratio);
     debug_msg(msg);
   }
}

void update_demands(network *net, double ctime) {

   for (int j = 0; j < net->nnodes; j++) {
     net->knot[j]->D_prev = net->knot[j]->D;
     net->knot[j]->D = gsl_spline_eval(spline_arr2[j], ctime, acc);
     sprintf(msg, "D%d: Demand %f\n", j, (net->knot[j])->D);
     debug_msg(msg);
   }
}









