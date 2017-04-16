#include "network.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

static const gsl_rng_type *T;
static gsl_rng *rand_gen;
static unsigned long int seed = 0; // must be non-zero for best results
static unsigned long int N, Nt;
static double *rarr, *new;
static double tau, delta, A, lcorr;
static double a1, a2, overN, L, cs;
static noise_ptr nse;
//static fftw_complex *ft; 	// stores DFT of the noise in time
static fftw_complex *fx; 	// stores DFT of the noise correlated in space 
static fftw_complex *gss;	// stores Gaussian kernel
static fftw_complex *gsk;	// stores DFT of the Gaussian kernel
//static fftw_complex *full; 	// stores DFT of the noise in time with zero padding
static fftw_plan p1, p2;

void init_noise(network_ptr net) {

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rand_gen = gsl_rng_alloc(T); 

  struct timeval tv;
  FILE *devrandom;

 //if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
 /*} else {
   fread(&seed,sizeof(seed),1,devrandom);
   fclose(devrandom);
 }*/


  gsl_rng_set(rand_gen, seed);

  printf ("generator type: %s\n", gsl_rng_name (rand_gen));
  printf ("seed = %lu\n", seed);
  printf ("first value = %lu\n", gsl_rng_get (rand_gen));


  nse = net->nse;
  tau = 60.*(nse->tau);
  lcorr = nse->lc; 
  //delta = 0.1*(nse->tau);  // must relate this to time step. will pass as argument
  //double dx = 1000;
  //double dt = dx/(net->c);
  delta = 1000./(net->npcent)/(net->c);
  L = (net->link[0]->L); 
  N = (net->link[0]->N); 
  cs = net->c;
  overN = 1./pow(N,0.5);
  Nt = round((net->tmax)/delta);
  printf("Noise Parameters:\nCorrelation length\t%.2f km\nCorrelation time\t%.2f mins\nNumber of subintervals\t%ld\n", lcorr, tau/60., N+1);
  printf("Timesteps\t\t%ld\n", Nt);
  A = nse->A;  
  a1 = exp(-delta/tau); a2 = sqrt(2*A*delta/tau);

  rarr = malloc(N*sizeof(double));
  new = malloc(N*sizeof(double));
  gss = fftw_malloc(N*sizeof(fftw_complex));
  gsk = fftw_malloc(N*sizeof(fftw_complex)); 		// stores Fourier of Gauss kernel
  
  fx = net->link[0]->fx;
  //fftw_malloc(N*sizeof(fftw_complex));



  //memset(new, 0, N*sizeof(double));
  //memset(rl, 0, N*Nt*sizeof(double));
  //memset(full, 0, N*Nt*sizeof(fftw_complex));
  
  memset(gsk, 0, N*sizeof(fftw_complex));

  for (int j = 0; j < N; j++) {
	new[j] = sqrt(A)*gsl_ran_gaussian_ziggurat(rand_gen, 1.);		// initial data for the noise processes in time
	gsk[j] = exp(-0.5*pow(2.*pi*j*lcorr/sqrt(2.)/L, 2))*sqrt(pi*sqrt(pi))*overN*sqrt(lcorr/L)*0.912;  			// L/N to have true convolution
        if (j > N/2) gsk[j] = exp(-0.5*pow(2.*pi*(j-N)*lcorr/sqrt(2.)/L, 2))*sqrt(pi*sqrt(pi))*overN*sqrt(lcorr/L)*0.912; 	// L/N

        //if (j == 1) gsk[j] = L/N;
	
  }
  p1 = fftw_plan_dft_1d(N, gsk, gss, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p1);
  FILE *fh = fopen("gauss.txt","w");
  fprintf(fh, "# 1. x 2. g(x)\n\n");
  for (int j = 0; j < N; j++) fprintf(fh, "%.12e\t%.12e\t%.12e\n", L*j/N, creal(gss[j]), cimag(gss[j]));
  fclose(fh);

  printf("G(k=0) = %.16e\n", creal(gsk[0])/overN);
  fftw_destroy_plan(p1);
  p1 = fftw_plan_dft_1d(N, fx, fx, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_1d(N, fx, fx, FFTW_BACKWARD, FFTW_ESTIMATE);

  gsl_rng_set (rand_gen, seed);

  for (int j = 0; j < N; j++) fx[j] = cexp(2.I*pi*j/N);	// proper coefficient is L/N
  fftw_execute(p1);
  for (int j = 0; j < N; j++) fx[j] = fx[j]*gsk[j];
  fftw_execute(p2);

  fh = fopen("convolution.txt","w");
  fprintf(fh, "# 1. x 2. g(x)\n\n");
  for (int j = 0; j < N; j++) fprintf(fh, "%.12e\t%.12e\t%.12e\n", L*j/N, creal(fx[j]), cimag(fx[j]));
  fclose(fh);



  /*int rank = 1;
  int n[] = {Nt};
  int howmany = N;
  int idist = 1, odist = 1;
  int istride = N, ostride = N;
  int *inembed = n, *onembed = n;

  p1 = fftw_plan_many_dft(rank, n, howmany,
                                  full, inembed,
                                  istride, idist,
                                  ft, onembed,
                                  ostride, odist,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

  p2 = fftw_plan_many_dft(rank, n, howmany,
                                  ft, onembed,
                                  ostride, odist,
				  full, inembed,
                                  istride, idist,
                                  FFTW_FORWARD, FFTW_ESTIMATE); */



}

void generate_noise_step() {
  for (int j = 0; j < N; j++) {
    rarr[j] = gsl_ran_gaussian_ziggurat(rand_gen, 1.); // return gauss random number with unit variance  
    new[j] = a1*new[j] + a2*rarr[j];
    fx[j] = new[j];
  }
  fftw_execute(p1);
  for (int j = 0; j < N; j++) {
     fx[j] = fx[j]*gsk[j];	// Parseval's Thm, probably is missing a bunch of coeff
  }
  fftw_execute(p2);
}

void generate_poor_mans_noise(double time) {
  for (int j = 0; j < N; j++) {
     //fx[j] = -A*(0.1 - cos(2.*pi*(j*L/N)/lcorr)*cos(2.*pi*time/tau));
     //fx[j] = A*(1. - cos(2.*pi*(time/tau)));
     //fx[j] = A*(1. - exp(-2.*pi*pow(time/tau,2)));
     fx[j] = A*sin(2.*pi*(time/tau)); 
     //fx[j] = 1.;
  }
}


void play_with_noise(network_ptr net) {
  init_noise(net);
  
  FILE *fh = fopen("noise.txt","w");
  fprintf(fh, "# 1. t 2. xi[0...N-1]\n\n");
  for (int k = 0; k < Nt; k++) {
    generate_noise_step();  
    fprintf(fh, "%.15e", k*delta);
      for (int n = 0; n < N; n=n+1600) {
        fprintf(fh, "\t%.15e", creal(fx[n]));
      }
    fprintf(fh, "\n");
  }
  fclose(fh);

  fh = fopen("spcorr.txt","w");
  fprintf(fh, "# 1. x 2. zeta(x,Tmax)\n\n");
  for(int j = 0; j < N; j++) fprintf(fh, "%.15e\t%.15e\n",L*j/N, creal(fx[j]));
  fclose(fh);

  //fh = fopen("convolution.txt","w");
  //fprintf(fh, "# 1. x 2. g(x)\n\n");
  //for (int j = 0; j < N; j++) fprintf(fh, "%.12e\t%.12e\t%.12e\n", L*j/N, creal(fx[j]), cimag(fx[j]));
  //fclose(fh);

  double *mean = malloc(N*sizeof(double));
  double *SD = malloc(N*sizeof(double));
  double *aut = malloc(N*Nt*sizeof(double));
  /*printf("\nGenerated Noise Params:\nMean/SD/Autocorrelation:");
  for (int k = 0; k < N; k++) {
    mean[k] = gsl_stats_mean(&rl[k], N, Nt);
    SD[k] = gsl_stats_sd(&rl[k], N, Nt);
    printf("\t%.2f/%.2f/%.2f  ", mean[k], SD[k], aut[k]);    
  }*/
  printf("\nFinal Time\t%.2f mins\n", Nt*delta/60.);

}



int noise_status(noise_p *noise) {
  char msg[80];
  debug_msg("\nNoise parameters:\n\n");
  sprintf(msg, "Tau = %.8e\nL_c = %.8e\nAmp = %.8e\n", noise->tau, noise->lc, noise->A);
  debug_msg(msg);
  return 0;
}





