#include "network.h"

int main (int argc, char *argv[])
{
  if (argc != 2) {
    printf("Usage:\n\t%s input.cfg\n", argv[0]);
    exit(1);
  }
  network syst;  
  init_network(argv[1], &syst);
  init_links(&syst);
  err_msg("here");

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  int i,n;
  double gauss,gamma;  

  n=atoi(argv[1]);
  for (i=0;i<n;i++)
    {
      gauss=gsl_ran_gaussian(r,2.0);
      gamma=gsl_ran_gamma(r,2.0,3.0);
      printf("%2.4f %2.4f\n", gauss,gamma);
    }
  return(0);
}

int noiseg(pipe_ptr in) {

return 0;
}
