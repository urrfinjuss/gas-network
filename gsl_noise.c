#include "network.h"

static network syst;

int init_network(char *filename, network *net) {
  FILE* fh = fopen(filename,"r");
  char line[80], param[80], value[80];
  noise_p *nse = net->noise;
  nse = malloc(sizeof(noise_p));

  if (fh == NULL) err_msg("Cannot open configuration file.");
  else while (fgets(line, 80, fh) != NULL) {
    sscanf(line, "%s\t%s", param, value);
    if (strcmp(param,"#name=") == 0) sprintf(net->fname,"%s", value);
    if (strcmp(param,"#nodes=") == 0) net->nnodes = atoi(value);
    if (strcmp(param,"#links=") == 0) net->nlinks = atoi(value);
    if (strcmp(param,"#corr_time=") == 0) nse->tc = atof(value);
    if (strcmp(param,"#corr_dist=") == 0) nse->lc = atof(value);
    if (strcmp(param,"#amplitude=") == 0) nse->A = atof(value);
  }
  fclose(fh);
  noise_status(nse);
  err_msg("Success");
  return 1;
}


int main (int argc, char *argv[])
{
  if (argc != 2) {
    printf("Usage:\n\t%s input.cfg\n", argv[0]);
    exit(1);
  }
  network syst;  
  init_network(argv[1], &syst);
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
