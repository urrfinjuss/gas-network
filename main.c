#include "network.h"


int main (int argc, char *argv[])
{
  char cwd[1024], outdir[1024];
  char msg[1024], title[80];
  network syst;  

  if (argc != 2) {
    printf("Usage:\n\t%s input.cfg\n", argv[0]);
    exit(1);
  }
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    printf("Current working dir: %s\n", cwd);
    sprintf(syst.current_dir, "%s", cwd);
  }
  else err_msg("getcwd() error");
  
  init_network(argv[1], &syst);
  init_links(&syst);
  init_data(&syst);
  init_evolve(&syst);
  mgl_init_draw(&syst);
  
  for (int n = 0; n < syst.nlinks; n++){ 
    sprintf(outdir, "%s/figures_%03d", cwd, n);
    debug_msg("Saving gas pipe data to:");
    debug_msg(syst.current_dir);
    debug_msg("\n");
    printf("Saving gas pipe data to:%s\n", outdir);
    mkdir(outdir, 0777);
    sprintf(msg, "%s/%s_%03d.png", outdir, syst.dname, 0);
    sprintf(title, "t = %f", 0.);
    mgl_draw_pipe(syst.link[n], &syst, msg, title);
  }
  /*mgl_draw_pressure(&syst, "pressure.png", "");
  mgl_draw_flux(&syst, "flux.png", "");*/


  evolve_network(&syst);

  /*mgl_draw_pressure(&syst, "pressure2.png", "");
  mgl_draw_flux(&syst, "flux2.png", "");*/

  err_msg("Complete");

  /*gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  int i,n;
  double gauss,gamma;  

  n=atoi(argv[1]);
  for (i=0;i<n;i++)
    {
      gauss=gsl_ran_gaussian(r,2.0);
      gamma=gsl_ran_gamma(r,2.0,3.0);
      printf("%2.4f %2.4f\n", gauss,gamma);
    }
  */
  //return(0);
}


