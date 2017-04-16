#include "network.h"


int main (int argc, char *argv[])
{
  char cwd[1024], outdir[1024], outdir2[1024];
  char msg[1024], msg2[1024], title[80];
  network syst; 

#if BALANCECHECK
  network systb;
#endif 

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
  // naive 			to here
  init_compressors(&syst);  

  init_data(&syst);
  init_evolve(&syst);
  //init_compressors(&syst);   moved from here
  init_demands(&syst);
#if NOISE
  init_noise(&syst);
#endif
  //mgl_init_draw(&syst);

#if BALANCECHECK
  init_network(argv[1], &systb);
  init_links(&systb);
  init_data(&systb);
  init_evolve(&systb);
  init_compressors(&systb);  
  init_demands(&systb);
  //mgl_init_draw(&systb);
#endif

  prepare_arrays_adiabatic(&syst); // initialize adiabatic code arrays
  
  for (int n = 0; n < syst.nlinks; n++){ 
    sprintf(outdir, "%s/pipe_%03d", cwd, n); 
    mkdir(outdir, 0777);
    debug_msg("Saving gas pipe data to:");
    debug_msg(syst.current_dir);
    debug_msg("\n");

    //sprintf(outdir, "%s/figures_%03d", cwd, n);
    //printf("Saving gas pipe data to:%s\n", outdir);
    //mkdir(outdir, 0777);
    //sprintf(msg, "%s/%s_%03d.png", outdir, syst.dname, 0);
    //sprintf(title, "t = %f", 0.);
    //mgl_draw_pipe(syst.link[n], &syst, msg, title);
  }

  sprintf(outdir2, "%s/network", cwd);
  mkdir(outdir2, 0777);
  sprintf(msg2, "%s/%s_%03d.png", outdir2, syst.dname, 0);
  //mgl_draw_network(&syst, msg2, title);

  //play_with_noise(&syst);
  //exit(1);

#if BALANCECHECK
  //printf("Generating balance of terms data\n");
  //evolve_network_balance(&syst, &systb);
#else
  evolve_network(&syst);
#endif
  //sync_nodes(&syst);
  //hyperbolic_step(&syst);

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


