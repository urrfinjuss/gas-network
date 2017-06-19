#include "header.h"


int call_init_temporal() {
	dmesg("call_init_temporal:\n", 0);
	long int Nf = par.nsteps;
	par.nsteps = floor(par.tmax/par.dt);
	par.tmax = par.nsteps*par.dt;
	if (DEBUG_MODE) {
		printf("Simulation time is %.10f hours\n", par.tmax/3600.); 
		printf("The simulation time step is %.16f secs\n", par.dt); 
		printf("Total time steps to be made %u\n", par.nsteps);
	}
	dmesg("call_init_temporal: Complete\n", 0);
	return 0;
}
