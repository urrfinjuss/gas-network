#include "header.h"

static char text[80];
static FILE *fstatus;

void hyperbolic_step() {
	dmesg("hyperbolic_step:\n", 0);
	forward_nodes();
	forward_interior();
	sync_nodes();
	dmesg("hyperbolic_step:\tPassed\n", 0);
}

void evolve_network() {
	dmesg("evolve_network:\n", 0);

	long int cntr = 0;
	FTYPE progress = 0.0;
	fstatus = fopen("progress.log","w");
	fprintf(fstatus, "# 1. Progress Bar\n\n");
	fprintf(fstatus, "%.2f\n", progress);
	fclose(fstatus);

	for (long int j = 0; j < par.nsteps; j++) {
		if ( j % par.skip == 0) {
			network_snapshot(cntr);
			cntr++;
		}
		hyperbolic_step();
		progress = 100.*(j+1)/par.nsteps;
		if ( (j+1) % (int)(par.nsteps/100) ) {
			fstatus = fopen("progress.log","a");
			fprintf(fstatus, "%.2f\n", progress);
			fclose(fstatus);
		}
	}
	network_snapshot(cntr);
	fstatus = fopen("progress.log","a");
	fprintf(fstatus, "%.2f\n", 100.0);
	fclose(fstatus);
	dmesg("evolve_network:\tPassed\n", 0);
}

