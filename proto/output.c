#include "header.h"


void network_snapshot() {
  char line[80];
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(line, "out/data_%03ld.txt", j);
		dump_pipe(&net.pipe[j], line);
	}
	dmesg("network_snapshot:\tPassed\n",0);
}

void dump_pipe(gpipe_ptr in, char *line) {
		FILE	*fh = fopen(line, "w");
		FTYPE	P,F;
		fprintf(fh, "# 1. x, km 2. Pressure, MPa 3. Flux, kg/m^2/s\n");
		fprintf(fh, "# N = %ld\n\n", in->N);
		for (long int j = 0; j < in->N; j++) {
			P = 1.0e-6*(in->y[0][j]);
			F = (in->y[1][j])/par.sound;
			fprintf(fh, "%.12e\t%.12e\t%.12e\n", 1.0e-3*j*in->dx, P, F);
		}
		fclose(fh);
}

