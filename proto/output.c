#include "header.h"


void network_snapshot() {
  char line[80];
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(line, "out/data_%03ld.txt", j);
		dump_pipe(&net.pipe[j], line);
	}
}

void dump_pipe(gpipe_ptr in, char *line) {
		FILE *fh = fopen(line, "w");
		fprintf(fh, "# 1. x, km 2. Pressure, MPa 3. Flux, kg/m^2/s\n");
		fprintf(fh, "# N = %ld\n\n", in->N);
		for (long int j = 0; j < in->N; j++) {
			fprintf(fh,"%.12e\t%.12e\t%.12e\n",j*in->dx,1.0e-6*in->y[0][j],in->y[1][j]/par.sound);
		}
		exit(1);
		fclose(fh);
}

