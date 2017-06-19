#include "header.h"

static char text[80];

int read_pipe(char *fname, gpipe_ptr in) {
	FILE *fh = fopen(fname, "r");
	if (fh) {
		char line[256];
		FTYPE x;
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		long int n = 0;
		dmesg(line, 0);
		while ( fgets(line, 256, fh) != NULL ) n++;
		sprintf(text, "Found %ld points\n", n-1);
		dmesg(text, 0);
		fclose(fh);
		in->N = n;
		in->dx = in->len/(n-1);
		in->y[0] = malloc(n*sizeof(FTYPE));
		in->y[1] = malloc(n*sizeof(FTYPE));
		fh = fopen(fname, "r");
		fgets(line, 256, fh); printf("%s", line);
		fgets(line, 256, fh); printf("%s", line);
		fgets(line, 256, fh); printf("%s", line);
		long int j = 0;
		while ( fgets(line, 256, fh) != NULL) {
			sscanf(line, "%lf\t%lf\t%lf", &x, &in->y[0][j], &in->y[1][j]);
			in->y[0][j] = 1.0e+6*(in->y[0][j]);
			in->y[1][j] = par.sound*(in->y[1][j]);
			j++;
			//printf("%.8e\t%.8e\n", in->y[0][j], in->y[1][j]);
		}
		fclose(fh);
		sprintf(text, "Disc Step %f m (L = %.1f m N = %ld)\n", in->dx, in->len, n-1);
		dmesg(text, 0);
		par.dx = in->dx;
		par.dt = in->dx/par.sound;
		return 0;
	} else return 1;
}

int load_initial_data() {
	char	fname[80];
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(fname, "./initialize/ic/data_%03ld.txt", j);
		FILE *fh = fopen(fname, "r");
		if (fh) {
			sprintf(text, "Found data for pipe %ld:", j);
			dmesg(text, 0);
			fclose(fh);
		} else dmesg("Missing data for one of the pipes\n", 1);
		if (!read_pipe(fname, &net.pipe[j])) {
			printf("Setting the pipe %ld at %p:\t", j, &net.pipe[j]);
			dmesg("Loaded\n", 0);
		}
	}
	return 0;
}
