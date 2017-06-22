#include "header.h"

static char text[80];
static gpipe_ptr pipe;
static gcomp_ptr cmp;
static FTYPE **wpi, **wmi;
static FTYPE **wpf, **wmf;
static FTYPE boost;
static long int num;

void init_forward_interior() {
	wpi = malloc(net.pipes*sizeof(FTYPE *));
	wpf = malloc(net.pipes*sizeof(FTYPE *));
	wmi = malloc(net.pipes*sizeof(FTYPE *));
	wmf = malloc(net.pipes*sizeof(FTYPE *));
	for (int j = 0; j < net.pipes; j++) {
		wpi[j] = malloc(net.pipe[j].N*sizeof(FTYPE));
		wpf[j] = malloc(net.pipe[j].N*sizeof(FTYPE));
		wmi[j] = malloc(net.pipe[j].N*sizeof(FTYPE));
		wmf[j] = malloc(net.pipe[j].N*sizeof(FTYPE));
	}
}

void free_forward_interior() {
	for (int j = net.pipes-1; j > -1; j--) {
		free(wpi[j]); 
		free(wpf[j]); 
		free(wmi[j]); 
		free(wmf[j]); 
	}
	free(wpi);
	free(wpf);
	free(wmi);
	free(wmf);
}

void forward_interior() {
	dmesg("forward_interior:\n", 0);
	for (int j = 0; j < net.pipes; j++) {
		printf("Forward interior of pipe %d\n", j);
		forward_pipe(j);
	}
	dmesg("forward_interior:\tPassed\n", 0);
}

void forward_pipe(int j) {
	num = net.pipe[j].N;
	pipe = &net.pipe[j];
	cmp = pipe->cmp;
	for (int k = 0; k < num; k++) {
		wpi[j][k] = pipe->y[0][k] + pipe->y[1][k];
		wmi[j][k] = pipe->y[0][k] - pipe->y[1][k];
	}
	memmove(&wpf[j][1], &wpi[j][0], (num-1)*sizeof(FTYPE));
	memmove(&wmf[j][0], &wmi[j][1], (num-1)*sizeof(FTYPE));
	if (cmp) boost = cmp->boost[par.curr];
	else boost = 1.0;
	printf("Taking temporal data from step %d\n", par.curr);
	wpf[j][0] = boost*(pipe->left->p[0]) + pipe->fs[0];  
	wmf[j][num-1] = (pipe->right->p[0]) - pipe->fd[0];
	for (int k = 0; k < num; k++) {
		pipe->y[0][k] = 0.5*(wpf[j][k] + wmf[j][k]);
		pipe->y[1][k] = 0.5*(wpf[j][k] - wmf[j][k]);
	}
	if (DEBUG_MODE) {
		sprintf(text, "debug/chars_%03d.log", j);
		FILE *fh = fopen(text, "w");
		fprintf(fh, "# 1. slot 2. w+ 3. w-\n\n");
		for (int k = 0; k < num; k++) {
			fprintf(fh, "%4d\t%.12e\t%.12e\n", k, wpf[j][k], wmf[j][k]);
		}
		fclose(fh);
	}
}

