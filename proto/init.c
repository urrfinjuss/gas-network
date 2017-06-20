#include "header.h"

static char text[80];

void call_init_network(char *filename) {
	read_input_file(filename);
	read_network_list(par.name);
	call_init_nodes();
	call_init_comps();
	call_init_pipes();
}

void call_init_nodes() { 
	dmesg("call_init_nodes:\n", 0);
	for (long int j = 0; j < net.nodes; j++) {
		net.node[j].var = malloc(par.nsteps*sizeof(FTYPE));  
		sprintf(text, "initialize/bc/nodes/node_%03ld.txt", j);
		long int nbase = count_data(text);
		FTYPE *Tbase = malloc((nbase+1)*sizeof(FTYPE));
		FTYPE *Vbase = malloc((nbase+1)*sizeof(FTYPE));
		FILE *fh = fopen_set(text, "r", 3);
		char line[80];
		for (long int k = 0; k < nbase+1; k++) {
			fgets(line, 80, fh);
			sscanf(line, "%lf\t%lf", &Tbase[k], &Vbase[k]);
			if (net.node[j].type) Vbase[k] = 1.0e+03*Vbase[k];
		}
		fclose(fh);
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(par.iType, nbase+1);
		gsl_spline_init(spline, Tbase, Vbase, nbase+1);
		for (long int n = 0; n < par.nsteps; n++) {
		  net.node[j].var[n] = gsl_spline_eval(spline, n*par.dt, acc);
		}
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
		if (DEBUG_MODE) verify_interp_nodes(j, nbase, Tbase, Vbase);
		free(Tbase);  
		free(Vbase);
	}
	dmesg("call_init_nodes:\tPassed\n",0);
}

void call_init_comps() { 
	dmesg("call_init_comps:\n", 0);
	for (long int j = 0; j < net.comps; j++) {
		net.comp[j].boost = malloc(par.nsteps*sizeof(FTYPE));  
		sprintf(text, "initialize/bc/comps/comp_%03ld.txt", j);
		long int nbase = count_data(text);
		FTYPE *Tbase = malloc((nbase+1)*sizeof(FTYPE));
		FTYPE *Vbase = malloc((nbase+1)*sizeof(FTYPE));
		FILE *fh = fopen_set(text, "r", 3);
		char line[80]; 
		for (long int k = 0; k < nbase+1; k++) {
			fgets(line, 80, fh);
			sscanf(line, "%lf\t%lf", &Tbase[k], &Vbase[k]);
		}
		fclose(fh); 
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(par.iType, nbase+1);
		gsl_spline_init(spline, Tbase, Vbase, nbase+1);
		for (long int n = 0; n < par.nsteps; n++) {
		  net.comp[j].boost[n] = gsl_spline_eval(spline, n*par.dt, acc);
		}
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
		if (DEBUG_MODE) verify_interp_comps(j, nbase, Tbase, Vbase); 
		free(Tbase);  
		free(Vbase);
	}
	dmesg("call_init_comps:\tPassed\n",0);
}

void call_init_pipes() {
	dmesg("call_init_pipes:\n", 0);
	for (long int j = 0; j < net.pipes; j++) {
		gpipe_ptr pipe = &net.pipe[j];
		sprintf(text, "initialize/ic/data_%03ld.txt", j);
		pipe->N = count_data(text)-1;
		pipe->y[0] = malloc(pipe->N*sizeof(FTYPE));
		pipe->y[1] = malloc(pipe->N*sizeof(FTYPE));
		FILE *fh = fopen_set(text, "r", 3);
		char line[80];
		FTYPE tmp;
		fgets(line, 80, fh);
		sscanf(line, "%lf\t%lf\t%lf", &tmp, &(pipe->left->p[0]), &pipe->fs[0]);
		(pipe->left)->p[0] =   1.0e+03*(pipe->left)->p[0];
		pipe->fs[0] = par.sound*pipe->fs[0];
		for (long int k = 0; k < pipe->N; k++) {
			fgets(line, 80, fh);
			sscanf(line, "%lf\t%lf\t%lf", &tmp, &(pipe->y[0][k]), &(pipe->y[1][k]));
			pipe->y[0][k] =   1.0e+03*pipe->y[0][k];
			pipe->y[1][k] = par.sound*pipe->y[1][k];
			if (k == 0) pipe->dx = 1.0e+03*tmp;
		}
		fgets(line, 80, fh);
		sscanf(line, "%lf\t%lf\t%lf", &tmp, &(pipe->right->p[0]), &(pipe->fd[0]));
		pipe->right->p[0] =   1.0e+03*(pipe->right->p[0]);
		pipe->fd[0] = par.sound*(pipe->fd[0]);
		fclose(fh);
		if (DEBUG_MODE) verify_set_pipes(j);
		if (pipe->dx != par.dx) {
			printf("Mismatch in pipe discretization.\n");
			exit(0);
		}
	}
	dmesg("call_init_pipes:\tPassed\n", 0);
}

