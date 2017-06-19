#include "header.h"

static char text[80];

void call_init_nodes() { 
	dmesg("call_init_nodes:\n", 0);
	const gsl_interp_type *iType;
	printf("Interpolation type selected: %s\n", par.intmethod);
	if (strcmp(par.intmethod,"linear")==0)						iType = gsl_interp_linear;
	else if (strcmp(par.intmethod,"cspline")==0)						iType = gsl_interp_cspline;
	else if (strcmp(par.intmethod,"akima")==0)							iType = gsl_interp_akima;
	else if (strcmp(par.intmethod,"steffen")==0)						iType = gsl_interp_steffen;
	else if (strcmp(par.intmethod,"cspline_periodic")==0)	iType = gsl_interp_cspline_periodic;
	else if (strcmp(par.intmethod,"akima_periodic")==0) 		iType = gsl_interp_akima_periodic;
	else dmesg("Unknown Interpolation Type.\n",1);
	//  --- Setting BC at nodes
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
		}
		fclose(fh);
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(iType, nbase+1);
		gsl_spline_init(spline, Tbase, Vbase, nbase+1);
		for (long int n = 0; n < par.nsteps; n++) {
		  net.node[j].var[n] = gsl_spline_eval(spline, n*par.dt, acc);
		}
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
		if (DEBUG_MODE) {
			sprintf(text, "debug/nodes/base_%03ld.txt", j);
			fh = fopen(text, "w");
			fprintf(fh, "# 1. base T 2. base V\n\n");
			for (long int n = 0; n < nbase + 1; n++) {
				fprintf(fh, "%.12e\t%.12e\n", Tbase[n], Vbase[n]);
			}
			fclose(fh);
			sprintf(text, "debug/nodes/%s_%03ld.txt", par.intmethod, j);
			fh = fopen(text, "w");
			fprintf(fh, "# 1. T 2. V (cspline)\n\n");
			for (long int n = 0; n < par.nsteps; n++) {
				fprintf(fh, "%.12e\t%.12e\n", n*par.dt, net.node[j].var[n]);
			}
			free(Tbase);  
			free(Vbase);
			fclose(fh);
		}
	}
	dmesg("call_init_nodes:\tComplete\n",0);
}
