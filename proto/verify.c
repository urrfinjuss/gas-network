#include "header.h"

static char text[80];

void verify_input_conf() {
	dmesg("verify_input_conf:\n", 0);
	printf("\t[Input Configuration]\n");
	printf("\tname=\t%s\n", par.name);
	printf("\tstep=\t%.8e\n", par.dx);
	printf("\tskip=\t%d\n", par.skip);
	printf("\ttmax=\t%.8e\n", par.tmax/3600.);
	printf("\tsounds=\t%.8e\n", par.sound);
	printf("\tDWdiss=\t%.8e\n", par.friction);
	printf("\tinterp=\t%s\n", par.intmethod);
	printf("\t[Simulation Configuration]\n");
	printf("\tTime step\t%.6f secs\n", par.dt); 
	printf("\tOutput every\t%.2f secs\n", par.skip*par.dt);
	printf("\tTotal steps\t%u\n", par.nsteps);
	printf("\tGSL Interp:\t%s\n", par.intmethod);
	dmesg("verify_input_conf:\tPassed\n", 0);
}

void verify_network_conf() {
	dmesg("verify_network_conf:\n", 0);
	verify_node_conf();
	verify_pipe_conf();
	verify_comp_conf();
	dmesg("verify_network_conf:\tPassed\n", 0);
}

void verify_node_conf() {
	dmesg("\tverify_node_conf:\n", 0);
	printf("\t[Node Configuration]\n");
	for (int j = 0; j < net.nodes; j++) {
		printf("\tNode %d at %p\n", j, &net.node[j]);
		if (net.node[j].type == 1) printf("\tPressure Type\n");
		else printf("\tTransport Type\n");
		printf("\tCompressors:\t%d\n", net.node[j].ncomp);
		if (net.node[j].ncomp) {
			printf("\tCompressor List:");
			for (int k = 0; k < net.node[j].ncomp; k++) {
				printf("\t%p", net.node[j].comp[k]);
			}
			printf("\n");
		}
		printf("\tIncoming pipes:\t%d\n", net.node[j].nleft);
		if (net.node[j].nleft) {
			printf("\tIncoming Pipe List:");
			for (int k = 0; k < net.node[j].nleft; k++) {
				printf("\t%p", net.node[j].left[k]);
			}
			printf("\n");
		}
		printf("\tOutgoing pipes:\t%d\n", net.node[j].nright);
		if (net.node[j].nright) {
			printf("\tOutgoing Pipe List:");
			for (int k = 0; k < net.node[j].nright; k++) {
				printf("\t%p", net.node[j].right[k]);
			}
			printf("\n");
		}
		printf("\n");
	}
	dmesg("\tverify_node_conf:\tPassed\n", 0);
}

void verify_pipe_conf() {
	dmesg("\tverify_pipe_conf:\n", 0);
	printf("\t[Pipe Configuration]\n");
	for (int j = 0; j < net.pipes; j++) {
		printf("\tPipe %d at %p\n", j, &net.pipe[j]);
		printf("\tSource Node:\t\t%p\n", net.pipe[j].left);
		printf("\tDestination Node:\t%p\n", net.pipe[j].right);
		printf("\tCompressor:\t\t%p\n", net.pipe[j].cmp);
		printf("\tLength:\t\t\t%e\n", net.pipe[j].len*1e-3);
		printf("\tWidth:\t\t\t%e\n\n", net.pipe[j].wid);
	}
	dmesg("\tverify_pipe_conf:\tPassed\n", 0);
}

void verify_comp_conf() {
	dmesg("\tverify_comp_conf:\n", 0);
	printf("\t[Compressor Configuration]\n");
	for (int j = 0; j < net.comps; j++) {
		printf("\tCompressor %d at %p\n", j, &net.comp[j]);
		printf("\tSource Node:\t%p\n", net.comp[j].loc);
		printf("\tBoost Pipe:\t%p\n\n", net.comp[j].dest);
	}
	dmesg("\tverify_comp_conf:\tPassed\n", 0);
}

void verify_interp_nodes(long int j, long int nbase, FTYPE *Tbase, FTYPE *Vbase) {
	sprintf(text, "debug/nodes/base_%03ld.txt", j);
	FILE *fh = fopen(text, "w");
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
	fclose(fh);
}

void verify_interp_comps(long int j, long int nbase, FTYPE *Tbase, FTYPE *Vbase) {
	sprintf(text, "debug/comps/base_%03ld.txt", j);
	FILE *fh = fopen(text, "w");
	fprintf(fh, "# 1. base T 2. base V\n\n");
	for (long int n = 0; n < nbase + 1; n++) {
		fprintf(fh, "%.12e\t%.12e\n", Tbase[n], Vbase[n]);
	}
	fclose(fh);
	sprintf(text, "debug/comps/%s_%03ld.txt", par.intmethod, j);
	fh = fopen(text, "w");
	fprintf(fh, "# 1. T 2. interp\n\n");
	for (long int n = 0; n < par.nsteps; n++) {
		fprintf(fh, "%.12e\t%.12e\n", n*par.dt, net.comp[j].boost[n]);
	}
	fclose(fh);
}

void verify_set_pipes(long int j) {
	sprintf(text, "debug/pipes/data_%03ld.txt", j);
	FILE *fh = fopen(text, "w");
	gpipe_ptr pipe = &net.pipe[j];
	gnode_ptr left = net.pipe[j].left;
	gnode_ptr right = net.pipe[j].right;
	gcomp_ptr cmp = net.pipe[j].cmp;
	FTYPE boost = 1.0;
	if (cmp) boost = cmp->boost[0];
	fprintf(fh, "# 1. x, km 2. Pressure, MPa 2. Flux, kg/m^2/s\n\n");
	fprintf(fh, "%.12e\t%.12e\t%.12e\n", 0., 1.0e-03*left->p[0]*boost, pipe->fs[0]/par.sound);
	for (long int n = 0; n < pipe->N; n++) {
		fprintf(fh, "%.12e\t%.12e\t%.12e\n", 1.0e-03*par.dx*(n+1), 1.0e-03*pipe->y[0][n], pipe->y[1][n]/par.sound);
	}
	fprintf(fh, "%.12e\t%.12e\t%.12e\n", 1.0e-03*pipe->len, 1.0e-03*right->p[0], pipe->fd[0]/par.sound);
	fclose(fh);
}
