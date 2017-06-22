#include "header.h"

static char line[80], text[80];
static FILE *fh_node;
static FILE *fh_pipe;
static FTYPE P,F;

void create_directory(char *path) {
	printf("Create Directory: %s\n", path);
#ifdef __linux
	mkdir(path, 0777);
#endif
#ifdef _WIN32
	_mkdir(path);
#endif
}

void init_node_output() {
	dmesg("init_node_output:\n", 0);
	create_directory("./result/nodes");
	for (long int j = 0; j < net.nodes; j++) {
		sprintf(line, "./result/nodes/data_%03ld.txt", j);
		fh_node = fopen(line, "w");
		fprintf(fh_node, "# 1. time, sec 2. Pressure, MPa\n");
		fprintf(fh_node, "# Auto-generated by Gas Toolbox\n\n");
		fclose(fh_node);
	}
	dmesg("init_node_output:\tPassed\n", 0);
}

void init_pipe_output() {
	dmesg("init_pipe_output:\n", 0);
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(line, "./result/pipe_%03ld", j);
		create_directory(line);
	}
	dmesg("init_pipe_output:\tPassed\n", 0);
}

void network_snapshot(long int k) {
	dmesg("network_snapshot:\n", 0);
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(line, "./result/pipe_%03ld/data_%03ld.txt", j, k);
		dump_pipe(&net.pipe[j], line);
	}
	for (long int j = 0; j < net.nodes; j++) {
		sprintf(line, "./result/nodes/data_%03ld.txt", j);
		dump_node(&net.node[j], line);
	}
	dmesg("network_snapshot:\tPassed\n",0);
}

void dump_pipe(gpipe_ptr in, char *filename) {
	fh_pipe = fopen(filename, "w");
	fprintf(fh_pipe, "# 1. x, km 2. Pressure, MPa 3. Flux, kg/m^2/s\n");
	fprintf(fh_pipe, "# T = %8.2f sec\tAuto-generated by Gas Toolbox\n\n", par.curr*par.dt);
	for (long int j = 0; j < in->N; j++) {
		P = 1.0e-6*(in->y[0][j]);
		F = (in->y[1][j])/par.sound;
		fprintf(fh_pipe, "%.12e\t%.12e\t%.12e\n", 1.0e-3*j*in->dx, P, F);
	}
	fclose(fh_pipe);
}

void dump_node(gnode_ptr in, char *line) {
	fh_node = fopen(line, "a");
	P = 1.0e-6*(in->p[0]);
	fprintf(fh_node, "%.12e\t%.12e\n", par.curr*par.dt, P);
	fclose(fh_node);
}