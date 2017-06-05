#include "header.h"

static char text[80];


int read_input_file(char *fname) {
	FILE *fh = fopen(fname, "r");
	if (fh) {
		char line[80], str[80], value[80];
		while (fgets(line, 80, fh) != NULL) {
			sscanf(line, "%s\t%s", str, value);
			if (strcmp(str, "name=") == 0)   sprintf(par.name, "%s", value);
			if (strcmp(str, "skip=") == 0)   par.skip = atoi(value);
			if (strcmp(str, "tmax=") == 0)   par.tmax = atof(value);
			if (strcmp(str, "sounds=") == 0) par.sound = atof(value);
			if (strcmp(str, "DWdiss=") == 0) par.friction = atof(value);
		}
		fclose(fh);
		return 0;
	} else {
		dmesg("Cannot open configuration file.\n" ,1);
		return 1;
	}
}


int read_network_list(char *fname) {
	FILE *fh = fopen(fname, "r");
	if (fh) {
	  char line[256], str[256];
		if ( fgets(line, 256, fh) != NULL) {
			sscanf(line, "%u\t%u\t%u", &net.nodes, &net.pipes, &net.compressors);
			sprintf(str, "%sNetwork:\nNodes\t%u\nPipes\t%u\n", line, net.nodes, net.pipes);
			dmesg(str, 0);
			net.node = malloc(net.nodes*sizeof(gnode));
			net.pipe = malloc(net.pipes*sizeof(gpipe));
		}
		for (int j = 0; j < net.nodes; j++) {
		  if (fgets(line, 256, fh) != NULL) dmesg(line, 0);
			net.node[j].left  = NULL;
			net.node[j].right = NULL;
			net.node[j].type = 0;
		}
		for (int j = 0; j < net.pipes; j++) {
			unsigned int L_Id, R_Id;
			unsigned int p_Id;
			FTYPE len, wid;
			if (fgets(line, 256, fh) != NULL) {
				sscanf(line, "%u\t%u\t%u\t%lf\t%lf", &p_Id, &L_Id, &R_Id, &wid, &len);
				dmesg(line,0);
			}
			net.pipe[j].left  = &net.node[L_Id];
			net.pipe[j].right = &net.node[R_Id];
			net.pipe[j].len   = 1000.*len;
			net.pipe[j].wid   = wid;
		}
		return 0;
	} else return 1;
}

int read_pipe(char *fname, gpipe_ptr pipe) {
	FILE *fh = fopen(fname, "r");
	if (fh) {
		char line[256];
		FTYPE x;
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		 int n = 0;
		dmesg(line, 0);
		while ( fgets(line, 256, fh) != NULL ) n++;
		//sprintf(text, "Found %ld points\n", n-1);
		//dmesg(text, 0);
		fclose(fh);
		pipe->N = n;
		pipe->dx = pipe->len/(n-1);
		pipe->y[0] = malloc(n*sizeof(FTYPE));
		pipe->y[1] = malloc(n*sizeof(FTYPE));
		fh = fopen(fname, "r");
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		fgets(line, 256, fh);
		long int j = 0;
		while ( fgets(line, 256, fh) != NULL) {
			sscanf(line, "%lf\t%lf\t%lf", &x, &pipe->y[0][j], &pipe->y[1][j]);
			pipe->y[0][j] = 1.0e+6*pipe->y[0][j];
			pipe->y[1][j] = par.sound*pipe->y[1][j];
			j++;
			//printf("%.8e\t%.8e\n", pipe->y[0][j], pipe->y[1][j]);
		}
		sprintf(text, "Discretization step %f\n", pipe->dx);
		dmesg(text, 0);
		return 0;
	} else return 1;
}

int install_compressors(char *fname) {
	FILE *fh = fopen(fname, "r");
	if (fh) {
		char line[256];
		unsigned int m,n,p;
		int counter = 0;
		while (counter < net.nodes+net.pipes+1) {
			fgets(line, 256, fh);
			counter++;
		}
		net.comp = malloc(net.compressors*sizeof(gcomp_ptr));
		for (long int j = 0; j < net.compressors; j++) {
		  if (fgets(line, 256, fh) != NULL) dmesg(line, 0);
			sscanf(line, "%u\t%u\t%u", &m, &n, &p);
			net.comp[m].loc  = &net.node[n];
			net.comp[m].dest = &net.pipe[p];
		}
		int nc;
		for (long int k = 0; k < net.nodes; k++) {
		  nc = 0;
			for (long int j = 0; j < net.compressors; j++) {
				if (net.comp[j].loc == &net.node[k]) nc++;
			}
			net.node[k].ncomp = nc;
			if (nc != 0) net.node[k].comp = malloc(nc*sizeof(gcomp_ptr));
		}
		printf("Building compressors into the network\n");
		for (long int k = 0; k < net.nodes; k++) {
			long int n = 0;
		  for (long int j = 0; j < net.compressors; j++) {
				if (net.comp[j].loc == &net.node[k]) {
					sprintf(line, "Compressor installed: node %p to %p\n", &net.node[k], net.comp[j].dest); 
					dmesg(line, 0);
					net.node[k].comp[n] = &net.comp[j];
					n++;
				}
			}
		}
		fclose(fh);
		return 0;
	} else return 1;
}

int load_initial_data() {
	char	fname[80];
	for (long int j = 0; j < net.pipes; j++) {
		sprintf(fname, "./ic/data_%03ld.txt", j);
		FILE *fh = fopen(fname, "r");
		if (fh) {
			sprintf(text, "Found data for pipe %ld:", j);
			dmesg(text, 0);
			fclose(fh);
		} else dmesg("Missing data for one of the pipes\n", 1);
		if (!read_pipe(fname, &net.pipe[j])) dmesg("Loaded\n", 0);
		
	}
	return 0;
}
