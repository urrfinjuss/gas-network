#include "header.h"

static char text[80];
static char line[256];
static char str[256];

void read_network_list(char *fname) {
	dmesg("read_network_list:\n", 0);
	FILE *fh = fopen(fname, "r");
	if (!fh) dmesg("Network List File missing\n", 1);
	if (fgets(line, 256, fh) != NULL) {
		dmesg(line, 0);
		sscanf(line, "%u\t%u\t%u", &net.nodes, &net.pipes, &net.comps);
		net.node = malloc(net.nodes*sizeof(gnode));
		net.pipe = malloc(net.pipes*sizeof(gpipe));
		net.comp = malloc(net.comps*sizeof(gcomp));
	}
	fclose(fh);
	read_nodes_list(fname);
	read_pipes_list(fname);
	read_comps_list(fname);
	update_nodes_comps();
	update_nodes_pipes();
	dmesg("read_network_list:\tPassed\n", 0);
}

void read_nodes_list(char *fname) {
	dmesg("\tread_nodes_list:\n", 0);
	FILE *fh = fopen_set(fname, "r", 1);
	unsigned int num;
	FTYPE v1, v2;
	for (int j = 0; j < net.nodes; j++) {
		if (fgets(line, 256, fh) != NULL) {
			dmesg(line, 0);
			sscanf(line, "%u\t%lf\t%lf\t%u", &num, &v1, &v2, &net.node[j].type);
		}
		net.node[j].left  = NULL;
		net.node[j].right = NULL;
		net.node[j].ncomp = 0;
		net.node[j].nleft = 0;
		net.node[j].nright = 0;
	}
	fclose(fh);
	dmesg("\tread_nodes_list:\tPassed\n", 0);
}

void read_pipes_list(char *fname) {
	dmesg("\tread_pipes_list:\n", 0);
	FILE *fh = fopen_set(fname, "r", net.nodes+1);
	unsigned int L_Id, R_Id, p_Id;
	FTYPE length, width;
	for (int j = 0; j < net.pipes; j++) {
		if (fgets(line, 256, fh) != NULL) {
			dmesg(line, 0);
			sscanf(line, "%u\t%u\t%u\t%lf\t%lf", &p_Id, &L_Id, &R_Id, &width, &length);
		}
		net.pipe[j].left  = &net.node[L_Id];
		net.pipe[j].right = &net.node[R_Id];
		net.pipe[j].len   = 1.0e+3*length;
		net.pipe[j].wid   = width;
		net.node[L_Id].nright++;
		net.node[R_Id].nleft++;
	}
	fclose(fh);
	dmesg("\tread_pipes_list:\tPassed\n", 0);
}

void read_comps_list(char *fname) {
	dmesg("\tread_comps_list:\n", 0);
	FILE *fh = fopen_set(fname, "r", net.nodes+net.pipes+1);
	unsigned int m, n, p;
	for (int j = 0; j < net.comps; j++) {
		if (fgets(line, 256, fh) != NULL) {
			dmesg(line, 0);
			sscanf(line, "%u\t%u\t%u", &m, &n, &p);
		}
		net.comp[j].loc  = &net.node[n];
		net.comp[j].dest = &net.pipe[p];
		net.node[n].ncomp++;
	}
	fclose(fh);
	dmesg("\tread_comps_list:\tPassed\n", 0);
}

void update_nodes_comps() {
	dmesg("\tupdate_nodes_comps:\n", 0);
	for (int j = 0; j < net.nodes; j++) {
		if (net.node[j].ncomp) {
			net.node[j].comp = malloc(net.node[j].ncomp*sizeof(gcomp_ptr));
			for (int k = 0; k < net.node[j].ncomp; k++) net.node[j].comp[k] = NULL;
		}
	}
	for (int j = 0; j < net.comps; j++) {
		int k = 0;
		while (k < net.comp[j].loc->ncomp) {
			if (net.comp[j].loc->comp[k] == NULL) {
				printf("Compressor %d placed to slot %d on node %p\n", j, k, net.comp[j].loc);
				net.comp[j].loc->comp[k] = &net.comp[j];
				break;
			}
			k++;
		}
	}
	dmesg("\tupdate_nodes_comps:\tPassed\n", 0);
}

void update_nodes_pipes() {
	dmesg("\tupdate_nodes_pipes:\n", 0);
	for (int j = 0; j < net.nodes; j++) {
		if (net.node[j].nleft) {
			net.node[j].left = malloc(net.node[j].nleft*sizeof(gpipe_ptr));
			for (int k = 0; k < net.node[j].nleft; k++) net.node[j].left[k] = NULL;
		}
		if (net.node[j].nright) {
			net.node[j].right = malloc(net.node[j].nright*sizeof(gpipe_ptr));
			for (int k = 0; k < net.node[j].nright; k++) net.node[j].right[k] = NULL;
		}
	}
	for (int j = 0; j < net.pipes; j++) {
		int l = 0, r = 0;
		while (l < net.pipe[j].right->nleft) {
			if (net.pipe[j].right->left[l] == NULL) {
				printf("Pipe %d incoming to slot %d on node %p\n", j, l, net.pipe[j].right);
				net.pipe[j].right->left[l] = &net.pipe[j];
				break;
			}
			l++;
		}
		while (r < net.pipe[j].left->nright) {
			if (net.pipe[j].left->right[r] == NULL) {
				printf("Pipe %d outgoing from slot %d on node %p\n", j, r, net.pipe[j].left);
				net.pipe[j].left->right[r] = &net.pipe[j];
				break;
			}
			r++;
		}
	}
	dmesg("\tupdate_nodes_pipes:\tPassed\n", 0);
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
			printf("%s", line);
		}
		net.comp = malloc(net.comps*sizeof(gcomp_ptr));
		for (long int j = 0; j < net.comps; j++) {
		  if (fgets(line, 256, fh) != NULL) dmesg(line, 0);
			sscanf(line, "%u\t%u\t%u", &m, &n, &p);
			net.comp[m].loc  = &net.node[n];
			net.comp[m].dest = &net.pipe[p];
			printf("Compressor %d (%p) located at %p\n", m, &net.comp[m], &net.node[n]);
		}
		for (long int j = 0; j < net.comps; j++) {
			printf("Cmp %ld (%p): Located at %p\n", j, &net.comp[j], net.comp[j].loc);
		}
		int nc;
		for (long int k = 0; k < net.nodes; k++) {
		  nc = 0;
			for (long int j = 0; j < net.comps; j++) {
				if (net.comp[j].loc == &net.node[k]) nc++;
			}
			net.node[k].ncomp = nc;
			printf("Node %ld has %d compressors\n", k, nc);
			if (nc != 0) {
				net.node[k].comp = malloc(nc*sizeof(gcomp_ptr));
				//for (int l = 0; l < nc; l++) (net.node[k].comp[l]) = NULL;
			} else (net.node[k].comp) = NULL;
		}
		for (long int j = 0; j < net.comps; j++) {
			printf("Cmp %ld (%p): Located at %p\n", j, &net.comp[j], net.comp[j].loc);
		}
		printf("Building compressors into the network\n");
		long int nn;
		for (long int k = 0; k < net.nodes; k++) {
			nn = 0;
		  for (long int j = 0; j < net.comps; j++) {
				printf("node %p    comp.loc %p\n", &net.node[k], net.comp[j].loc);
				if (net.comp[j].loc == &net.node[k]) {
					sprintf(line, "Compressor %p installed at node %p to %p\n", &net.comp[j], &net.node[k], net.comp[j].dest); 
					dmesg(line, 0);
					net.node[k].comp[nn] = &net.comp[j];
					nn++;
				}
			}
			if (nn != net.node[k].ncomp) dmesg("Error building the network\n",1);
			else { 
				sprintf(text, "Node %ld at %p passed successfully\n", k, &net.node[k]);
				dmesg(text, 0);
			}
		}
		fclose(fh);
		return 0;
	} else return 1;
}


FILE* fopen_set(char* fname, char* flag, const unsigned skip_lines) {
  char line[80];
  FILE *fh = fopen(fname, flag);
	if (fh) {
		for (int j = 0; j < skip_lines; j++) fgets(line, 80, fh);
	} else dmesg("Cannot open file\n" ,1);
  return fh;
}

long int count_data(char *fname) {
	char line[80];
	FILE *fh = fopen_set(fname, "r", 3);
	long int nbase = -1;
	while(fgets(line, 80, fh)) nbase++;
	fclose(fh);
	return nbase;
}

