#include "header.h"

static char text[80];


int call_init_network(char *filename) {
	if (!read_input_file(filename))     dmesg("read_input_file\tPassed\n", 0);
	if (!read_network_list(par.name))   dmesg("read_network_list:\tPassed\n", 0);
	if (!build_network())               dmesg("build_network:\t\tPassed\n", 0);
	if (!install_compressors(par.name)) dmesg("install_compressors:\t\tPassed\n", 0);
	return 0;
}

int verify_consistency() {
	dmesg("Verifying Nodes:\n", 0);
	for (long int j = 0; j < net.nodes; j++) {
		sprintf(text, "Node %ld:\n", j); dmesg(text, 0);
		sprintf(text, "node_type = %s\n", net.node[j].node_type); dmesg(text, 0);
		sprintf(text, "type = %u\n", net.node[j].type); dmesg(text, 0);
		sprintf(text, "ncomp = %u\n", net.node[j].ncomp); dmesg(text, 0);
		sprintf(text, "comp = %p\n", net.node[j].comp); dmesg(text, 0);
		if (net.node[j].comp != NULL) {
			for (long int k = 0; k < net.node[j].ncomp; k++) {
				//verify_compressor();
				sprintf(text, "compressor (%ld of %u) at node %ld\n", k+1, net.node[j].ncomp, j); 
				dmesg(text, 0);
			}
		}
		sprintf(text, "nleft = %u\n", net.node[j].nleft); dmesg(text, 0);
		if (net.node[j].left != NULL) {
			for (long int k = 0; k < net.node[j].nleft; k++) {
				sprintf(text, "incoming (%ld at %p) at node %ld\n", k, net.node[j].left[k], j); 
				dmesg(text, 0);
			}
		}
		sprintf(text, "nright = %u\n", net.node[j].nright); dmesg(text, 0);
		if (net.node[j].right != NULL) {
			for (long int k = 0; k < net.node[j].nright; k++) {
				sprintf(text, "outgoing (%ld at %p) at node %ld\n", k, net.node[j].right[k], j); 
				dmesg(text, 0);
			}
		}
		dmesg("\n", 0);
	}
	dmesg("Nodal Information Complete\n", 0);
	return 0;
}

int build_network() {
	unsigned int left, right;
	printf("Scanning for connections\n");
	for (long int k = 0; k < net.nodes; k++) {
		left = 0;
		right = 0;
		for (long int j = 0; j < net.pipes; j++) {
			if ( net.pipe[j].left  == &net.node[k]) left++;
			if ( net.pipe[j].right == &net.node[k]) right++;
		}
		if (left != 0) {
			net.node[k].right = malloc(left*sizeof(gpipe_ptr));
			net.node[k].nright = left;
		} else net.node[k].nright = 0;
		if (right != 0) {
			net.node[k].left = malloc(right*sizeof(gpipe_ptr));
			net.node[k].nleft = right;
		} else net.node[k].nleft = 0;
	}
	printf("Building nodal connections.\n");
	for (long int k = 0; k < net.nodes; k++) {
		long int n = 0, m = 0;
		for (long int j = 0; j < net.pipes; j++) {
			if (net.pipe[j].left == &net.node[k]) {
				net.node[k].right[n] = &net.pipe[j];
				n++;
			}
			if (net.pipe[j].right == &net.node[k]) {
				net.node[k].left[m] = &net.pipe[j];
				m++;
			}
		}
	}
	if (DEBUG_MODE) {
		printf("Network nodal list:\n");
		for (long int k = 0; k < net.nodes; k++) {
			printf("Node %ld at %p:\n", k, &net.node[k]);
			printf("Inbound connections:\t%u to", net.node[k].nleft);
			for (long int j = 0; j < net.node[k].nleft; j++) {
				printf(" %p", net.node[k].left[j]);
			}
			printf("\n");
			printf("Outbound connections:\t%u from", net.node[k].nright);
			for (long int j = 0; j < net.node[k].nright; j++) {
				printf(" %p", net.node[k].right[j]);
			}
			printf("\n");
		}
	}
	return 0;
}

