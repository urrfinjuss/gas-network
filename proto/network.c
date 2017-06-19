#include "header.h"

static char text[80];



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

