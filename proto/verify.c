#include "header.h"

static char text[80];

void verify_input_conf() {
	dmesg("verify_input_conf:\n", 0);
	printf("\t[Configuration]\n");
	printf("\tname=\t%s\n", par.name);
	printf("\tskip=\t%d\n", par.skip);
	printf("\ttmax=\t%.8e\n", par.tmax/3600.);
	printf("\tsounds=\t%.8e\n", par.sound);
	printf("\tDWdiss=\t%.8e\n", par.friction);
	printf("\tinterp=\t%s\n", par.intmethod);
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
void verify_consistency() {
	dmesg("verify_consistency:\n", 0);
	for (long int j = 0; j < net.nodes; j++) {
		sprintf(text, "Node %ld:\n", j); dmesg(text, 0);
		sprintf(text, "type = %u\n", net.node[j].type); dmesg(text, 0);
		sprintf(text, "ncomp = %u\n", net.node[j].ncomp); dmesg(text, 0);
		for (int jk = 0; jk < net.node[j].ncomp; jk++) {
		  sprintf(text, "comp = %p\n", &net.node[j].comp[jk]); dmesg(text, 0);
		}
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
	dmesg("verify_consistency:\tPassed\n", 0);
}
