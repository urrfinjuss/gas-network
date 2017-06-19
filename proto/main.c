#include "header.h"

network	net;
params	par;
static char text[80];

int main(int argc, char** argv) {
	printf("Gas Network Toolbox\n");
	if (argc != 2) {
		printf("Usage:\n%s conf.cfg\n", argv[0]);
		exit(0);
	} else call_init_network(argv[1]);

	/*
	load_initial_data();
	verify_consistency();
	network_snapshot();

	call_init_temporal();
	call_init_nodes();
	call_init_comps();
	*/

	return 0;
}
