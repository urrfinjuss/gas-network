#include "header.h"

network	net;
params	par;
static char text[80];

int main(int argc, char** argv) {
	printf("Gas Network Toolbox v2.0.0: ");
	if (!OS) printf("Linux detected.\n");
	else printf("Windows detected.\n");
	if (argc != 2) {
		printf("Usage:\n%s conf.cfg\n", argv[0]);
		exit(0);
	} 
	call_init_network(argv[1]);

	evolve_network();

	free_forward_interior();
	return 0;
}
