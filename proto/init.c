#include "header.h"

static char text[80];

int call_init_network(char *filename) {
	read_input_file(filename);
	if (DEBUG_MODE) verify_input_conf();

	read_network_list(par.name);
	if (DEBUG_MODE) verify_network_conf();

  exit(0);
	if (!build_network())               dmesg("build_network:\t\tPassed\n", 0);
	if (!install_compressors(par.name)) dmesg("install_compressors:\t\tPassed\n", 0);
	return 0;
}

void read_input_file(char *fname) {
	dmesg("read_input_file:\n", 0);
	FILE *fh = fopen(fname, "r");
	if (fh) {
		char line[80], str[80], value[80];
		while (fgets(line, 80, fh) != NULL) {
			sscanf(line, "%s\t%s", str, value);
			if (strcmp(str, "name=") == 0)   sprintf(par.name, "%s", value);
			if (strcmp(str, "skip=") == 0)   par.skip = atoi(value);
			if (strcmp(str, "tmax=") == 0)   par.tmax = 3600.*atof(value);
			if (strcmp(str, "sounds=") == 0) par.sound = atof(value);
			if (strcmp(str, "DWdiss=") == 0) par.friction = atof(value);
			if (strcmp(str, "interp=") == 0) sprintf(par.intmethod, "%s", value);
		}
		fclose(fh);
	} else dmesg("Cannot open configuration file.\n", 1);
	dmesg("read_input_file:\tPassed\n", 0);
}
