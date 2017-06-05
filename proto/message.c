#include "header.h"

void dmesg(char *line, unsigned int flag) {
	if (DEBUG_MODE) {
		printf("%s", line);
  	if (flag) exit(0);
	}
}
