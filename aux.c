#include "network.h"

void err_msg(char* msg) {
  printf("%s\n", msg);
  exit(1);
}

void debug_msg(char* msg) {
  if (DEBUGMODE) printf("%s", msg);
  
}

int noise_status(noise_p *noise) {
  char msg[80];
  debug_msg("\nNoise parameters:\n\n");
  sprintf(msg, "Tau = %.8e\nL_c = %.8e\nAmp = %.8e\n", noise->tc, noise->lc, noise->A);
  debug_msg(msg);
  return 0;
}
