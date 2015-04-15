#include "network.h"

int err_msg(char* msg) {
  printf("%s\n", msg);
  exit(1);
}

int noise_status(noise_p *noise) {
  printf("\nNoise parameters:\n\n");
  printf("Tau = %.8e\nL_c = %.8e\nAmp = %.8e\n", noise->tc, noise->lc, noise->A);

}
