#include "network.h"

static network syst;

int init_network(char *filename, network *net) {
  FILE* fh = fopen(filename,"r");
  char line[80], param[80], value[80];
  noise_p *nse = net->noise;
  nse = malloc(sizeof(noise_p));

  if (fh == NULL) err_msg("Cannot open configuration file.");
  else while (fgets(line, 80, fh) != NULL) {
    sscanf(line, "%s\t%s", param, value);
    if (strcmp(param,"#name=") == 0) sprintf(net->fname,"%s", value);
    if (strcmp(param,"#nmatr=") == 0) sprintf(net->matname,"%s", value);
    if (strcmp(param,"#nodes=") == 0) net->nnodes = atoi(value);
    if (strcmp(param,"#links=") == 0) net->nlinks = atoi(value);
    if (strcmp(param,"#corr_time=") == 0) nse->tc = atof(value);
    if (strcmp(param,"#corr_dist=") == 0) nse->lc = atof(value);
    if (strcmp(param,"#amplitude=") == 0) nse->A = atof(value);
  }
  fclose(fh);
  noise_status(nse);
  return 1;
}



