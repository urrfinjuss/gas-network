#include "network.h"

void err_msg(char* msg) {
  printf("%s\n", msg);
  exit(1);
}

void debug_msg(char* msg) {
  if (DEBUGMODE) printf("%s", msg);
  
}

void copy_network(network_ptr ns, network_ptr nd) {
  gpipe_ptr tmp1, tmp2;
  if (ns->nlinks != nd->nlinks) err_msg("Function copy_network got incompatible networks.\n");
  nd->curr_T = ns->curr_T;
  for (int j = 0; j < ns->nlinks; j++) {
	tmp1 = ns->link[j];
	tmp2 = nd->link[j];
	for (int k = 0; k < tmp1->N; k++){
          tmp2->p[k] = tmp1->p[k];
	  tmp2->f[k] = tmp1->f[k];
	}
  }
}

void view_network(network_ptr ns){
  printf("Network name %s:\n", ns->fname);
  printf("Nodes in network %d:\t", ns->nnodes);
  for (int j = 0; j < ns->nnodes; j++) printf("%p\t", ns->knot[j]);
  printf("\n");
  printf("Links in network %d:\t", ns->nlinks);
  for (int j = 0; j < ns->nlinks; j++) printf("%p\t", ns->link[j]);
  printf("\n");
  printf("Compressors in network %d\n", ns->ncomps);
  for (int j = 0; j < ns->nlinks; j++) if (ns->link[j]->c_id != NULL ) printf("%p\t", ns->link[j]->c_id);
  printf("\n");
  printf("Number of points per km %d\n", ns->npcent);
  
}

#if 0
struct network {
  char fname[80], nname[80], dname[80], current_dir[1024];
  char incname[256];
  int nnodes;
  int nlinks;
  int ncomps;
  int npcent;
  int mglf;
  int nskip;
  double c, diss, tmax;
  noise_p *noise;
  node_ptr *knot;
  gpipe_ptr *link;
  compressor_ptr cssr;
} *network_ptr;
#endif


