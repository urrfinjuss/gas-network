#include "network.h"

void init_links(network *net) {
   FILE *fh = fopen(net->matname,"r");
   size_t *nbytes;
   char line[256];
   int *link_matrix, rnum, m = 0;
   link_matrix = malloc(sizeof(int)*pow(net->nnodes,2));
   while(!feof(fh)){
	rnum = fscanf(fh,"%s", line);
	while(1){
          sscanf(line,"%d",&link_matrix[m]);  //reading a single digit
	  m++;
          /*Store or use float_num*/                                                              
          if(strchr(line,' ')) strcpy(line , strchr(line,' ')+1 );  
            else break;  //exiting once the digits is empty
        }
   }
   printf("\nLink Matrix %d elements\n", m-1);
   int nlinks2 = 0;
   for (int j = 0; j < (net->nnodes)*(net->nnodes); j++) nlinks2 += link_matrix[j];
   if (nlinks2 != net->nlinks) err_msg("Number of links in the matrix contradicts input file.");
   allocate_memory(net, link_matrix);
   /*for (int j = 0; j < net->nnodes; j++) {
	for (int k = 0; k < net->nnodes; k++) printf("%d\t", link_matrix[j*(net->nnodes)+k]);
	printf("\n");
   }*/
   printf("\nTotal links in network:  %d\n", net->nlinks);
   err_msg("Complete");	
}

void allocate_memory(network *net, int *lm) {
  net->knot = malloc(sizeof(node_ptr)*(net->nnodes));
  net->link = malloc(sizeof(pipe_ptr)*(net->nlinks));
  node_ptr p_k;
  int l;
  for (int i = 0; i < net->nnodes; i++) {
    net->knot[i] = malloc(sizeof(node));
    p_k = net->knot[i];
    p_k->adj_n = 0;
    for (int j = 0; j < net->nnodes; j++) {
	(p_k->adj_n) += lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j];
    }
    (net->knot[i])->adj_k = malloc(sizeof(node_ptr)*(p_k->adj_n));
    //printf("Node %d has %d adjacent nodes\n", i, p_k->adj_n);
  }
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    printf("Node %d neigbors:\t", i);
    l = 0;
    for( int j = 0; j < net->nnodes; j++) {
	if ((lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j]) == 1) {
	  (net->knot[i])->adj_k[l] = net->knot[j];
	  printf("%p\t", net->knot[j]);
	  l++;
	}
    }
    printf("\n");	
  }
  printf("\n");

  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    printf("Node %d my neighbors:\t", i);
    for (int k = 0; k < (p_k->adj_n); k++){
      printf("%p\t", p_k->adj_k[k]);
    }
    printf("\n");
  }




















  /*int my_neighbors;
  for (int i = 0; i < net->nlinks; i++) {
    p_k = net->knot[i];
    my_neighbors = p_k->adj_n;
    (net->knot[i])->adj_p = malloc(sizeof(pipe_ptr)*my_neighbors);
    printf("Generated %d pipes connected to node %d\n", my_neighbors, i);

  }*/
}











