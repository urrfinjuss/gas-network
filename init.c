#include "network.h"

void init_links(network *net) {
   FILE *fh = fopen(net->matname,"r");
   size_t *nbytes;
   char line[256], msg[80];
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
   sprintf(msg, "\nLink Matrix %d elements\n", m-1);
   debug_msg(msg);
   int nlinks2 = 0;
   for (int j = 0; j < (net->nnodes)*(net->nnodes); j++) nlinks2 += link_matrix[j];
   if (nlinks2 != net->nlinks) err_msg("Number of links in the matrix contradicts input file.");
   allocate_memory(net, link_matrix);
   sprintf(msg, "\nTotal links in network:  %d\n", net->nlinks);
   debug_msg(msg);
   err_msg("Complete");	
}

void allocate_memory(network *net, int *lm) {
  char msg[80];
  node_ptr p_k, o_k;
  int l, m, n;
  net->knot = malloc(sizeof(node_ptr)*(net->nnodes));
  net->link = malloc(sizeof(pipe_ptr)*(net->nlinks));
  debug_msg("Pointers to pipes:\t");
  for (int i = 0; i < net->nlinks; i++) {
    sprintf(msg, "%p\t", net->link[i]);
    debug_msg(msg); 
  }
  debug_msg("\n"); 
  for (int i = 0; i < net->nnodes; i++) {
    net->knot[i] = malloc(sizeof(node));
    p_k = net->knot[i];
    p_k->adj_n = 0;
    for (int j = 0; j < net->nnodes; j++) {
	(p_k->adj_n) += lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j];
    }
    (net->knot[i])->adj_k = malloc(sizeof(node_ptr)*(p_k->adj_n));
  }
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    //sprintf(msg, "Node %d neigbors:\t", i);
    //debug_msg(msg);
    l = 0;
    for( int j = 0; j < net->nnodes; j++) {
	o_k = net->knot[j];
	if ((lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j]) == 1) {
	  (net->knot[i])->adj_k[l] = net->knot[j];
	  //sprintf(msg, "%p\t", net->knot[j]);
   	  //debug_msg(msg);
	  l++;
	}
        p_k->adj_p = malloc(l*sizeof(pipe_ptr));
        //(net->knot[j])->adj_p = malloc(l*sizeof(pipe_ptr));
    }
    
    sprintf(msg, "Knot %d connected via %d pipes (adj_n = %d)\n", i, l, p_k->adj_n);
    debug_msg(msg);
  }


  /*for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    sprintf(msg, "Node %d my neighbors:\t", i);
    debug_msg(msg);
    for (int k = 0; k < (p_k->adj_n); k++){
      sprintf(msg, "%p\t", p_k->adj_k[k]);
      debug_msg(msg);
    }
    sprintf(msg, "\n");
    debug_msg(msg);
  }*/

  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    for (int j = 0; j < p_k->adj_n; j++) {
      //sprintf(msg, "%p\n", p_k->adj_p[j]);
      //debug_msg(msg);
      if ((p_k->adj_p[j]) == NULL) {
	p_k->adj_p[j] = malloc(sizeof(pipe));
        l = 0;
	while(1) {
	  //sprintf(msg, "\n%p\n", (p_k->adj_k[j])->adj_p[l] );
	  //debug_msg(msg);
          if ((p_k->adj_k[j])->adj_p[l] == NULL) {
	    (p_k->adj_k[j])->adj_p[l] = p_k->adj_p[j];
            break;
	  }
	  l++;
	}
      }	
      sprintf(msg, "Node (%d) at %p connects to node %p via link %p\n", i, p_k, p_k->adj_k[j], p_k->adj_p[j]);
      debug_msg(msg);
    }
    debug_msg("\n\n");
  }


















  /*int my_neighbors;
  for (int i = 0; i < net->nlinks; i++) {
    p_k = net->knot[i];
    my_neighbors = p_k->adj_n;
    (net->knot[i])->adj_p = malloc(sizeof(pipe_ptr)*my_neighbors);
    printf("Generated %d pipes connected to node %d\n", my_neighbors, i);

  }*/
}











