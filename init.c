#include "network.h"

int init_network(char *filename, network *net) {
  FILE* fh = fopen(filename,"r");
  char line[80], param[80], value[80], msg[256];
  noise_p *nse = net->noise;
  nse = malloc(sizeof(noise_p));
  strcpy(net->nname, filename);

  if (fh == NULL) err_msg("Cannot open configuration file.");
  else while (fgets(line, 80, fh) != NULL) {
    sscanf(line, "%s\t%s", param, value);
    if (strcmp(param,"#name=") == 0) sprintf(net->fname,"%s", value);
    if (strcmp(param,"#nmatr=") == 0) sprintf(net->matname,"%s", value);
    if (strcmp(param,"#nprof=") == 0) sprintf(net->dname,"%s", value);
    if (strcmp(param,"#nodes=") == 0) net->nnodes = atoi(value);
    if (strcmp(param,"#npts=") == 0) net->npcent = atoi(value);
    if (strcmp(param,"#links=") == 0) net->nlinks = atoi(value);
    if (strcmp(param,"#corr_time=") == 0) nse->tc = atof(value);
    if (strcmp(param,"#corr_dist=") == 0) nse->lc = atof(value);
    if (strcmp(param,"#amplitude=") == 0) nse->A = atof(value);
    if (strcmp(param,"#sound_speed=") == 0) net->c = atof(value);
    if (strcmp(param,"#dissip_coef=") == 0) net->diss = atof(value);
    if (strcmp(param,"#simul_time=") == 0) net->tmax = 3600.*atof(value);

  }
  fclose(fh);
  noise_status(nse);
  sprintf(msg, "%s net->fname\n%s net->matname\n%s net->dname\n%d net->nnodes\n%d net->npcent\n%d net->nlinks\n", 
		net->fname, net->matname, net->dname, net->nnodes, net->npcent, net->nlinks);  
  debug_msg(msg);

  return 1;
}


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
   sprintf(msg, "Link Matrix %d elements\n", m-1);
   debug_msg(msg);
   int nlinks2 = 0;
   for (int j = 0; j < (net->nnodes)*(net->nnodes); j++) nlinks2 += link_matrix[j];
   if (nlinks2 != net->nlinks) err_msg("Number of links in the matrix contradicts input file.");
   allocate_memory(net, link_matrix);
   init_arrays(net);
   sprintf(msg, "\nTotal links in network:  %d\n", net->nlinks);
   debug_msg(msg);
   //err_msg("Complete");	
}

void allocate_memory(network *net, int *lm) {
  char msg[256];
  node_ptr p_k, o_k;
  gpipe_ptr l_p;
  int l, m, n, rnum;
  net->knot = malloc(sizeof(node_ptr)*(net->nnodes));
  net->link = malloc(sizeof(gpipe_ptr)*(net->nlinks));
  memset(net->link, 0, (net->nlinks)*sizeof(gpipe_ptr));
  memset(net->knot, 0, (net->nnodes)*sizeof(node_ptr));
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
    //p_k->nr = 0;
    //p_k->nl = 0;
    for (int j = 0; j < net->nnodes; j++) {
	(p_k->adj_n) += lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j];
	//(p_k->nl) += lm[j*(net->nnodes)+i];
        //(p_k->nr) += lm[i*(net->nnodes)+j];
    }
    (net->knot[i])->adj_k = malloc(sizeof(node_ptr)*(p_k->adj_n));
    (net->knot[i])->W  = malloc(sizeof(double)*(p_k->adj_n));
    (net->knot[i])->Wb = malloc(sizeof(double)*(p_k->adj_n));
  }
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    l = 0;
    for( int j = 0; j < net->nnodes; j++) {
	o_k = net->knot[j];
	if ((lm[j*(net->nnodes)+i]+lm[i*(net->nnodes)+j]) == 1) {
	  p_k->adj_k[l] = o_k;
	  //(net->knot[i])->adj_k[l] = net->knot[j];
	  l++;
	}
        //p_k->adj_p = malloc(l*sizeof(gpipe_ptr));
        //(net->knot[j])->adj_p = malloc(l*sizeof(gpipe_ptr));
    }
    p_k->adj_p = malloc((p_k->adj_n)*sizeof(gpipe_ptr));
    memset(p_k->adj_p, 0, (p_k->adj_n)*sizeof(gpipe_ptr));
    sprintf(msg, "Knot %d connected via %d pipes (adj_n adjacent nodes  = %d) at %p\n", i, l, p_k->adj_n, p_k);
    debug_msg(msg);
    //sprintf(msg, "Knot %d has %d connections from the left and %d from the right\n", i, p_k->nl, p_k->nr);
    //debug_msg(msg);
  }
  int nlinks = 0;
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    for (int j = 0; j < p_k->adj_n; j++) {
      p_k->W[j] = 0.;
      p_k->Wb[j] = 0.;
      //printf("%d of %d\t%p\n", j, p_k->adj_n, net->link[nlinks]);
      if ((p_k->adj_p[j]) == NULL) {
	net->link[nlinks] = malloc(sizeof(gpipe));
	p_k->adj_p[j] = net->link[nlinks];
	//p_k->adj_p[j] = malloc(sizeof(gpipe));
        l = 0;
	while(1) {
          if ((p_k->adj_k[j])->adj_p[l] == NULL) {
	    (p_k->adj_k[j])->adj_p[l] = p_k->adj_p[j];


	    (net->link[nlinks])->left = p_k;
	    (net->link[nlinks])->right = p_k->adj_k[j];
		
	    //((p_k->adj_k[j])->adj_p[l])->left = net->knot[i];
	    //((p_k->adj_k[j])->adj_p[l])->right = p_k->adj_k[j];
	    //printf("%p\n", ((p_k->adj_k[j])->adj_p[l])->right);
            break;
	  }
	  l++;
	}
	nlinks++;
      }	
      sprintf(msg, "Node (%d) at %p connects to node %p via link %p\n", i, p_k, p_k->adj_k[j], p_k->adj_p[j]);
      debug_msg(msg);
    }
    debug_msg("\n");
  }
  sprintf(msg, "Total number of links in the network %d\n", nlinks);
  debug_msg(msg);

  for (int i = 0; i < net->nlinks; i++) {
    l_p = net->link[i];
    sprintf(msg, "Link %2d (%p) left side connected to %p right side to %p\n", i, l_p, l_p->left, l_p->right);
    debug_msg(msg);
  }

  FILE *fh = fopen(net->nname,"r");
  char *dm;
  char msg2[80];

  debug_msg(net->nname);
  debug_msg("\n");
  while(1) {
    if(fgets(msg, 256, fh)==NULL) break;
    if( strcmp(msg, "#lengths\n")==0) {
	dm = fgets(msg, 256, fh);
        sprintf(msg2, "Found pipe lengths\n");
	debug_msg(msg2);
	debug_msg(msg);
	(net->link[0])->L = strtod(msg, &dm);
	sprintf(msg2, "Length of link %2d is %e\n", 0, (net->link[0])->L);
	debug_msg(msg2);
        for (int j = 1; j < net->nlinks; j++) {
	  if (strcmp(dm, "\n") == 0) err_msg("Input file has incorrect number of pipe lengths");
	  (net->link[j])->L = strtod(dm, &dm);
	  sprintf(msg2, "Length of link %2d is %e\n", j, (net->link[j])->L);
	  debug_msg(msg2);
	}
    }
    if( strcmp(msg, "#diameters\n")==0) {
	dm = fgets(msg, 256, fh);
        sprintf(msg2, "Found pipe diameters\n");
	debug_msg(msg2);
	debug_msg(msg);
	(net->link[0])->d = strtod(msg, &dm);
	sprintf(msg2, "Diameter of link %2d is %e\n", 0, (net->link[0])->d);
	debug_msg(msg2);
        for (int j = 1; j < net->nlinks; j++) {
	  if (strcmp(dm, "\n") == 0) err_msg("Input file has incorrect number of pipe diameters");
	  (net->link[j])->d = strtod(dm, &dm);
	  sprintf(msg2, "Diameter of link %2d is %e\n", j, (net->link[j])->d);
	  debug_msg(msg2);
	}
    }
    if( strcmp(msg, "#compressions\n")==0) {
	dm = fgets(msg, 256, fh);
        sprintf(msg2, "Found compressions in nodes\n");
	debug_msg(msg2);
	debug_msg(msg);
	(net->knot[0])->cratio = strtod(msg, &dm);
	sprintf(msg2, "Compression in node %2d is %e\n", 0, (net->knot[0])->cratio);
	debug_msg(msg2);
        for (int j = 1; j < net->nnodes; j++) {
	  if (strcmp(dm, "\n") == 0) err_msg("Input file has incorrect number of compression ratios in nodes");
	  (net->knot[j])->cratio = strtod(dm, &dm);
	  sprintf(msg2, "Compression in node %2d is %e\n", j, (net->knot[j])->cratio);
	  debug_msg(msg2);
	}
    }
  }
  fclose(fh);
}

void init_arrays(network *net) {
  int npts;
  gpipe_ptr p_k;

  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];	
    p_k->N = (p_k->L)*(net->npcent);
    p_k->p = malloc((p_k->N)*sizeof(double));
    p_k->f = malloc((p_k->N)*sizeof(double));
    p_k->q = malloc((p_k->N)*sizeof(double));				
  }

}

void init_data(network *net) {
  char str[256], msg[256];
  FILE *fh;
  for (int j = 0; j < net->nlinks; j++) {
	sprintf(str, "%s_%03d.txt", net->dname, j);
	//debug_msg(str);
	fh = fopen(str, "r");
	if (fh == NULL) {
		sprintf(str, "Missing initial data along pipe %d, %d points expected\n", j, (net->link[j])->N);
		debug_msg(str);
	} else {
		if (load_data(fh, net->link[j])) {
		   sprintf(msg, "Data for pipe %d has wrong number of lines.\n%d data points expected\n", j, (net->link[j])->N);
		   printf("%s", msg);
		   fclose(fh);
		   err_msg("Complete");
 		}
		fclose(fh);
	}
  }
  rescale_data(net);
}


int load_data(FILE *fh, gpipe_ptr lnk) {
  char line[256], val0[64], val1[64], val2[64], val3[64];
  char msg[80];
  char *dm;
  dm = fgets(line, 80, fh);
  dm = fgets(line, 80, fh);
  dm = fgets(line, 80, fh);

  int k = 0;
  while (fgets(line, 80, fh) != NULL) {
	sscanf(line,"%s\t%s\t%s\t%s\n", val0, val1, val2, val3);
	lnk->p[k] = strtod(val1, NULL);
	lnk->f[k] = strtod(val2, NULL);
	lnk->q[k] = strtod(val3, NULL);
	k++;
	if (k == lnk->N) {
	   sprintf(msg, "Read %d lines\n", lnk->N);
	   debug_msg(msg);
	   break;
	}
  }
  if (k != lnk->N) return 1;
  else return 0;
}

void rescale_data(network *net) {
  gpipe_ptr p_k;
  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];
    for(int k = 0; k < p_k->N; k++) {
      p_k->p[k] = 1000.*(p_k->p[k]);
      p_k->f[k] = (p_k->f[k])*(net->c);
      p_k->q[k] = p_k->q[k];
    }
  }
}

void save_data(FILE* fh, gpipe_ptr lnk) {
  if (fh == NULL) err_msg("Cannot write to file.");
}










