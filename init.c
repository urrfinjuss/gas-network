#include "network.h"
static double ovsq2, sq2, cs;
//static compressor_ptr cssr;

int init_network(char *filename, network *net) {
  FILE* fh = fopen(filename,"r");
  char line[80], param[80], value[80], msg[256];
  noise_p *nse = net->noise;
  net->nse = malloc(sizeof(noise_p));
  strcpy(net->nname, filename);
  ovsq2 = 1./sqrt(2.);	
  sq2 = sqrt(2.);

  if (fh == NULL) err_msg("Cannot open configuration file.");
  else while (fgets(line, 80, fh) != NULL) {
    sscanf(line, "%s\t%s", param, value);
    if (strcmp(param,"#name=") == 0) sprintf(net->fname,"%s", value);
    if (strcmp(param,"#nincd=") == 0) sprintf(net->incname,"%s", value);
    if (strcmp(param,"#nprof=") == 0) sprintf(net->dname,"%s", value);
    if (strcmp(param,"#nodes=") == 0) net->nnodes = atoi(value);
    if (strcmp(param,"#ncmps=") == 0) net->ncomps = atoi(value);
    if (strcmp(param,"#drawpng=") == 0) net->mglf = atoi(value);
    if (strcmp(param,"#skipnum=") == 0) net->nskip = atoi(value);
    if (strcmp(param,"#npts=") == 0) net->npcent = atoi(value);
    if (strcmp(param,"#links=") == 0) net->nlinks = atoi(value);
    if (strcmp(param,"#corr_time=") == 0) (net->nse)->tau = atof(value);
    if (strcmp(param,"#corr_dist=") == 0) (net->nse)->lc = atof(value);
    if (strcmp(param,"#amplitude=") == 0) (net->nse)->A = atof(value);
    if (strcmp(param,"#sound_speed=") == 0) net->c = atof(value);
    if (strcmp(param,"#dissip_coef=") == 0) net->diss = atof(value);	
    //if (strcmp(param,"#distr_gamma=") == 0) net->gamma = atof(value);
    if (strcmp(param,"#fix_pressure=") == 0) net->LFLAG = atoi(value);
    if (strcmp(param,"#simul_time=") == 0) net->tmax = 3600.*atof(value);

  }
  fclose(fh);
  cs = net->c;
  noise_status(net->nse);
  sprintf(msg, "%s net->fname\n%s net->incname\n%s net->dname\n%d net->nnodes\n%d net->npcent\n%d net->nlinks\n", 
		net->fname, net->incname, net->dname, net->nnodes, net->npcent, net->nlinks);  
  debug_msg(msg);

  return 1;
}


void init_links(network *net) {
   FILE *fh = fopen(net->incname,"r");
   size_t count = 1000;
   char *line = malloc(1000);
   char msg[80];
   int rnum, m = 0;
   int nlinks2 = 0;
   double *link_matrix = malloc(sizeof(double)*(net->nnodes)*(net->nnodes));


   getline(&line, &count, fh);
   int read = -1, cur = 0, cCount = 0;
   while( sscanf(line+cur, "%lf%n", &link_matrix[cCount], &read) == 1) {
	cur+=read;
	cCount++;
   }
   int rCount = 1;
   while(getline(&line, &count, fh)!=-1) {
	rCount++;
   }
   rewind(fh);
   printf("%d\n%d\n%d\n", rCount, cCount, net->nnodes);
   int i = 0;
   while(getline(&line, &count, fh)!=-1)
   {
     read = -1;
     cur  =  0;
     while(sscanf(line+cur, "%lf%n", &link_matrix[i], &read) == 1) {
	cur+=read;
	i=i+1;
     }
   }
   fclose(fh);
   sprintf(msg, "Link Matrix %d elements\n", cCount*rCount);
   debug_msg(msg);

   for (int j = 0; j < (net->nnodes)*(net->nnodes); j++) {
	if (link_matrix[j] > 0.) nlinks2 ++;
   }
   if (nlinks2 != net->nlinks) err_msg("Number of links in the matrix contradicts input file.\n");
   allocate_memory(net, link_matrix);
   init_arrays(net);
   sprintf(msg, "\nTotal links in network:  %d\n", net->nlinks);
   debug_msg(msg);
}

void allocate_memory(network *net, double *lm) {
  char msg[1024];
  node_ptr p_k, o_k;
  gpipe_ptr l_p;
  int l, il, ir;

  printf("Allocate Mem started\n");
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
    p_k->nr = 0;
    p_k->nl = 0;
    for (int j = 0; j < net->nnodes; j++) {
	if (lm[j*(net->nnodes)+i] > 0) (p_k->nr)++;
	if (lm[j*(net->nnodes)+i] < 0) (p_k->nl)++;
    }
    (p_k->adj_n) = (p_k->nl) + (p_k->nr);
    p_k->adj_k = malloc(sizeof(node_ptr)*(p_k->adj_n));
    p_k->adj_p = malloc(sizeof(gpipe_ptr)*(p_k->adj_n));    
    memset(p_k->adj_k, 0, (p_k->adj_n)*sizeof(node_ptr));
    memset(p_k->adj_p, 0, (p_k->adj_n)*sizeof(gpipe_ptr));
    p_k->outg = &(p_k->adj_p[0]);
    p_k->incm = &(p_k->adj_p[p_k->nr]);
    p_k->right = &(p_k->adj_k[0]);
    p_k->left = &(p_k->adj_k[p_k->nr]);
  }
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    p_k->key = i;
    il = 0; ir = 0;
    for( int j = 0; j < net->nnodes; j++) {
	o_k = net->knot[j];
	if (lm[j*(net->nnodes)+i] > 0.) {
	  p_k->right[ir] = o_k;
	  ir++;
	} else if (lm[j*(net->nnodes)+i] < 0.) {
	  p_k->left[il] = o_k;
	  il++;
	} 
    }

    sprintf(msg, "Knot %d address %p\n", i, p_k);
    debug_msg(msg);
    sprintf(msg, "Knot %d has %d outgoing and %d incoming connections\n", i, p_k->nr, p_k->nl);
    debug_msg(msg);
  }


  int idx, nlinks = 0;
  for (int i = 0; i < net->nnodes; i++) {
    p_k = net->knot[i];
    for (int j = 0; j < p_k->nr; j++) {
      o_k = p_k->right[j];
      p_k->outg[j] = malloc(sizeof(gpipe));
      net->link[nlinks] = p_k->outg[j];
      net->link[nlinks]->left = p_k;

      idx = 0;
      while (1) {
        if ((o_k->incm[idx])== NULL) {
	  o_k->incm[idx] = p_k->outg[j];
          net->link[nlinks]->right = o_k;
	  break;
	} else idx++;
      }   
      nlinks++;
    }
  }
  for (int i = 0; i < net->nlinks; i++) {
    l_p = net->link[i];
    l_p->key = i;
    l_p->c_id = NULL;
    sprintf(msg, "Link %2d (%p) left side connected to %p right side to %p\n", i, l_p, l_p->left, l_p->right);
    debug_msg(msg);
  }

  //err_msg("Complete");


  FILE *fh = fopen(net->nname,"r");
  char *dm;
  char msg2[80];

  debug_msg(net->nname);
  debug_msg("\n");
  while(1) {
    if(fgets(msg, 1024, fh)==NULL) break;
    if( strcmp(msg, "#lengths\n")==0) {
	dm = fgets(msg, 1024, fh);
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
	dm = fgets(msg, 1024, fh);
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
	dm = fgets(msg, 1024, fh);
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

  printf("Generating Compressor Data\n");
  fclose(fh);
  fh = fopen(net->fname,"r");
  for (int j = 0; j < net->nlinks+net->nnodes+1; j++) dm = fgets(msg, 1024, fh);
  net->cssr = malloc(sizeof(compressor)*net->ncomps);
  int key[net->ncomps], dummy;
  printf("Compressors: %d\n", net->ncomps);
  for (int j = 0; j < net->ncomps; j++) {
	dm = fgets(msg, 1024, fh);
	printf("Comp %d: %s -> %s", j, net->fname, msg);
	for (int l = 0; l < 5; l++) {
	  dummy = strtol(dm, &dm, 10);
	  if (l == 2) key[j] = dummy;
	}
  }
  for (int m = 0; m < net->ncomps; m++) {
    for (int j = 0; j < net->nlinks; j++) {
      l_p = net->link[j];
      if ((l_p->key) == key[m]) {
        l_p->c_id = &(net->cssr[m]);
	sprintf(msg2, "Added Compressor %p to link %d (%p)\n", &(net->cssr[m]), key[m], l_p);
	debug_msg(msg2);
      }
    }
  }
  printf("Compressor Data Passed\n");
}

void init_arrays(network *net) {
  gpipe_ptr p_k;

  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];	
    p_k->N = (p_k->L)*(net->npcent)-1;
    p_k->p = malloc((p_k->N)*sizeof(double));
    p_k->f = malloc((p_k->N)*sizeof(double));
    p_k->q = malloc((p_k->N)*sizeof(double));
    p_k->gamma = malloc((p_k->N)*sizeof(double));		
    p_k->fx = malloc((p_k->N)*sizeof(fftw_complex));
    //printf("Expected data points for pipe %d is %d\n", j, p_k->N);							
  }

}

void init_data(network *net) {
  char str[1024], msg[1024];
  FILE *fh;
  //printf("Init Data started\n");
  for (int j = 0; j < net->nlinks; j++) {
	sprintf(str, "%s_%03d.txt", net->dname, j);
	//debug_msg(str);
	fh = fopen(str, "r");
	if (fh == NULL) {
		sprintf(str, "Missing initial data along pipe %d, %d points expected\n", j, (net->link[j])->N  + 2);
		debug_msg(str);
		err_msg("Complete");
	} else {
		printf("Link %03d:\t", j);
		if (load_data(fh, net->link[j])) {
		   sprintf(msg, "Data for pipe %d has wrong number of lines.\n%d data points expected\n", j, (net->link[j])->N + 2);
		   //printf("%s", msg);
		   fclose(fh);
		   err_msg("Complete");
 		}
		printf("Interior points loaded. Loading Boundary ...\t");
		if ( j==0 ) {
			if (net->link[j]->c_id == NULL) {
				net->P0 = ovsq2*((net->link[j])->W_l + (net->link[j])->Wb_l);
				net->F0 = (net->link[0])->Fl;
			} else {
				net->P0 = ovsq2*((net->link[j])->W_l + (net->link[j])->Wb_l)/(net->link[j]->c_id->cratio);	
				net->F0 = (net->link[0])->Fl;
			}
		}
		fclose(fh);
		printf("Done\n");
	}
  }

  printf("Load Data Complete: Left Flux is %.12e\n", net->F0/cs);
  //net->P0 = (net->link[0])->Pl;
  //rescale_data(net);
  //printf("Init Data passed\n");
  //exit(1);
}


int load_data(FILE *fh, gpipe_ptr lnk) {
  char line[512], val0[64], val1[64], val2[64], val3[64];
  char msg[80];
  char *dm;
  dm = fgets(line, 512, fh);
  dm = fgets(line, 512, fh);
  dm = fgets(line, 512, fh);
  sprintf(msg, "%s", dm);
  int k = 0;
  if (fgets(line, 512, fh) != NULL) {
    sscanf(line,"%s\t%s\t%s\t%s\n", val0, val1, val2, val3);
    //printf("%s\t %p\n", line, lnk->c_id);

    //(lnk->left)->F = strtod(val2, NULL)*cs;  // dumb


    (lnk->W_l) = ovsq2*(1e6*strtod(val1, NULL) + strtod(val2, NULL)*cs);
    (lnk->Wb_l) = ovsq2*(1e6*strtod(val1, NULL) - strtod(val2, NULL)*cs);

    if (lnk->c_id != NULL) {
       (lnk->left)->P = 1e6*strtod(val1, NULL)/(lnk->c_id)->cratio;
       (lnk->W_l) = ovsq2*(1e6*strtod(val1, NULL) + strtod(val2, NULL)*cs);
       (lnk->Wb_l) = ovsq2*(1e6*strtod(val1, NULL) - strtod(val2, NULL)*cs);
       //printf("c = %.12e\nP = %.12e\n",(lnk->c_id)->cratio, (lnk->left)->P);  
    } else {
       (lnk->left)->P = 1e6*strtod(val1, NULL);
       (lnk->W_l) = ovsq2*(1e6*strtod(val1, NULL) + strtod(val2, NULL)*cs);
       (lnk->Wb_l) = ovsq2*(1e6*strtod(val1, NULL) - strtod(val2, NULL)*cs);
    }
    (lnk->left)->G = pow(strtod(val2, NULL)*cs/(1e6*strtod(val1, NULL)), 2);
#if ZERO_GAMMA
    (lnk->left)->G = 0.;
#endif
    lnk->Fl = strtod(val2, NULL)*cs;	     // better
    (lnk->left)->Q = strtod(val3, NULL)*cs;
    //printf("%.15e\n",lnk->left->F);
    k = 1;
  } else err_msg("Error reading input.");
  while (fgets(line, 512, fh) != NULL) {
	sscanf(line,"%s\t%s\t%s\t%s\n", val0, val1, val2, val3);
	k++;
	if (k == (lnk->N+2)) {
           (lnk->right)->P = 1e6*strtod(val1, NULL);
           lnk->Fr = strtod(val2, NULL)*cs;
           (lnk->right)->Q = strtod(val3, NULL)*cs;
	   printf("element %d of %d", k-2, lnk->N);
	   (lnk->right)->G = pow(strtod(val2, NULL)*cs/(1e6*strtod(val1, NULL)), 2);
#if ZERO_GAMMA 
	   (lnk->right)->G = 0.;
#endif
	   (lnk->W_r) = ovsq2*(1e6*strtod(val1, NULL) + strtod(val2, NULL)*cs);
    	   (lnk->Wb_r) = ovsq2*(1e6*strtod(val1, NULL) - strtod(val2, NULL)*cs);
	   sprintf(msg, "Read %d(%d) lines\n", k, lnk->N+2);
	   debug_msg(msg);
	   break;
	}
	lnk->p[k-2] = 1e6*strtod(val1, NULL);
	lnk->f[k-2] = strtod(val2, NULL)*cs;
	lnk->gamma[k-2] = pow(lnk->f[k-2]/lnk->p[k-2], 2);  // 
	lnk->q[k-2] = strtod(val3, NULL)*cs;
	//printf("%.15e\n",lnk->f[k-2]);
        //if (k == 5) exit(1);
  }
  (lnk->Pl) = (lnk->left)->P;

  if (lnk->c_id != NULL) {
       lnk->Pl = ((lnk->left)->P)*(lnk->c_id)->cratio;
       printf("P = %.12e\t", lnk->Pl);  
       printf("c = %.12e\t",(lnk->c_id)->cratio);  
  } else printf("P = %.12e\t", lnk->Pl);  
  (lnk->Pr) = (lnk->right)->P;

  /*printf("Left: W = %e\tWb = %e\n", lnk->W_l, lnk->Wb_l);
  for (int j = 0; j < lnk->N; j++) printf("%2d of %d: W = %e\tWb = %e\n", j, lnk->N, ovsq2*(lnk->p[j] + lnk->f[j]), ovsq2*(lnk->p[j] + lnk->f[j]));
  printf("Right: W = %e\tWb = %e\n", lnk->W_r, lnk->Wb_r);
  err_msg("");*/
  printf("Read %d lines of %d expected from initial data, Check %d \n", k, lnk->N+2, k == lnk->N+2);
  if (k == lnk->N+2) {
	return 0;
  } else {
	return 1;
  }
#if ZERO_GAMMA
  memset(lnk->gamma, 0, (lnk->N)*sizeof(double)); 
#endif
}

/*void rescale_data(network *net) {
  gpipe_ptr p_k;
  for (int j = 0; j < net->nlinks; j++) {
    p_k = net->link[j];
    for(int k = 0; k < p_k->N; k++) {
      p_k->p[k] = 1000.*(p_k->p[k]);
      p_k->f[k] = (p_k->f[k])*(net->c);
      p_k->q[k] = p_k->q[k];
    }
  }
}*/

void save_data(FILE* fh, gpipe_ptr lnk, network* net, double time) {
  if (fh == NULL) err_msg("Cannot write to file.");
  fprintf(fh, "# 1. x, km 2. P, MPa 3. Flux, kg/m^2/s\n# N = %d\ttime = %.12e\n\n", lnk->N+2, time);
#if 1  
  if (lnk->c_id == NULL) {
    fprintf(fh, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 0., 1e-6*(lnk->left->P), (lnk->Fl)/net->c, lnk->q[0]/(net->c), creal(lnk->fx[0])/(net->c)); //ovsq2*(lnk->W_l - lnk->Wb_l)/(net->c));  ovsq2*(lnk->W_l + lnk->Wb_l)/1000000.
    //printf( "%.15e\t%.15e\t%.15e\n", 0., ovsq2*(lnk->W_l + lnk->Wb_l)/1000., ovsq2*(lnk->W_l - lnk->Wb_l)/(net->c)); 
    //exit(1);
  } else {
    fprintf(fh, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 0., 1e-6*(lnk->c_id)->cratio*(lnk->left->P), (lnk->Fl)/net->c, lnk->q[0]/(net->c), creal(lnk->fx[0])/(net->c)); //ovsq2*((lnk->c_id)->cratio*(lnk->W_l + lnk->Wb_l))/1000000., ovsq2*(lnk->W_l - lnk->Wb_l)/(net->c));
    //printf("%.15e\t%.15e\t%.15e\n", 0., 1e-6*(lnk->c_id)->cratio*(lnk->left->P), (lnk->Fl)/net->c); //ovsq2*((lnk->c_id)->cratio*(lnk->W_l + lnk->Wb_l))/1000000., ovsq2*(lnk->W_l - lnk->Wb_l)/(net->c));
    //exit(1);
  }
#endif  
  for (int j = 0; j < lnk->N; j++) {
    fprintf(fh, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", (j+1)*(lnk->L)/((lnk->N)+1), (lnk->p[j])*1e-6, (lnk->f[j])/(net->c), (lnk->q[j])/(net->c), creal(lnk->fx[j])/(net->c) );
  }
  fprintf(fh, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", lnk->L,  1e-6*(lnk->right->P), (lnk->Fr)/net->c, lnk->q[lnk->N-1]/(net->c), creal(lnk->fx[lnk->N-1])/(net->c)); //ovsq2*(lnk->W_r + lnk->Wb_r)/1000000., ovsq2*(lnk->W_r - lnk->Wb_r)/(net->c)); 
}


void init_save_balance(network* net) {
  FILE *fh;
  char filename[80];
  for (int j = 0; j < net->nlinks; j++){
    sprintf(filename, "balance_%03d.txt", j);
    fh = fopen(filename, "w");
    fprintf(fh, "# 1. x, km 2. Dphi/Dt 3. D p/Dx 4. D(phi^2/rho)/ Dx\n\n");
    fclose(fh);
  }
}

void save_balance(network* net, network* netb, double dt, double dx) {
  gpipe_ptr tmp1, tmp2;
  double f1, f2, dp1, dp2, d1, d2;
  double T1, T2, T3;
  FILE *fh;
  char filename[80];
  for (int j = 0; j < net->nlinks; j++){
    tmp1 = net->link[j];
    tmp2 = netb->link[j];
    sprintf(filename, "balance_%03d.txt", j);
    fh = fopen(filename, "a");
    T1 = 0.;    T2 = 0.;    T3 = 0.;
    for (int k = 1; k < tmp1->N-1; k++) {
      f1 = tmp2->f[k]; dp1 = (tmp2->p[k+1]-tmp2->p[k-1])/dx; d1 = pow(tmp2->f[k+1],2)/(tmp2->p[k+1]) - pow(tmp2->f[k-1],2)/(tmp2->p[k-1]);
      f2 = tmp1->f[k]; dp2 = (tmp1->p[k+1]-tmp1->p[k-1])/dx; d2 = pow(tmp1->f[k+1],2)/(tmp1->p[k+1]) - pow(tmp1->f[k-1],2)/(tmp1->p[k-1]);
      T1 += pow((f2-f1)/dx,2);
      T2 += pow(0.5*(dp1+dp2),2);
      T3 += pow(0.5*(d1+d2)/dx,2);
    }
    fprintf(fh, "%.15e\t%.15e\t%.15e\t%.15e\n", (net->curr_T) + 0.5*dt, 0.001*sqrt(T1)*dx/(tmp1->L),  0.001*sqrt(T2)*dx/(tmp1->L), 0.001*sqrt(T3)*dx/(tmp1->L));
    fclose(fh);
  } 


}










