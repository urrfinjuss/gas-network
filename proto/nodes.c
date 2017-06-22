#include "header.h"

static char text[80];
static unsigned int step;
static unsigned int ntype;
static long int N;
static FTYPE boost, wbr, w;
static FTYPE FlowI, FlowJ;
static FTYPE SectI, SectJ;
static FTYPE node_flow;
static gcomp_ptr cmp;
static gpipe_ptr left, right;

void forward_nodes() {
	dmesg("forward_nodes:\n", 0);
	step = par.curr;
	for (int j = 0; j < net.nodes; j++) {
		ntype = net.node[j].type;
		if (ntype == 1) forward_pressure_node(&net.node[j]);
		else            forward_transport_node(&net.node[j]); 
	}
	dmesg("forward_nodes:\tPassed\n", 0);
}

void forward_pressure_node(gnode_ptr in) {
	printf("Processing Pressure Node %p\n", in);
	in->p[1] = in->var[step+1];
	forward_incoming_pipes(in);
	forward_outgoing_pipes(in);
}

void forward_incoming_pipes(gnode_ptr in) {
	printf("Loop over incoming pipes (%u)\n", in->nleft);
	for (int k = 0; k < in->nleft; k++) {
		left = in->left[k];
		N = (left->N) - 1;
		w = left->y[0][N] + left->y[1][N];
		left->fd[1] = w - in->p[1];
		printf("Init Pressure/Flux:\t%.12e\t%.12e\n", in->p[0], left->fd[0]/par.sound);
		printf("Scnd Pressure/Flux:\t%.12e\t%.12e\n", in->p[1], left->fd[1]/par.sound);
	}
}

void forward_outgoing_pipes(gnode_ptr in) {
	printf("Loop over outgoing pipes (%u)\n", in->nright);
	for (int k = 0; k < in->nright; k++) {
		right = in->right[k];
		cmp = right->cmp;
		wbr = right->y[0][0] - right->y[1][0];
		if (cmp) boost = cmp->boost[step+1];
		else boost = 1.0;
		right->fs[1] = boost*in->p[1]-wbr;
		printf("Init Pressure/Flux:\t%.12e\t%.12e\n", in->p[0], right->fs[0]/par.sound);
		printf("Scnd Pressure/Flux:\t%.12e\t%.12e\n", in->p[1], right->fs[1]/par.sound);
	}
}

void forward_transport_node(gnode_ptr in) {
	calculate_pressure(in);
	forward_incoming_pipes(in);
	forward_outgoing_pipes(in);
}


void calculate_pressure(gnode_ptr in) {
	SectI = 0.0;
	FlowI = 0.0;
	for (int k = 0; k < in->nleft; k++) {
		left = in->left[k];
		N = (left->N) - 1;
		w = left->y[0][N] + left->y[1][N];
		FlowI += 0.25*pi*pow(left->wid,2)*w;
		SectI += 0.25*pi*pow(left->wid,2);
	}
	SectJ = 0.0;
	FlowJ = 0.0;
	for (int k = 0; k < in->nright; k++) {
		right = in->right[k];
		cmp = right->cmp;
		wbr = right->y[0][0] - right->y[1][0];
		if (cmp) boost = cmp->boost[step+1];
		else boost = 1.0;
		FlowJ += 0.25*pi*pow(right->wid,2)*wbr;
		SectJ += 0.25*pi*pow(right->wid,2)*boost;
	}
	node_flow = in->var[step+1];
	in->p[1] = (FlowI + FlowJ) - par.sound*node_flow;
	in->p[1] = in->p[1]/(SectI + SectJ);
	//printf("Transport Node Pressure: %.12e\n", in->p[1]);
}










