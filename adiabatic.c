#include "network.h"

static double *Z, *sqrZ, *Q, *tmp1, *tmp2, *tmp3, *Dphi; 
static double *tmp4;
static double *phi_A, P_aux, mean_xi; 
static double *P_A, T;
static double gam, cs, dx, L, alpha;
static gpipe_ptr p_k;
static int N, savenum;

#define CACHE_LINE_IN_DOUBLES	8

void save_array(char* fname, double *in) {
  FILE *fh = fopen(fname, "w");
  for (int j = 0; j < N+2; j++) {
    fprintf(fh, "%.15e\t%.15e\n", 1e-3*j*dx, in[j]);
  } 
  fclose(fh);
}


void integrate_array_backward (double *input, double *output) {
	unsigned long int i, j, lines_number = (unsigned long int)((1.0*(N+1)-1)/(CACHE_LINE_IN_DOUBLES)), remainder = (N+2) - lines_number*(CACHE_LINE_IN_DOUBLES);
	double half_DX = 0.5*dx, test_value1;//, test_value2;
	double *temp_ptr1, *temp_ptr2;

	output [(N+2)-1] = 0.0;
	/* Handling all whole cache lines */
	for (j = 0; j < lines_number; ++j) {
		temp_ptr1 = &(input [(N+2)-1 - j*(CACHE_LINE_IN_DOUBLES)]);
		temp_ptr2 = &(output [(N+2)-1 - j*(CACHE_LINE_IN_DOUBLES)]);
		
		test_value1 = *(temp_ptr1 - (CACHE_LINE_IN_DOUBLES));
		//test_value2 = *(temp_ptr2 - (CACHE_LINE_IN_DOUBLES));

		for (i = 1; i < (CACHE_LINE_IN_DOUBLES); ++i) {
			*(temp_ptr2 - i) =  *(temp_ptr2 - i + 1) + (*(temp_ptr1 - i) + *(temp_ptr1 - i + 1))*half_DX;
		}
		*(temp_ptr2 - (CACHE_LINE_IN_DOUBLES)) =  *(temp_ptr2 - (CACHE_LINE_IN_DOUBLES) + 1) + (test_value1 + *(temp_ptr1 - (CACHE_LINE_IN_DOUBLES) + 1))*half_DX;
	}
	/* Handling remainder of the array */
	test_value1 = input[0];
	//test_value2 = output[0];
	
	temp_ptr1 = (input + remainder);
	temp_ptr2 = (output + remainder);
	for (i = 1; i < remainder; ++i) {
		*(temp_ptr2 - i) =  *(temp_ptr2 - i + 1) + (*(temp_ptr1 - i) + *(temp_ptr1 - i + 1))*half_DX;
	}
	*(temp_ptr2 - remainder) =  *(temp_ptr2 - remainder + 1) + (test_value1 + *(temp_ptr1 - remainder + 1))*half_DX;
}

void integrate_array (double *input, double *output) {
	unsigned long int i;
	double half_DX = 0.5*dx;

	output [0] = 0.0;
	for (i = 1; i < (N+2); ++i) {
		output[i] = output [i-1] + (input[i] + input[i-1])*half_DX;
	}
}

double compute_integral_noise (fftw_complex *input) {
	unsigned long int i;
	double sum = 0.0;

	for (i = 1; i < N-1; ++i) {
		sum += creal(input[i]);
	}

	return creal(dx*(0.5*(input[0] + input[N-1]) + sum));
}

double compute_integral (double *input) {
	unsigned long int i;
	double sum = 0.0;

	for (i = 1; i < N+1; ++i) {
		sum += input[i];
	}

	return dx*(0.5*(input[0] + input[N+1]) + sum);
}


void calc_pressure() {
  for (int j = 0; j < N+2; j++) {
    tmp1[j] = alpha*(phi_A[j] + 0.*Dphi[j])*cs*cs*fabs(phi_A[j] + 0.*Dphi[j])/(Z[j]*P_A[N+1]*P_A[N+1]);  // 1 is new adiabatic, 0 is old
  }

  //save_array("Z.txt", Z);
  //save_array("tmp1.txt", tmp1);


  P_aux = P_A[N+1];
  integrate_array_backward(tmp1, P_A);
  //save_array("secondtermI.txt", P_A);

  P_A[N+1] = P_aux;
  //save_array("sqrZ.txt", sqrZ);  
  for (int j = 0; j < N+1; j++) {
    P_A[j] = P_aux*sqrZ[j]*sqrt(1 + P_A[j] ); 
  }
  P_A[0] = 2.*P_A[1] - P_A[2];
  //save_array("pressure.txt", P_A);  
  //exit(1);
  //printf("%.15e\n", P_A[0]*cs*cs/(P_aux*P_aux));
}

double rhs(fftw_complex *in) {
  mean_xi = compute_integral_noise(in);

  for (int j = 0; j < N+2; j++) {
    tmp2[j] = Z[j]/P_A[j];
    tmp3[j] = creal(in[j]);
  }
  double mean_Zp = compute_integral(tmp2);
  integrate_array(tmp3, tmp4);	
  integrate_array(tmp2, Dphi);
  for (int j = 0; j < N+2; j++) {
  	Dphi[j] = (tmp4[j] - 1.*Dphi[j]*mean_xi/mean_Zp)/cs;
  }   
  return 2.*cs*mean_xi/mean_Zp;
}



void adiabatic_rk4(fftw_complex *in, double dt) {
  double k1, k2, k3, k4, p1, p2, p3, p4;
  double P0 = P_A[N+1]*P_A[N+1];
  //save_array("P_A.txt", P_A);
  k1 = rhs(in);
  
  //printf("%.12e\n", k1); //exit(1); 
  p1 = P0 + 0.5*dt*k1;
  P_A[N+1] = sqrt(p1);

  calc_pressure();
  //printf("P0 = %.12e\tP_right = %.12e\n", P_A[0], P_A[N+1]);
  k2 = rhs(in);
  p2 = P0 + 0.5*dt*k2;
  P_A[N+1] = sqrt(p2);
  calc_pressure();

  k3 = rhs(in);
  p3 = P0 + dt*k3;
  P_A[N+1] = sqrt(p3);
  calc_pressure();

  k4 = rhs(in);
  p4 = P0 + dt*(k1 + 2.*k2 + 2.*k3 + k4)/6.;
  P_A[N+1] = sqrt(p4);
  
  k1 = rhs(in);
  calc_pressure();
  //printf("k1 = %.12e\tk2 = %.12e\tk3 = %.12e\tk4 = %.12e\nP0 = %.12e\tP_right = %.12e\n", k1, k2, k3, k4, P0, P_A[N+1]);
  //printf("P0 = %.12e\tP_right = %.12e\n", P_A[0], P_A[N+1]);
  //save_adiabatic(); 
  //save_array("P_A2.txt", P_A);
  //exit(1);
  //if ( !(savenum % 10)) save_adiabatic();
  T = T + dt;
}



void save_adiabatic() {
  char filename[256];
  
  sprintf(filename, "./pipe_000/adiabatic_%03d.txt", savenum);
  savenum++;
  FILE *fh = fopen(filename, "w");
  for (int j = 0; j < N+2; j++) fprintf(fh, "%.15e\t%.15e\t%.15e\n", 0.001*j*dx, P_A[j]*1.e-6, phi_A[j]+1.*Dphi[j]);
  fclose(fh);
}

void save_adiabatic_temporal() {
  
  FILE *fh = fopen("adiabatic_pressure.txt", "a");
  fprintf(fh, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n", T, 
	P_A[      0]*1.e-6, P_A[(N+1)/8]*1.e-6, P_A[(N+1)/4]*1.e-6, P_A[3*(N+1)/8]*1.e-6,
	P_A[(N+1)/2]*1.e-6, P_A[5*(N+1)/8]*1.e-6, P_A[3*(N+1)/4]*1.e-6, P_A[7*(N+1)/8]*1.e-6,	
        P_A[N+1]*1.e-6);
  fclose(fh);

  fh = fopen("adiabatic_flux.txt", "a");
  fprintf(fh, "%.12e\t%.12e\t%.12e\t%.12e\n", T, phi_A[0]+Dphi[0], phi_A[(N+1)/2]+Dphi[(N+1)/2], phi_A[N+1]+Dphi[N+1]); 
  fclose(fh);
}


void prepare_arrays_adiabatic(network_ptr net){
  p_k = net->link[0];
  double phi0;
  N = (p_k->N); cs = net->c;  alpha = (net->diss)/(p_k->d);
  L = 1000.*(p_k->L); 
  dx = L/(N+1);
  savenum = 1;
  T = 0.;

  FILE *fh;
  fh = fopen("adiabatic_pressure.txt", "w");
  fprintf(fh, "# 1. t 2. Inlet 3. Center 4. Outlet\n\n");
  fclose(fh);
  fh = fopen("adiabatic_flux.txt", "w");
  fprintf(fh, "# 1. t 2. Inlet 3. Center 4. Outlet\n\n");
  fclose(fh);

  tmp1 = malloc((N+2)*sizeof(double));
  tmp2 = malloc((N+2)*sizeof(double));
  tmp3 = malloc((N+2)*sizeof(double));
  tmp4 = malloc((N+2)*sizeof(double));
  Z = malloc((N+2)*sizeof(double));
  sqrZ = malloc((N+2)*sizeof(double));
  Q = malloc((N+2)*sizeof(double));
  phi_A = malloc((N+2)*sizeof(double));
  Dphi = malloc((N+2)*sizeof(double));
  P_A = malloc((N+2)*sizeof(double));
  P_A[N+1] = (p_k->right)->P;
  memset(Dphi, 0, (N+2)*sizeof(double));

  double S = 0.25*pi*(p_k->d)*(p_k->d); // below divide by S
  Q[0] = (p_k->left->Q);
  for (int j = 1; j < N+1; j++) Q[j] = p_k->q[j-1];
  Q[N+1] = p_k->right->Q;
  
  for (int j = 1; j < N+1; j++) {
	tmp1[j] = 0.5*alpha*p_k->gamma[j-1];
  }
  tmp1[0] = 0.5*alpha*p_k->gamma[0];
  tmp1[N+1] = 0.5*alpha*p_k->gamma[N-1];
  //save_array("gamma.txt", tmp1);
  
  // generate Z 
  double tmp = 0.;
  integrate_array_backward (tmp1, Z);
  for (int j = 0; j < N+2; j++) {
	tmp = Z[j];
	Z[j]    = exp(-2.0*tmp);
	sqrZ[j] = exp(-1.0*tmp);
  }

  // generate stationary flux
  integrate_array (Q, phi_A); 
  //save_array("integralofinjection.txt", phi_A);
  phi0 = (p_k->W_l - p_k->Wb_l)/sqrt(2)/cs; 
  printf("phi0 = %.12e\n", phi0);
  for (int j = 0; j < N+2; j++) phi_A[j] = phi_A[j]/cs + phi0;

  // generate stationary pressure

  /*for (int j = 0; j < N+2; j++) tmp1[j] = alpha*phi_A[j]*fabs(phi_A[j])/Z[j];
  integrate_array_backward(tmp1, P_A);
  for (int j = 0; j < N+2; j++) P_A[j] = ((p_k->right)->P)*sqrt( Z[j]*(1 +  P_A[j]*cs*cs/pow((p_k->right)->P,2) )); */
  calc_pressure();
  //save_array("injection.txt", Q);
  //save_array("flux.txt", phi_A);
  //save_array("Z.txt", Z);  
  //exit(1);

  /*FILE *fh = fopen("verify.txt", "w");
  for (int j = 0; j < N+2; j++) {
    fprintf(fh, "%.12e\t%.12e\t%.12e\n", 0.001*j*dx, P_A[j]*1.e-6, phi_A[j]);
  }
  fclose(fh);*/
  
}







/*void adiabatic() {

	compute_Z();
*/


