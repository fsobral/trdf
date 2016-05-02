#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void c_calobjf(int n, double *x, double *f, int *flag) {

  // Problem HS37 from Hock-Schittkowski collection

  *flag = 0;

  *f = - x[0] * x[1] * x[2];

  return;

}

void c_calcon(int n, double *x, int ind, double *c, int *flag) {

  // Problem HS37 from Hock-Schittkowski collection

  *flag = 0;

  if ( ind == 0 ) {

    *c = - 72.0 + x[0] + 2.0 * x[1] + 2.0 * x[2];

  }
  else if ( ind == 1 ) {

    *c = - x[0] - 2.0 * x[1] - 2.0 * x[2];

  }
  else {

    *flag = -1;

  }

  return;

}

void c_caljac(int n, double *x, int ind, int *jcvar, double *jcval,
	    int *jcnnz, int lim, _Bool *lmem, int *flag) {

  // Problem HS37 from Hock-Schittkowski collection

  *flag = 0;

  *lmem = 0;

  if ( ind == 0 ) {

    *jcnnz = 3;
    jcvar[0] = 0;
    jcval[0] = 1.0;
    jcvar[1] = 1;
    jcval[1] = 2.0;
    jcvar[2] = 2;
    jcval[2] = 2.0;

  } else if ( ind == 1 ) {

    *jcnnz = 3;
    jcvar[0] = 0;
    jcval[0] = - 1.0;
    jcvar[1] = 1;
    jcval[1] = - 2.0;
    jcvar[2] = 2;
    jcval[2] = - 2.0;

  } else {
    
    *flag = -1;

  }

  return;

}

void c_calhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
	   int *hcnnz, int lim, _Bool *lmem, int *flag) {
  

  // Problem HS37 from Hock-Schittkowski collection

   *flag = 0;

   *lmem = 0;

   *hcnnz = 0;

   return;

}


int main()
{

  int i, fcnt, m, n;

  double f, feas, *x, *l, *u;

  _Bool *equatn, *linear, ccoded[2];

  // Problem HS37 from Hock-Schittkowski collection

  // Number of variables
  n = 3;

  // Allocates initial point and bounds

  x = (double*) malloc(n * sizeof(double));
  l = (double*) malloc(n * sizeof(double));
  u = (double*) malloc(n * sizeof(double));
  
  if ( x == NULL || l == NULL || u == NULL ) {

    printf("\nMemory allocation error.\n");

    return 1;

  }

  for (i = 0; i < n; i++) {

    l[i] = 0.0;
    u[i] = 42.0;
    x[i] = 10.0;

  }

  // Number of constraints
  m = 2;

  equatn = (_Bool*) malloc(m * sizeof(_Bool));
  linear = (_Bool*) malloc(m * sizeof(_Bool));

  if ( equatn == NULL || linear == NULL ) {

    printf("\nMemory allocation error.\n");

    return 1;

  }

  // Type of constraints

  for (i = 0; i < m; i++) equatn[i] = 0;
  for (i = 0; i < m; i++) linear[i] = 1;

  // Coded subroutines for constraints' derivatives

  ccoded[0] = 1;
  ccoded[1] = 1;


  // Calls the solver
  easytrdf(n,x,l,u,m,equatn,linear,ccoded,&c_calobjf,&c_calcon,
           &c_caljac,&c_calhc,&f,&feas,&fcnt);

  // Frees memory
  free(x);
  free(l);
  free(u);
  free(equatn);
  free(linear);

}
