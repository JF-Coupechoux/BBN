#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <locale.h>
#include "nrutil.c"
#include "nrutil.h"

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 1000000000
#define TINY 1.0e-30
#define SORTIE "M5.out"

/*---------------------------------------------------------------------------------------------------------------*/

#define N_VARIABLES 3

#define pi 3.14159265359
#define G 6.67259e-11
#define c 2.997924562e8
#define h_ 1.05457266913e-34
#define e 1.60217733e-19
#define k 1.380658e-23     /* Quelques constantes. */

int test1=0;
double rho_0_c, m, omega_rad, omega_b, omega_phi, domega_dx, omega0_phi, omega0_b, omega0_rad, H0,omega_rad_adapt, omega_b_adapt, omega_phi_adapt, Q,omega_de;

/*---------------------------------------------------------------------------------------------------------------*/

double V_pot(double x)
{
	return 0.5*m*m*x*x;
}


double dV_pot(double x)
{
	return m*m*x;
}



void fonction_diff(double t, double val[], double dval_dx[])

{ 
        double V, dV_dsigma;
	
	V= V_pot(val[2]);
	dV_dsigma= dV_pot(val[2]);

	/*
	val[1] = a/a0 */
	  
	dval_dx[1] = H0 * val[1] * sqrt(omega0_rad * pow(val[1],-4) + omega0_b
	* pow(val[1],-3)+ omega_de + omega0_phi*( 0.5 * val[3]*val[3] + V));
	
	/* val[2] = sigma
	val[3] = dval[2]/dx */
	
	dval_dx[2] = val[3];
	
	dval_dx[3] = - 3. * dval_dx[1]/val[1] * val[3] - dV_dsigma;
}

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
double yerr[], void (*derivs)(double, double [], double []))
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
  ak2=dvector(1,n);
  ak3=dvector(1,n);
  ak4=dvector(1,n);
  ak5=dvector(1,n);
  ak6=dvector(1,n);
  ytemp=dvector(1,n);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4); 
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5); 
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6); 
  for (i=1;i<=n;i++) 
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=1;i<=n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  free_dvector(ytemp,1,n);
  free_dvector(ak6,1,n);
  free_dvector(ak5,1,n);
  free_dvector(ak4,1,n);
  free_dvector(ak3,1,n);
  free_dvector(ak2,1,n);
}

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], 
double *hdid, double *hnext, void (*derivs)(double, double [], double []))
{
  int i;
  double errmax,h,htemp,xnew,*yerr,*ytemp;
  yerr=dvector(1,n);
  ytemp=dvector(1,n);
  h=htry;
  for (;;) {
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs); 
    errmax=0.0; 
    for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps; 
    if (errmax <= 1.0) break; 
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
    xnew=(*x)+h;
    if (xnew == *x) nrerror("stepsize underflow in rkqs");
  }
  if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=1;i<=n;i++) y[i]=ytemp[i];
  free_dvector(ytemp,1,n);
  free_dvector(yerr,1,n);
}


double calcul(double y[], double x1, double x2, int nvar, void (*derivs)(double, double [], double []))
{
  time_t T1;
  T1 = time(0);
  int nstp,i;
  double x,hnext,hdid,h,eps,h1,hmin,atemp=y[1];
  double *yscal, *dydx, *yout;
  double rho_phi, omega_rad, omega_phi, omega_b;
  FILE *output;
  yscal=dvector(1,nvar);
  dydx=dvector(1,nvar);
  yout=dvector(1,nvar);
  eps=1.0e-10;
  h1=1.e-25;
  hmin=0.0;
  x=x1;
  if ((x2-x1)>0) h=h1; else h=-h1;
  for (nstp=1;nstp<=MAXSTP;nstp++) { 
    nstp=1;
    (*derivs)(x,y,dydx);
    if(test1==0) {
      output=fopen(SORTIE,"w");
      fprintf(output,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", x, y[1], y[2], y[3], m*h_/e,Q,0.,0.,0.);
      fclose(output);
      test1=1;
    }
    else {
      if(fabs((atemp-y[1])/atemp)>0.001) {
	atemp=y[1];
	rho_phi= ( 0.5 * y[3]*y[3] + V_pot(y[2]));
	omega_rad = omega0_rad * pow(y[1],-4.);
	omega_b = omega0_b * pow(y[1],-3.);
	omega_phi= omega0_phi *rho_phi;
	output=fopen(SORTIE,"a");
	fprintf(output,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", x, y[1], y[2], y[3], omega_rad, omega_b, omega_phi, omega_de, omega_b/omega_phi);
	fclose(output);
      }
    }
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (difftime(time(0),T1)>3)
      { omega0_phi=omega_phi*pow(y[1],3.);
	while(y[1]<1) {
		y[1]=y[1]*pow(10.,0.05);
		omega_rad = omega0_rad * pow(y[1],-4.);
		omega_b = omega0_b * pow(y[1],-3.);
		omega_phi= omega0_phi * pow(y[1],-3.);
		output=fopen(SORTIE,"a");
		fprintf(output,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", x, y[1], y[2], y[3], omega_rad, omega_b, omega_phi, omega_de,omega_b/omega_phi);
		fclose(output);
	}
	y[1]=1.;
	omega_rad = omega0_rad * pow(y[1],-4.);
	omega_b = omega0_b * pow(y[1],-3.);
	omega_phi= omega0_phi * pow(y[1],-3.);
	output=fopen(SORTIE,"a");
	fprintf(output,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", x, y[1], y[2], y[3], omega_rad, omega_b, omega_phi,omega_de, omega_b/omega_phi);
	fclose(output);	
	free_dvector(dydx,1,nvar);
	free_dvector(yscal,1,nvar);
	free_dvector(yout,1,nvar);
	return 1;
      }
    rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
    h=hnext;
  }

  printf("ca n'a pas ete jusqu'a x2 \n");
  free_dvector(dydx,1,nvar);
  free_dvector(yscal,1,nvar);
  free_dvector(yout,1,nvar);
  return 1;
}



/*---------------------------------------------------------------------------------------------------------------*/
int main()

{
  double h,T0,x0,x1,val [N_VARIABLES+1],pas,erreur,junk;
  int test;
  FILE *ext;
  
  
  	x0=1.e-25 ;                /* valeur initiale de t */
  	val[1]= 1.e-12; 	   /* valeur initiale du parametre d'expansion */  
  	val[2]= -1e7;  	   /* valeur initiale de sigma */
  	val[3]= 1.e+10; 		   /* valeur initiale de dsigma/dt */
  m=1.e-24 * e/h_;      /* voir article 2 */
  Q=0.*6.96e-18;              /* tuning pour trouver omega_b/omega_phi = 0.2 */ 


  erreur =1.e-10;
  
  h=0.7;         /* Constante de Hubble. */
  T0=2.728 ;      /* Parametre cosmologique. */

  
  H0= h * 100. / 3.0856e19;
  
  rho_0_c =  3.0 * H0 * H0 / (8.0 * pi * G);
  
  omega0_rad = (1+3.*7./8.*pow(4./11.,4./3.))*pi * pi / 15. * pow(k * T0, 4.0) / pow(h_,3.0)/pow(c,5.0) / rho_0_c ;
  omega0_b = 0.024/(h*h);
  omega0_phi = 1. / rho_0_c ;
  omega_de= 0.69;   /* on veut omega_phi^0 = 0.96 */
  
  printf("a_eq_exp=%.10e \n",4.15e-5/(0.3*h*h));
  x1=100.e9*3600*24*365;     /* age de l'Univers < 100 milliards d'annees */
  
  calcul(val,x0,x1,3,fonction_diff);



  FILE * f;
  f = popen("gnuplot", "w");
  fprintf(f, "load \"trace_m.gnu\" \n");
  fflush(f);
  sleep(5000);
  pclose(f);
  return 0;
}
