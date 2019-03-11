#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <locale.h>
#include "nrutil.c"
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define RM 0.61803399 
#define C (1.0-RM)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 1000000000
#define TINY 1.0e-30

/*---------------------------------------------------------------------------------------------------------------*/

#define N_VARIABLES 3

#define pi 3.14159265359
#define G 6.67259e-11
#define c 2.997924562e8
#define h_ 1.05457266913e-34
#define e 1.60217733e-19
#define k 1.380658e-23     /* Quelques constantes. */

int test1=0;
double rho_0_c, m, omega_rad, omega_b, omega_phi, domega_dx, omega0_phi, omega0_b, omega0_rad, H0,omega_rad_adapt, omega_b_adapt, omega_phi_adapt, Q, a_eq_exp, cut, alpha, beta, lambda,V0;

/*---------------------------------------------------------------------------------------------------------------*/
double V_pot(double x)
{
	return V0+0.5*m*m*pow(x,2)+0.25*lambda*pow(x,4);
}


double dV_pot(double x)
{
	return m*m*x+lambda*pow(x,3);
}



void fonction_diff(double t, double val[], double dval_dx[])

{ 
        double V, dV_dsigma;
	
	V= V_pot(val[2]);
	dV_dsigma= dV_pot(val[2]);

	/*
	val[1] = a/a0 */
	  
	dval_dx[1] = H0 * val[1] * sqrt(omega0_rad * pow(val[1],-4) +  omega0_b
	* pow(val[1],-3) + omega0_phi*( 0.5 * val[3]*val[3] + 0.5*Q*Q/val[2]/val[2]/pow(val[1],6.) + V));
	
	/* val[2] = sigma
	val[3] = dval[2]/dx */
	
	dval_dx[2] = val[3];
	
	dval_dx[3] = - 3. * dval_dx[1]/val[1] * val[3] + pow(Q,2.)/pow(val[2],3.)/pow(val[1],6.) - dV_dsigma;
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

double f(double CI, double val[], double x1, int nvar, void (*derivs)(double, double [], double []))
{
  time_t T1;
  T1 = time(0);
  double *y;
  y=dvector(1,nvar);
  y[1]=val[1];
  y[2]=CI;
  y[3]=val[3];
  int nstp,i;
  double x,hnext,hdid,h,eps,hmin,atemp=y[1],a_eq;
  double *yscal, *dydx, *yout;
  yscal=dvector(1,nvar);
  dydx=dvector(1,nvar);
  yout=dvector(1,nvar);
  eps=1.0e-10;
  h=1.e-25;
  hmin=0.0;
  x=x1;
  for (nstp=1;nstp<=MAXSTP;nstp++) { 
    nstp=1;
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (difftime(time(0),T1)>4||y[1]>0.01)
      { a_eq=omega0_rad/(omega0_b+omega0_phi*(0.5*y[3]*y[3]+ 0.5*Q*Q/y[2]/y[2]/pow(y[1],6.)+V_pot(y[2]))*pow(y[1],3.));
	free_dvector(dydx,1,nvar);
	free_dvector(yscal,1,nvar);
	free_dvector(yout,1,nvar);
	return a_eq;
      }
    rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
    h=hnext;
  }
  printf("ca n'a pas ete jusqu'a x2 \n");
  free_dvector(dydx,1,nvar);
  free_dvector(yscal,1,nvar);
  free_dvector(yout,1,nvar);
  return 1.;
}


double g(double CI, double val[], double x1, int nvar, void (*derivs)(double, double [], double []))
{
  time_t T1;
  T1 = time(0);
  double *y;
  y=dvector(1,nvar);
  y[1]=val[1];
  y[2]=CI;
  y[3]=val[3];
  int nstp,i;
  double x,hnext,hdid,h,eps,hmin,atemp=y[1],a_eq;
  double *yscal, *dydx, *yout;
  yscal=dvector(1,nvar);
  dydx=dvector(1,nvar);
  yout=dvector(1,nvar);
  eps=1.0e-10;
  h=1.e-25;
  hmin=0.0;
  x=x1;
  for (nstp=1;nstp<=MAXSTP;nstp++) { 
    nstp=1;
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (difftime(time(0),T1)>4||y[1]>0.01)
      { 
        a_eq=omega0_rad/(omega0_b+omega0_phi*(0.5*y[3]*y[3]+ 0.5*Q*Q/y[2]/y[2]/pow(y[1],6.)+V_pot(y[2]))*pow(y[1],3.));
	free_dvector(dydx,1,nvar);
	free_dvector(yscal,1,nvar);
	free_dvector(yout,1,nvar);
	return fabs(a_eq-a_eq_exp);
      }
    rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
    h=hnext;
  }
  printf("ca n'a pas ete jusqu'a x2 \n");
  free_dvector(dydx,1,nvar);
  free_dvector(yscal,1,nvar);
  free_dvector(yout,1,nvar);
  return 1.;
}


void golden(double *ax, double bx, double *cx, double y[], double xin, int nvar, void (*derivs)(double, double [], double []), double (*f)(double, double [], double, int, void (*)(double, double [], double [])), double tol, double *xmin)                                    
{
  double f1,f2,x0,x1,x2,x3;
  x0=*ax; 
  x3=*cx; 
  if (fabs(*cx-bx) > fabs(bx-*ax)) { 
    x1=bx;
    x2=bx+C*(*cx-bx); 
  } else {
    x2=bx;
    x1=bx-C*(bx-*ax);
  }
  f1=(*f)(x1,y,xin,nvar,derivs); 
  f2=(*f)(x2,y,xin,nvar,derivs);
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) { 
      SHFT3(x0,x1,x2,RM*x1+C*x3) 
      SHFT2(f1,f2,(*f)(x2,y,xin,nvar,derivs)) 
      *ax=x0;
      *cx=x3;
    } else { 
      SHFT3(x3,x2,x1,RM*x2+C*x0)
      SHFT2(f2,f1,(*f)(x1,y,xin,nvar,derivs))
      *ax=x0;
      *cx=x3;	
    }
  } 
  if (f1 < f2) {
    *xmin=x1; 
    return;
  } else {
    *xmin=x2;
    return;
  }
}

/*---------------------------------------------------------------------------------------------------------------*/
int main()

{
  double h,T0,x0,x1,val [N_VARIABLES+1],pas,erreur,junk,a_eq;
  int test;
  FILE *ext;
  
  
  	x0=1.e-25 ;                /* valeur initiale de t */
  	val[1]= 1.e-12; 	   /* valeur initiale du parametre d'expansion */  
  	val[2]= 1.e4;  	   /* valeur initiale de sigma */
  	val[3]= 1.e7; 		   /* valeur initiale de dsigma/dt */
 
  erreur =1.e-10;
  
  h=0.7;         /* Constante de Hubble. */
  T0=2.728 ;      /* Parametre cosmologique. */

  
  H0= h * 100. / 3.0856e19;
  m=1.e-24 * e/h_;      /* voir article 2 */
  Q=1.e-25;              /* tuning pour trouver omega_b/omega_phi = 0.2 */
  lambda= 1.e-98*pow(c,5)/h_;
  rho_0_c =  3.0 * H0 * H0 / (8.0 * pi * G);
  
  omega0_rad = (1+3.*7./8.*pow(4./11.,4./3.))*pi * pi / 15. * pow(k * T0, 4.0) / pow(h_,3.0)/pow(c,5.0) / rho_0_c ;
  omega0_b = 0.024/(h*h);
  omega0_phi = 1. / rho_0_c ;
  
  V0= 0.69 * rho_0_c ;   /* on veut omega_phi^0 = 0.96 */ 
  a_eq_exp=1./3392.;
  a_eq=1./3392.;

  x1=100.e9*3600*24*365;     /* age de l'Univers < 100 milliards d'annees */
  cut = 8.e-4;


  double ax=2.7e3,bx=3.e+03,cx=5.e4,xmin=0.;
  printf("%.10e \n",f(ax,val,x0,3,fonction_diff));
  printf("a_eq_exp=%.10e \n",a_eq);
  printf("%.10e \n",f(cx,val,x0,3,fonction_diff));
  if ( (f(ax,val,x0,3,fonction_diff) < a_eq) || (f(cx,val,x0,3,fonction_diff) > a_eq) ) {
  printf("Pb avec les CI pour la minimisation !! \n");
  return 0;
  }
  
  golden(&ax,bx,&cx,val,x0,3,fonction_diff,g,0.01,&xmin);
  printf("a_eq (calculée par le programme)=%.10e \n",f(xmin,val,x0,3,fonction_diff));
	printf("%.10e \n",xmin);
  	val[2]= xmin; 	   /* valeur initiale de sigma */


  return 0;
}
