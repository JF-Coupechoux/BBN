#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <locale.h>
#include "nrutil.c"
#include "nrutil.h"

#define SORTIE "MLL.out"
#define pi 3.14159265359
#define G 6.67259e-11
#define c 2.997924562e8
#define h_ 1.05457266913e-34
#define e 1.60217733e-19
#define k 1.380658e-23     /* Quelques constantes. */

int main()
{
double h,T0,rho_0_c,V0,omega0_rad,omega0_b,omega_phi,y0,y,omega_rad,omega_b,H0;
FILE *output;
  h=0.673;         /* Constante de Hubble. */
  T0=2.728 ;      /* Parametre cosmologique. */
  H0= h * 100. / 3.0856e19;
  rho_0_c =  3.0 * H0 * H0 / (8.0 * pi * G);
  V0=0.69;
  omega0_rad = (1+3.*7./8.*pow(4./11.,4./3.))*pi * pi / 15. * pow(k * T0, 4.0) / pow(h_,3.0)/pow(c,5.0) / rho_0_c ;
  omega0_b = 0.024/(h*h);
  y0=2.5011935007e-04;
  y=2.5011935007e-04;
while(y<1.1) {
		omega_rad = omega0_rad * pow(y,-4.);
		omega_b = omega0_b * pow(y,-3.);
		omega_phi= 2.1649960949e+10*pow(y/y0,-3);
		output=fopen(SORTIE,"a");
		fprintf(output,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e \n", 0., y, 0., 0., omega_rad, omega_b, omega_phi);
		fclose(output);
		y=y*pow(10.,0.05);	
	}
}
