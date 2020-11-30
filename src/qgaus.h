/*
 * qgaus.h
 *
 *  Created on: 26 íîÿá. 2020 ã.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>

void gauleg(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w);

struct GaussIntegrator {
	int n;
	std::vector<double> x, w;
	explicit GaussIntegrator(const int n_): n(n_), x(n), w(n) {
		gauleg(-1., 1., x, w);
	}
	template <class T>
	double Integrate(T& func, const double a, const double b) {
		double xm=0.5*(b+a);
		double xr=0.5*(b-a);
		double s=0;
		int jmax = x.size()/2;
		for (int j=0;j<jmax;j++) {
			double dx=xr*x[j];
			s += w[j]*(func(xm+dx)+func(xm-dx));
		}
		return s *= xr;
	}
};

template <class T>
double qgaus(T& func, const double a, const double b)
//Returns the integral of the function or functor func between a and b, by ten-point Gauss-
//Legendre integration: the function is evaluated exactly ten times at interior points in the range
//of integration.
{
//Here are the abscissas and weights:
	static const double x[]={0.1488743389816312,0.4333953941292472,
							 0.6794095682990244,0.8650633666889845,0.9739065285171717};
	static const double w[]={0.2955242247147529,0.2692667193099963,
								0.2190863625159821,0.1494513491505806,0.0666713443086881};
	double xm=0.5*(b+a);
	double xr=0.5*(b-a);
	double s=0; //Will be twice the average value of the function, since the
				//ten weights (five numbers above each used twice)
				//sum to 2.
	for (int j=0;j<5;j++) {
		double dx=xr*x[j];
		s += w[j]*(func(xm+dx)+func(xm-dx));
	}
	return s *= xr; //Scale the answer to the range of integration.
}

void gauleg(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w)
//Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
//and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point
//quadrature formula.
{
	const double EPS=1.0e-14; //EPS is the relative precision.
	double z1,z,xm,xl,pp,p3,p2,p1;
	int n=x.size();
	int m=(n+1)/2; //The roots are symmetric in the interval, so
	xm=0.5*(x2+x1); //we only have to find half of them.
	xl=0.5*(x2-x1);
	for (int i=0;i<m;i++) { //Loop over the desired roots.
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		//Starting with this approximation to the ith root, we enter the main loop of refinement
		//by Newton’s method.
		do {
			p1=1.0;
			p2=0.0;
			for (int j=0;j<n;j++) { //Loop up the recurrence relation to get the
				p3=p2;				//Legendre polynomial evaluated at z.
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			//p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp; //Newton’s method.
		} while (abs(z-z1) > EPS);
		x[i]=xm-xl*z; 	  //Scale the root to the desired interval,
		x[n-1-i]=xm+xl*z; //and put in its symmetric counterpart.
		w[i]=2.0*xl/((1.0-z*z)*pp*pp); //Compute the weight
		w[n-1-i]=w[i]; //and its symmetric counterpart.
	}
}
