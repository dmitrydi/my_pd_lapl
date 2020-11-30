/*
 * quadrature.h
 *
 *  Created on: 26 íîÿá. 2020 ã.
 *      Author: Dmitry_Di
 */

#pragma once

#include "interp_1d.h"
#include <cmath>

struct Quadrature {
// Abstract base class for elementary quadrature algorithms
	int n;
	virtual double next() = 0;
};

template <class T>
struct Trapzd: Quadrature {
	// Routine implementing the extended trapezoidal rule.
	T& func;
	double a, b, s;

	Trapzd(T& funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb), s(0) {n=0;}

	double next() {
		double x,tnm,sum,del;
		int it,j;
		n++;
		if (n == 1) {
			s=0.5*(b-a)*(func(a)+func(b));
			return s;
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
};

template <class T>
double qtrap(T& func, const double a, const double b, const double eps=1.0e-10) {
	//Returns the integral of the function or functor func from a to b. The constants EPS can be
	//set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
	//allowed number of steps. Integration is performed by the trapezoidal rule.
	const int JMAX=100;
	double s,olds=0.0; 						//Initial value of olds is arbitrary.
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		s=t.next();
		if (j > 5) 	{						//Avoid spurious early convergence.
			if (abs(s-olds) < eps*abs(olds) ||
					(s == 0.0 && olds == 0.0))
				return s;
		}
		olds=s;
	}
	throw("Too many steps in routine qtrap");
}

template<class T>
double qsimp(T &func, const double a, const double b, const double eps=1.0e-10) {
//Returns the integral of the function or functor func from a to b. The constants EPS can be
//set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
//allowed number of steps. Integration is performed by Simpson’s rule.
	const int JMAX=20;
	double s,st,ost=0.0,os=0.0;
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		st=t.next();
		s=(4.0*st-ost)/3.0;
		if (j > 5) { //Avoid spurious early convergence.
			if (abs(s-os) < eps*abs(os) ||
					(s == 0.0 && os == 0.0))
				return s;
		}
		os=s;
		ost=st;
	}
	throw("Too many steps in routine qsimp");
}

template <class T>
double qromb(T &func, double a, double b, const double eps=1.0e-10) {
//Returns the integral of the function or functor func from a to b. Integration is performed by
//Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.
	const int JMAX=20, JMAXP=JMAX+1, K=5;
	//Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
	//JMAX limits the total number of steps; K is the number of points used in the
	//extrapolation.
	vector<double> s(JMAX),h(JMAXP); //These store the successive trapezoidal approxi
	Poly_interp polint(h,s,K);     //mations and their relative stepsizes.
	h[0]=1.0;
	Trapzd<T> t(func,a,b);
	for (int j=1;j<=JMAX;j++) {
		s[j-1]=t.next();
		if (j >= K) {
			double ss=polint.rawinterp(j-K,0.0);
			if (abs(polint.dy) <= eps*abs(ss)) return ss;
		}
		h[j]=0.25*h[j-1];
		//This is a key step: The factor is 0.25 even though the stepsize is decreased by only
		//0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1),
		//not just a polynomial in h.
	}
	throw("Too many steps in routine qromb");
}
