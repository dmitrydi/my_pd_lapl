/*
 * series.h
 *
 *  Created on: 26 но€б. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "series_wrappers.h"
#include <cmath>
#include <vector>
#include <limits>
#include <exception>
#include <string>
#include <utility>
#include <cstdint>
#include <iostream>
#include <array>

static const int64_t MAX_LONG_VAL = std::numeric_limits<int64_t>::max();
static const int MAX_POWER_2 = static_cast<int>(std::log2(MAX_LONG_VAL));
static const int WIJN_MAX_R = 100;
static const int WIJN_STEP = 5;
static const int MAX_WIJN_ITER = static_cast<int>(std::log2(MAX_LONG_VAL/WIJN_MAX_R))/WIJN_STEP;
static const double WIJN_EPS = 1e-30;
static const double SMALL_CONST = std::numeric_limits<double>::min()*10.0;
static const double BIG_CONST = std::numeric_limits<double>::max();

struct Epsalg {
//Convergence acceleration of a sequence by the  algorithm. Initialize by calling the constructor
//with arguments nmax, an upper bound on the number of terms to be summed, and epss, the
//desired accuracy. Then make successive calls to the function next, with argument the next
//partial sum of the sequence. The current estimate of the limit of the sequence is returned by
//next. The flag cnvgd is set when convergence is detected.
	std::vector<double> e; //Workspace.
	int n,ncv;
	bool cnvgd;
	double eps,small,big,lastval,lasteps; //Numbers near machine underflow and
										  //overflow limits.
	Epsalg(int nmax, double epss) : e(nmax), n(0), ncv(0),
		cnvgd(0), eps(epss), lastval(0.) {
			small = SMALL_CONST;
			big = BIG_CONST;
		}
	double next(double sum) {
		double diff,temp1,temp2,val;
		e[n]=sum;
		temp2=0.0;
		for (int j=n; j>0; j--) {
			temp1=temp2;
			temp2=e[j-1];
			diff=e[j]-temp2;
			if (std::abs(diff) <= small)
			e[j-1]=big;
			else
			e[j-1]=temp1+1.0/diff;
		}
		n++;
		val = (n & 1) ? e[0] : e[1]; //Cases of n even or odd.
		if (std::abs(val) > 0.01*big) val = lastval;
		lasteps = std::abs(val-lastval);
		if (lasteps > eps) ncv = 0;
		else ncv++;
		if (ncv >= 3) cnvgd = 1;
		return (lastval = val);
	}
};

struct Levin {
	//Convergence acceleration of a sequence by the Levin transformation. Initialize by calling the
	//constructor with arguments nmax, an upper bound on the number of terms to be summed, and
	//epss, the desired accuracy. Then make successive calls to the function next, which returns
	//the current estimate of the limit of the sequence. The flag cnvgd is set when convergence is
	//detected.
	std::array<double, WIJN_MAX_R> numer,denom; //Numerator and denominator computed via (5.3.16).
	int n,ncv;
	bool cnvgd;
	double small,big; //Numbers near machine underflow and overflow limits.
	double eps,lastval,lasteps;
	// numer(nmax), denom(nmax),
	Levin(int nmax, double epss) : n(0), ncv(0),
			cnvgd(0), eps(epss), lastval(0.) {
		small=SMALL_CONST;
		big=BIG_CONST;
	}
	double next(double sum, double omega, double beta=1.) {
	//Arguments: sum, the nth partial sum of the sequence; omega, the nth remainder estimate
	//!n, usually from (5.3.19); and the parameter beta, which should usually be set to 1, but
	//sometimes 0.5 works better. The current estimate of the limit of the sequence is returned.
		int j;
		double fact,ratio,term,val;
		term=1.0/(beta+n);
		denom[n]=term/omega;
		numer[n]=sum*denom[n];

		if (n > 0) {
			ratio=(beta+n-1)*term;
			for (j=1;j<=n;j++) {
				fact=(n-j+beta)*term;
				numer[n-j]=numer[n-j+1]-fact*numer[n-j];
				denom[n-j]=denom[n-j+1]-fact*denom[n-j];
				term=term*ratio;
			}
		}
		n++;
		if (std::abs(denom[0]) < small) {
			val = lastval;
		} else {
			val = numer[0]/denom[0];
		}
		lasteps = std::abs(val-lastval);
		if (lasteps <= eps) ncv++;
		if (ncv >= 2) cnvgd = 1;
		lastval = val;
		return lastval;
	}
};

struct RetFormBack {
	double sum;
	double err;
	int cntr;
};

template <typename F>
RetFormBack BackSum(SeriesWrapper<F>& sw, const double eps, const int64_t step, const int64_t maxit) {
// Summates series wraped in sw
// Summation in performed by chunks of size step, each is summated in backward direction (over interval of indices [i, i+ step) )
	double sum = 0., old_sum = 0., d;
	const double small = std::numeric_limits<double>::min()*10.0;
	int i;
	for (i = 0; i < maxit; ++i) {
		sum += PartSum(sw, step);
		d = std::abs((sum - old_sum)/sum);
		if ((std::abs(sum - old_sum) < eps*std::abs(sum)) || (std::abs(sum) < small)) {
			return {sum, d, i};
		}
		old_sum = sum;
	}
	throw std::runtime_error("In BackSum: series did not converge in" + std::to_string(step*maxit));
	return {};
}

template <typename F>
double PartSum(SeriesWrapper<F>& sw, const int64_t steps) {
// Calculate partial sum from sw(cntr) to sw(cntr + step) in backward direction
	double psum = 0.;
	int64_t next_counter;
	try {
		next_counter = sw.jump(steps); // increment counter by step and remember incremented value;
	}
	catch (const std::underflow_error& ex) {
		std::cout<<ex.what();
		throw;
	}
	for (int i = 0; i < steps; ++i) {
		sw.dec_cntr();
		psum += sw.get_val();
	}
	sw.assign_cntr(next_counter); // assign counter with incremented value;

	return psum;
}

template <typename F>
double LevinDSum(F& func, const bool from_zero_i, const double eps = 1e-9) {
	Levin lev(WIJN_MAX_R, eps);
	double sum = from_zero_i ? func(0): 0.; // for series starting with 0 index calculate zero-member
											// as in WijnWaarden transformation seriaes start from index 1
	double omega;
	double ans;
	const double beta = 1.;
	int sign = 1;
	for (int64_t r = 1; r < WIJN_MAX_R; ++r) {
		auto wrap = Wrappers::Wijn<F>(func, r);
		omega = sign * BackSum(wrap, WIJN_EPS, WIJN_STEP, MAX_WIJN_ITER).sum;
		ans = lev.next(sum, omega, beta);
		if (lev.cnvgd) {
			return ans;
		}
		sum += omega;
		sign *= -1;
	}
	throw std::runtime_error("Levin did not converge in: " + std::to_string(WIJN_MAX_R));
}

template <typename F>
double EpsSum(F& func, const bool from_zero_i, const double eps = 1e-9) {
	Epsalg epsalg(WIJN_MAX_R, eps);
	double sum = from_zero_i ? func(0): 0.;
	double next_member;
	double ans;
	int sign = 1;
	for (int64_t r = 1; r < WIJN_MAX_R; ++r) {
		auto wrap = Wrappers::Wijn<F>(func, r);
		next_member = sign * BackSum(wrap, WIJN_EPS, WIJN_STEP, MAX_WIJN_ITER).sum;
		ans = epsalg.next(sum);
		if (epsalg.cnvgd) {
			return ans;
		}
		sum += next_member;
		sign *= -1;
	}
	throw std::runtime_error("Levin did not converge in: " + std::to_string(WIJN_MAX_R));
}
