/*
 * series.h
 *
 *  Created on: 26 но€б. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <cmath>
#include <vector>
#include <limits>
#include <exception>
#include <string>
#include <utility>
#include <cstdint>

static int64_t MAX_LONG_VAL = std::numeric_limits<int64_t>::max();
static int MAX_POWER_2 = static_cast<int>(std::log2(MAX_LONG_VAL));
static int WIJN_MAX_R = 100;
static int WIJN_STEP = 5;
static int MAX_WIJN_ITER = static_cast<int>(std::log2(MAX_LONG_VAL/WIJN_MAX_R))/WIJN_STEP;
static double WIJN_EPS = 1e-12;

void PrintConstants() {
	std::cout << "MAX_LONG_VAL " << MAX_LONG_VAL << std::endl;
	std::cout << "MAX_POWER_2 " << MAX_POWER_2 << std::endl;
	std::cout << "WIJN_MAX_R " << WIJN_MAX_R << std::endl;
	std::cout << "WIJN_STEP " << WIJN_STEP << std::endl;
	std::cout << "MAX_WIJN_ITER " << MAX_WIJN_ITER << std::endl;
}



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
			small = std::numeric_limits<double>::min()*10.0;
			big = std::numeric_limits<double>::max();
		}
	double next(double sum) {
		double diff,temp1,temp2,val;
		e[n]=sum;
		temp2=0.0;
		for (int j=n; j>0; j--) {
			temp1=temp2;
			temp2=e[j-1];
			diff=e[j]-temp2;
			if (abs(diff) <= small)
			e[j-1]=big;
			else
			e[j-1]=temp1+1.0/diff;
		}
		n++;
		val = (n & 1) ? e[0] : e[1]; //Cases of n even or odd.
		if (abs(val) > 0.01*big) val = lastval;
		lasteps = abs(val-lastval);
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
	std::vector<double> numer,denom; //Numerator and denominator computed via (5.3.16).
	int n,ncv;
	bool cnvgd;
	double small,big; //Numbers near machine underflow and overflow limits.
	double eps,lastval,lasteps;
	Levin(int nmax, double epss) : numer(nmax), denom(nmax), n(0), ncv(0),
			cnvgd(0), eps(epss), lastval(0.) {
		small=std::numeric_limits<double>::min()*10.0;
		big=std::numeric_limits<double>::max();
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
		val = abs(denom[0]) < small ? lastval : numer[0]/denom[0];
		lasteps = abs(val-lastval);
		if (lasteps <= eps) ncv++;
		if (ncv >= 2) cnvgd = 1;
		return (lastval = val);
	}
};

template <typename Func>
class SeriesWrapper {
public:
	virtual void inc_cntr() = 0;
	virtual void dec_cntr() = 0;
	virtual int64_t jump(int64_t steps) = 0;
	virtual double get_val() = 0;
	virtual void assign_cntr(const int64_t val) = 0;
};

template <typename F>
class ForwardWrapper: public SeriesWrapper<F> {
public:
	ForwardWrapper(F& fun, int i_min): f(fun), i(i_min) {};
	void inc_cntr() { ++i; };
	void dec_cntr() { --i; };
	int64_t jump(int64_t steps) { i += steps; return i;};
	void assign_cntr(const int64_t val) { i = val; };
	double get_val() {return f(i); };
private:
	F& f;
	int64_t i;
};

template <typename F>
class BackwardWrapper: public SeriesWrapper<F> {
public:
	BackwardWrapper(F& fun, int i_min): f(fun), i(i_min) {};
	void inc_cntr() { --i; };
	void dec_cntr() { ++i; };
	int64_t jump(int64_t steps) { i -= steps; return i;};
	void assign_cntr(const int64_t val) { i = val; };
	double get_val() {return f(i); };
private:
	F& f;
	int64_t i;
};


template <typename F>
class WijnWrapper: public SeriesWrapper<F> {
public:
	WijnWrapper(F& fun, int r_): f(fun), r(r_) {mult = 1;};
	void inc_cntr() {
		mult *= 2;
	};
	void dec_cntr() {
		if (mult > 1) {
			mult /= 2;
		} else {
			throw std::underflow_error("error in decrement in WijnWrapper\n");
		}
	};
	int64_t jump(int64_t steps) {
		int64_t before = mult;
		for (int i = 0; i < steps; i++) {
			mult *= 2;
		}
		if (mult < 1) {
			throw std::underflow_error(std::to_string(before));
		}

		return mult;
	}
	void assign_cntr(const int64_t val) {
		mult = val;
	}
	double get_val() {
		return static_cast<double>(mult) * f(r*mult);
	};
private:
	F& f;
	int64_t r;
	int64_t mult;
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
		d = abs((sum - old_sum)/sum);
		if ((abs(sum - old_sum) < eps*abs(sum)) || (abs(sum) < small)) {
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
	Levin lev(MAX_WIJN_ITER, eps);
	double sum = from_zero_i ? func(0): 0.; // for series starting with 0 index calculate zero-member
											// as in WijnWaarden transformation seriaes start from index 1
	double omega;
	double ans;
	const double beta = 1.;
	int sign = 1;
	for (int64_t r = 1; r < WIJN_MAX_R; ++r) {
		auto wrap = WijnWrapper<F>(func, r);
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
