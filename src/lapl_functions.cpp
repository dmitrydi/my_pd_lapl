/*
 * lapl_functions.cpp
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "lapl_functions.h"

using namespace std;

static const double PI = 3.141592653589793;
static const double EPS_SEXP = 1e-30;

namespace LaplFunc {

ek::ek(const double u_, const double ksi_, const double alpha_): u(u_), ksi(ksi_), alpha(alpha_){};

double ek::operator()(const int64_t k) const {
	double dum = static_cast<double>(k)*PI/ksi;
	return sqrt(u + dum*dum + alpha*alpha);
};

sexp::sexp(const double yed_, const double eps_, const SumAlgo algo_): yed(yed_), eps(eps_), cntr(0), algo(algo_) {};

double sexp::operator()(const double e) const {
	auto fun = [e, this](int64_t m) { return this->get_k(e, m); };
	if (algo == SumAlgo::levin) return LevinDSum(fun, false, eps);
	if (algo == SumAlgo::epsalg) return EpsSum(fun, false, eps);
	throw;
}

double sexp::get_k(const double e, const int64_t k) const {
	return exp(-2*k*yed*e);
}

int64_t sexp::count() const {
	return cntr;
}

}



