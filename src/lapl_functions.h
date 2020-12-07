/*
 * lapl_functions.h
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "series.h"
#include "series_wrappers.h"
#include <cmath>
#include <cstdint>


namespace LaplFunc {

enum SumAlgo {
	levin,
	epsalg
};

class ffunc {
public:
	virtual double operator()(const double k) const = 0;

};

class ifunc {
public:
	virtual double operator()(const int64_t k) const = 0;
};

class ek: public ifunc {
public:
	ek(const double u_, const double ksi_, const double alpha_);
	double operator()(const int64_t k) const override;
	int64_t count() const;
private:
	const double u, ksi, alpha;
};

class sexp: public ffunc {
public:
	sexp(const double yed_, const double eps_, const SumAlgo algo_ = SumAlgo::epsalg);
	double operator()(const double e) const override;
	double get_k(const double e, const int64_t k) const;
	int64_t count() const;
private:
	const double yed;
	const double eps;
	int64_t cntr;
	const SumAlgo algo;
};

}
