/*
 * series_wrappers.h
 *
 *  Created on: 30 но€б. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <cstdint>
#include <exception>
#include <string>
#include <stdexcept>

template <typename Func>
class SeriesWrapper {
public:
	virtual void inc_cntr() = 0;
	virtual void dec_cntr() = 0;
	virtual int64_t jump(int64_t steps) = 0;
	virtual double get_val() = 0;
	virtual void assign_cntr(const int64_t val) = 0;
};

namespace Wrappers {
template <typename F>
class Forward: public SeriesWrapper<F> {
public:
	Forward(F& fun, int64_t i_min): f(fun), i(i_min) {};
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
class Backward: public SeriesWrapper<F> {
public:
	Backward(F& fun, int64_t i_min): f(fun), i(i_min) {};
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
class Wijn: public SeriesWrapper<F> {
public:
	Wijn(F& fun, int64_t r_): f(fun), r(r_) {mult = 1;};
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
}
