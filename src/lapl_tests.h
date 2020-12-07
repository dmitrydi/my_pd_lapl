/*
 * lapl_tests.h
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */
#pragma once

#include "lapl_functions.h"
#include "test_runner.h"
#include "series.h"
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>

using namespace std;

template <typename T>
bool ParseLine(ifstream& input, const int pos, T& val) {
	string line;
	T ans;
	if (getline(input, line)) {
		istringstream is(line);
		for (int i = 0; i < pos; ++i) {
			is >> ans;
		}
		val = ans;
		return true;
	} else {
		return false;
	}
}

double SexpGuaranteed(const double yed, const double ek) {
// sums exp(-2*m*yed*ek) to almost numeric limit precision
	const double tiny = numeric_limits<double>::min();
	const double ln_tiny = log(tiny);
	const double dum = -0.5*ln_tiny/yed/ek;
	const int64_t max_member = static_cast<int>(dum);
	double sum = 0;
	for (int64_t i = max_member; i >= 1; --i) {
		sum += exp(-2*i*yed*ek);
	}
	return sum;
}

void Test_sexp() {
	const double TEST_SEXP_EPS = 1e-15;
	const double xed_max = 1000;
	const double u_min = 1e-9;
	const double yed_min = 1;
	auto ek = LaplFunc::ek(u_min, xed_max, 0);
	auto sexp = LaplFunc::sexp(yed_min, TEST_SEXP_EPS);
	ASSERT_CLOSE(sexp(ek(1)), SexpGuaranteed(yed_min, ek(1)), TEST_SEXP_EPS);
}

void Test_sexp_speed() {
	const double TEST_SEXP_EPS = 1e-15;
	const double xed_max = 1000;
	const double u_min = 1e-9;
	const double yed_min = 1;
	auto ek = LaplFunc::ek(u_min, xed_max, 0);


	double d;
	{
		auto sexp = LaplFunc::sexp(yed_min, TEST_SEXP_EPS);
		int n = 50'000;

		LOG_DURATION("Sexp with Epsalg " + to_string(n) + " iterations");
		for (int i = 0; i < n; ++i) {
			d = sexp(ek(1));
		}
	}
	{
		auto sexp = LaplFunc::sexp(yed_min, TEST_SEXP_EPS, LaplFunc::SumAlgo::levin);
		int n = 50'000;

		LOG_DURATION("Sexp with Levin " + to_string(n) + " iterations");
		for (int i = 0; i < n; ++i) {
			d = sexp(ek(1));
		}
	}
	{
		int n = 500;

		LOG_DURATION("Plain summation to numeric limit " + to_string(n) + " iterations");
		for (int i = 0; i < n; ++i) {
			d = SexpGuaranteed(yed_min, ek(1));
		}
	}
	cerr << d << endl;
}

void Test_sexp_get_k() {
	const double xed_max = 1000;
	const double u_min = 1e-9; // td = 1e+9;
	const double yed_min = 1;
	auto ek = LaplFunc::ek(u_min, xed_max, 0);
	auto sexp = LaplFunc::sexp(yed_min, 1e-15);
	auto func = [&sexp, ek](int64_t r) { return sexp.get_k(ek(1), r);};
	for (int64_t r = 1; r < WIJN_MAX_R; ++r) {
		auto wrap = Wrappers::Wijn(func, r);
		cerr << BackSum(wrap, 1e-30, WIJN_STEP, MAX_WIJN_ITER).sum << endl;
	}
}


void Test_ek() {
	const double TEST_EK_EPS = 1e-14;
	const string fname  = "src\\test_data\\ek_test.txt";
	ifstream input(fname);
	if (!input.is_open()) {
		cerr << "Could not open " + fname << endl;
	} else {
		double u, ksi, alpha;
		ParseLine(input, 1, u);
		ParseLine(input, 1, ksi);
		ParseLine(input, 1, alpha);
		vector<double> expected;
		for(; ;) {
			double v;
			if (ParseLine(input, 2, v)) {
				expected.push_back(v);
			} else {
				break;
			}
		}
		LaplFunc::ek ek(u, ksi, alpha);
		for (size_t i = 1; i <= expected.size(); ++i) {
			ASSERT(abs((ek(i) - expected[i-1])/expected[i-1]) < TEST_EK_EPS);
		}
	}
}

