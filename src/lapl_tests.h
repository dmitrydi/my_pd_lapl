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
#include "auxillary.h"
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <random>
#include <ctime>

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
	const double TEST_SEXP_EPS = 1e-14;
	const double NTESTS = 1000'000;
	const int N = 20;
	const int KMAX = 1000;
	const double xed_min = 1;
	const double xed_max  = 1000;
	const double yed_min = 1;
	const double yed_max = 1000;
	const double u_min = 1e-9;
	const double u_max = 1000;

	const auto xeds = LogSpaced(xed_min, xed_max, N);
	const auto yeds = LogSpaced(yed_min, yed_max, N);
	const auto us = LogSpaced(u_min, u_max, N);

	mt19937 rng;
	rng.seed(time(nullptr));
	uniform_int_distribution<std::mt19937::result_type> dist(0, N-1);
	uniform_int_distribution<std::mt19937::result_type> k_dist(1, KMAX);

	for (int i = 0; i < NTESTS; ++i) {
		auto x_ind = dist(rng);
		auto y_ind = dist(rng);
		auto u_ind = dist(rng);
		double u = us[u_ind];
		double xed = xeds[x_ind];
		double yed = yeds[y_ind];
		auto ek = LaplFunc::ek(u, xed, 0);
		auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS/100., LaplFunc::SumAlgo::epsalg);
		int64_t k = k_dist(rng);
		auto s = sexp(ek(k));
		ASSERT_CLOSE(s, SexpGuaranteed(yed, ek(k)), TEST_SEXP_EPS);
	}
}


void Test_sexp_speed() {
	const double TEST_SEXP_EPS = 1e-14;
		const int NTESTS = 1000'000;
		const int N = 20;
		const int KMAX = 10;
		const double xed_min = 1;
		const double xed_max  = 1000;
		const double yed_min = 1;
		const double yed_max = 1000;
		const double u_min = 1e-9;
		const double u_max = 1000;

		const auto xeds = LogSpaced(xed_min, xed_max, N);
		const auto yeds = LogSpaced(yed_min, yed_max, N);
		const auto us = LogSpaced(u_min, u_max, N);

		mt19937 rng;
		rng.seed(time(nullptr));
		uniform_int_distribution<std::mt19937::result_type> dist(0, N-1);
		uniform_int_distribution<std::mt19937::result_type> k_dist(1, KMAX);

		double s;

		ostringstream os;

		os << "Num tests: ";
		os << fixed << to_string(NTESTS);
		os << ", duration: ";

		{LOG_DURATION(os.str());
			for (int i = 0; i < NTESTS; ++i) {
				auto x_ind = dist(rng);
				auto y_ind = dist(rng);
				auto u_ind = dist(rng);
				double u = us[u_ind];
				double xed = xeds[x_ind];
				double yed = yeds[y_ind];
				auto ek = LaplFunc::ek(u, xed, 0);
				auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS/100., LaplFunc::SumAlgo::epsalg);
				int64_t k = k_dist(rng);
				s = sexp(ek(k));
			}
		}
		cerr << s << endl;
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

