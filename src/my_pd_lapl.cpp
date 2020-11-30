//============================================================================
// Name        : my_pd_lapl.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "bessel.h"
#include "profile.h"
#include "test_runner.h"
#include "quadrature.h"
#include "qgaus.h"
#include "series.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <cmath>
#include <limits>
#include <cstdint>

using namespace std;

struct F {
	double operator()(int64_t n) {
		++calls;
		return 1./(static_cast<double>(n)*static_cast<double>(n));
	}
	int calls = 0;
};

ostream& operator<<(ostream& os, pair<double, double> p) {
	os << scientific <<  p.first << ' ' << p.second;
	return os;
}


int main() {
	F f;
	double ans;
	double eps = 1e-9;
	double corr_ans = PI*PI/6.;
	{LOG_DURATION("Levin D-sum")
		for (int i = 0; i < 500000; i++) {
		ans = LevinDSum(f, false, eps);}
	}

	cout << ans << ' ' << corr_ans << ' ' << abs(ans-corr_ans)/corr_ans << endl;
	cout << f.calls << endl;
	return 0;
}
