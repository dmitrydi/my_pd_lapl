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
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <cmath>

using namespace std;

double fRand(double fMin, double fMax) {
	double f = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	return fMin + f * (fMax - fMin);
};

int main() {

	srand(time(NULL));
	const int n = 20'000'000;
	vector<double> v;
	v.reserve(n+1);
	for (int i = 0; i < n; i++) {
		double x = fRand(0.01, 10);
		if (x <= 0.0) {
			cout << x << endl;
		}
		v.push_back(x);

	}

	Bessik bess;
	double y;
	{ LOG_DURATION((to_string(n) + " numbers in custom lib"))
		for (auto x: v) {
			y = bess.k0(x);
		}
	}
	{ LOG_DURATION((to_string(n) + " numbers in standard lib"))
			for (auto x: v) {
				y = cyl_bessel_k(0, x);
			}
	}
	cout << y << endl;
	return 0;
}
