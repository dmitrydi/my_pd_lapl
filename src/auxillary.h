/*
 * auxillary.h
 *
 *  Created on: 7 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "test_runner.h"
#include <cmath>
#include <vector>
#include <random>

std::vector<double> LogSpaced(const double min, const double max, const int steps) {
	std::vector<double> ans(steps);
	double ln_min = std::log(min);
	double ln_max = std::log(max);
	double d = (ln_max - ln_min)/(steps - 1);
	for (int i = 0; i < steps; ++i) {
		ans[i] = std::exp(ln_min + i*d);
	}
	return ans;
}

void TestLogSpaced() {
	double min = 0.01;
	double max = 1000;
	int steps = 100;
	double eps = 1.e-12;
	auto ans = LogSpaced(min, max, steps);
	ASSERT_CLOSE(ans.front(), min, eps);
	ASSERT_CLOSE(ans.back(), max, eps);
	ASSERT_CLOSE(ans[10], 0.031992671378, eps);
	ASSERT_EQUAL(ans.size(), steps);
}
