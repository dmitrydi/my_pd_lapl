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
#include "lapl_functions.h"
#include "lapl_tests.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <cmath>
#include <limits>
#include <cstdint>
#include <random>

using namespace std;

int main()
{
	TestRunner tr;
	RUN_TEST(tr, Test_ek);
	RUN_TEST(tr, Test_sexp);
	RUN_TEST(tr, Test_sexp_speed);
}

