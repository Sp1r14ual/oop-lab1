#pragma once
#include "lib.h"
#include "huber_distribution.h"

using namespace std;

double mixture(double x, double v1, double v2, double K1, double K2, double p, double scale1 = 1., double shift1 = 0., double scale2 = 1., double shift2 = 0.);

double mixture_expected_value(double p, double shift1, double shift2);
