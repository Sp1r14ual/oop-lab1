#pragma once
#include "lib.h"
#include "huber_distribution.h"

using namespace std;

struct Mixture
{
	double p;
	HuberDistribution* HB1;
	HuberDistribution* HB2;
};

Mixture* init_mixture(double p, double v1, double v2, double scale1 = 1, double scale2 = 1, double shift1 = 0, double shift2 = 0);

double mixture(double x, Mixture* M);

double mixture_expected_value(Mixture* M);

double mixture_variance(Mixture* M);

double mixture_asymmetry(Mixture* M);

double mixture_kurtosis(Mixture* M);
