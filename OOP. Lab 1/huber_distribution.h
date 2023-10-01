#pragma once
#include "lib.h"

using namespace std;

struct HuberDistribution
{
	double v;
	double K;
	double scale;
	double shift;
};

HuberDistribution* init_huber_distribution(double v, double K, double scale = 1, double shift = 0);

double Huber(double x, HuberDistribution* HB);

double phi(double x);

double phi_lower(double x);

double huber_expected_value(HuberDistribution* HB);

double huber_variance(HuberDistribution* HB);

double huber_asymmetry(HuberDistribution* HB);

double huber_kurtosis(HuberDistribution* HB);

double P(HuberDistribution* HB);

double K(double v);

double calculate_x(HuberDistribution* HB);