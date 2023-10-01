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

struct Mixture
{
	double p;
	HuberDistribution HB1;
	HuberDistribution HB2;
};

double Huber(double x, double v, double K, double scale = 1., double shift = 0.);

double phi(double x);

double phi_lower(double x);

double huber_expected_value(double shift);

double huber_variance(double v, double K);

double huber_asymmetry();

double huber_kurtosis(double v, double K);

double P(double v, double K);

double K(double v);

double calculate_x(double v, double K, double scale = 1., double shift = 0.);

vector<double> generate_sequence(int n, double v, double K);