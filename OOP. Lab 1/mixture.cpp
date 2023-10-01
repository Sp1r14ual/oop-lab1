#include "mixture.h"

double mixture(double x, double v1, double v2, double K1, double K2, double p, double scale1, double shift1, double scale2, double shift2)
{
    return (1 - p) * Huber(x, v1, K1, scale1, shift1) + p * Huber(x, v2, K2, scale2, shift2);
}

double mixture_expected_value(double p, double shift1, double shift2)
{
	return (1 - p) * huber_expected_value(shift1) + p * huber_expected_value(shift2);
}

double mixture_variance(double v, double K, double p, double shift1, double shift2)
{
	return 0;
}

double mixture_asymmetry()
{
	return 0;
}

double mixture_kurtosis()
{
	return 0;
}