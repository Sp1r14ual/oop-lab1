#include "mixture_distribution.h"

Mixture* init_mixture(double p, double v1, double v2, double scale1, double scale2, double shift1, double shift2)
{
	Mixture* M = new Mixture();
	HuberDistribution* HB1 = new HuberDistribution();
	HuberDistribution* HB2 = new HuberDistribution();
	double K1 = K(v1);
	double K2 = K(v2);

	HB1->v = v1;
	HB1->K = K1;
	HB1->scale = scale1;
	HB1->shift = shift1;

	HB2->v = v2;
	HB2->K = K2;
	HB2->scale = scale2;
	HB2->shift = shift2;

	M->p = p;
	M->HB1 = HB1;
	M->HB2 = HB2;

	return M;
}

double mixture(double x, Mixture* M)
{
    return (1 - M->p) * Huber(x, M->HB1) + M->p * Huber(x, M->HB2);
}

double mixture_expected_value(Mixture* M)
{
	return (1 - M->p) * huber_expected_value(M->HB1) + M->p * huber_expected_value(M->HB2);
}

double mixture_variance(Mixture* M)
{
	return (1 - M->p) * (pow(huber_expected_value(M->HB1), 2) + huber_variance(M->HB1)) + M->p * (pow(huber_expected_value(M->HB2), 2) + 
		huber_variance(M->HB2)) - pow(mixture_expected_value(M), 2);
}

double mixture_asymmetry(Mixture* M)
{
	return ((1 - M->p) * (pow((huber_expected_value(M->HB1) - mixture_expected_value(M)), 3) + 3 * (huber_expected_value(M->HB1) -
		mixture_expected_value(M)) * huber_variance(M->HB1) + pow(huber_variance(M->HB1), 3 / 2) * huber_asymmetry(M->HB1)) +
		M->p * (pow((huber_expected_value(M->HB2) - mixture_expected_value(M)), 3) + 3 * (huber_expected_value(M->HB2) - mixture_expected_value(M)) * 
			huber_variance(M->HB2) + pow(huber_variance(M->HB2), 3 / 2) * huber_asymmetry(M->HB2))) /pow(mixture_variance(M), 3 / 2);
}

double mixture_kurtosis(Mixture* M)
{
	return ((1 - M->p) * (pow((huber_expected_value(M->HB1) - mixture_expected_value(M)), 4) + 6 * huber_variance(M->HB1) * 
		pow((huber_expected_value(M->HB1) - mixture_expected_value(M)), 2) +
		4 * (huber_expected_value(M->HB1) - mixture_expected_value(M)) * pow(huber_variance(M->HB1), 3 / 2) * 
		huber_asymmetry(M->HB1) + pow(huber_variance(M->HB1), 2) * huber_kurtosis(M->HB1)) +
		M->p * (pow((huber_expected_value(M->HB2) - mixture_expected_value(M)), 4) + 6 * huber_variance(M->HB2) * pow((huber_expected_value(M->HB2) - 
			mixture_expected_value(M)), 2) + 4 * (huber_expected_value(M->HB2) - mixture_expected_value(M)) * pow(huber_variance(M->HB2), 3 / 2) * 
			huber_asymmetry(M->HB2) + pow(huber_variance(M->HB2), 2) * huber_kurtosis(M->HB2)) - 3) / pow(mixture_variance(M), 2);
}