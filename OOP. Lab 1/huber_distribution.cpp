#include "huber_distribution.h"

HuberDistribution* init_huber_distribution(double v, double K, double scale, double shift)
{
    HuberDistribution* HB = new HuberDistribution();
    HB->v = v;
    HB->K = K;
    HB->scale = scale;
    HB->shift = shift;
    return HB;
}

double Huber(double x, HuberDistribution* HB)
{
    return (1. / (sqrt(2. * M_PI) * HB->K) * (abs((x - HB->shift) / HB->scale) <= HB->v ? exp(-pow((x - HB->shift) / HB->scale, 2.) / 2.) : exp(pow(HB->v, 2.) / 2. - HB->v * abs((x - HB->shift) / HB->scale)))) / HB->scale;
}

double phi(double x)
{
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

double phi_lower(double x)
{
    return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}

double huber_expected_value(HuberDistribution* HB)
{
    return HB->shift;
}

double huber_variance(HuberDistribution* HB)
{
    return 1. + 2. * phi_lower(HB->v) * (pow(HB->v, 2.) + 2.) / (pow(HB->v, 3.) * HB->K);
}

double huber_asymmetry(HuberDistribution* HB)
{
    return 0.;
}


double huber_kurtosis(HuberDistribution* HB)
{
    return (3. * (2. * phi(HB->v) - 1.) + 2. * phi_lower(HB->v) * (24. / pow(HB->v, 5.) + 24. / pow(HB->v, 3.) + 12. / HB->v + HB->v)) / (pow(huber_variance(HB), 2.) * HB->K) - 3.;
}

double P(HuberDistribution* HB)
{
    return (2. * phi(HB->v) - 1.) / HB->K;
}

double K(double v)
{
    return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}

double calculate_x(HuberDistribution* HB)
{
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<> d(0, 1);

    //шаг 1
    double r1 = d(gen);


    if (r1 <= P(HB))
    {
        //шаг 2
        double r2, r3, x1;

        do {
            r2 = d(gen);
            r3 = d(gen);
            x1 = sqrt(-2 * log(r2)) * cos(2 * M_PI * r3);
            //double x1 = sqrt(-2 * log(r2)) * sin(2 * M_PI * r3)
        } while (!(-HB->v <= x1 && x1 <= HB->v)); //шаг 3

        return x1 * HB->scale + HB->shift;
    }
    else
    {
        //шаг 4
        double r4 = d(gen);
        double x2 = HB->v - log(r4) / HB->v;

        //шаг 5
        return r1 < (1 + P(HB)) / 2 ? x2 * HB->scale + HB->shift : -x2 * HB->scale + HB->shift;
    }
}

