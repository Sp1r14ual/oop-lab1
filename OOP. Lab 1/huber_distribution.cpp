#include "huber_distribution.h"

double Huber(double x, double v, double K, double scale, double shift)
{
    return (1. / (sqrt(2. * M_PI) * K) * (abs((x - shift) / scale) <= v ? exp(-pow((x - shift) / scale, 2.) / 2.) : exp(pow(v, 2.) / 2. - v * abs((x - shift) / scale)))) / scale;
}

double phi(double x)
{
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

double phi_lower(double x)
{
    return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}

double huber_expected_value(double shift)
{
    return shift;
}

double huber_variance(double v, double K)
{
    return 1. + 2. * phi_lower(v) * (pow(v, 2.) + 2.) / (pow(v, 3.) * K);
}

double huber_asymmetry()
{
    return 0.;
}


double huber_kurtosis(double v, double K)
{
    return (3. * (2. * phi(v) - 1.) + 2. * phi_lower(v) * (24. / pow(v, 5.) + 24. / pow(v, 3.) + 12. / v + v)) / (pow(huber_variance(v, K), 2.) * K) - 3.;
}

double P(double v, double K)
{
    return (2. * phi(v) - 1.) / K;
}

double K(double v)
{
    return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}

double calculate_x(double v, double K, double scale, double shift)
{
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<> d(0, 1);

    //шаг 1
    double r1 = d(gen);


    if (r1 <= P(v, K))
    {
        //шаг 2
        double r2, r3, x1;

        do {
            r2 = d(gen);
            r3 = d(gen);
            x1 = sqrt(-2 * log(r2)) * cos(2 * M_PI * r3);
            //double x1 = sqrt(-2 * log(r2)) * sin(2 * M_PI * r3)
        } while (!(-v <= x1 && x1 <= v)); //шаг 3

        return x1 * scale + shift;
    }
    else
    {
        //шаг 4
        double r4 = d(gen);
        double x2 = v - log(r4) / v;

        //шаг 5
        return r1 < (1 + P(v, K)) / 2 ? x2 * scale + shift : -x2 * scale + shift;
    }
}

vector<double> generate_sequence(int n, double v, double K)
{
    vector<double> x_s;

    for (int i = 0; i < n; i++)
    {
        double x = calculate_x(v, K);
        x_s.push_back(x);
    }

    return x_s;
}