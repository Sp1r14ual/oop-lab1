#include "empirical_distribution.h"

double empirical_expected_value(int n, vector<double> x_s)
{
    double sum = 0;
    for (double& x : x_s)
        sum += x;
    return sum / n;
}

double empirical_variance(int n, vector<double> x_s)
{
    double expected_value = empirical_expected_value(n, x_s);
    double sum = 0;
    for (double& x : x_s)
        sum += pow(x - expected_value, 2);
    return sum / n;
}

double empirical_asymmetry(int n, vector<double> x_s)
{
    double expected_value = empirical_expected_value(n, x_s);
    double variance = empirical_variance(n, x_s);

    double sum = 0;
    for (double& x : x_s)
        sum += pow(x - expected_value, 3);

    return 1 / (n * pow(variance, 3 / 2)) * sum;
}

double empirical_kurtosis(int n, vector<double> x_s)
{
    double expected_value = empirical_expected_value(n, x_s);
    double variance = empirical_variance(n, x_s);

    double sum = 0;
    for (double& x : x_s)
        sum += pow(x - expected_value, 4);

    return 1 / (n * pow(variance, 2)) * sum - 3;
}

double empirical_huber(int n, double x, vector<double> x_s)
{

    vector<double> f_s;
    int k = (int)trunc(log2((double)n)) + 1;
    double min_x = *min_element(begin(x_s), end(x_s));
    double max_x = *max_element(begin(x_s), end(x_s));
    double delta = (max_x - min_x) / (double)k;


    for (int i = 0; i < k; i++)
    {

        if (min_x + delta * i <= x && x < min_x + delta * (i + 1))
        {
            int n_i = count_if(x_s.begin(), x_s.end(), [i, k, min_x, max_x, delta](double x) { return i == k - 1 ? min_x + delta * (double)i <= x && x <= min_x + delta * (double)(i + 1) : min_x + delta * (double)i <= x && x < min_x + delta * (double)(i + 1); });
            return n_i / (n * delta);
        }
    }
}

vector<double> generate_sequence(int n, HuberDistribution* HB)
{
    vector<double> x_s;

    for (int i = 0; i < n; i++)
    {
        double x = calculate_x(HB);
        x_s.push_back(x);
    }

    return x_s;
}