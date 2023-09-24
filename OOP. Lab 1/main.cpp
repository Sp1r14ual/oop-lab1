#define CATCH_CONFIG_RUNNER
#define _USE_MATH_DEFINES

#include "catch.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <string>

using namespace std;

//распределение Хьюбера
double Huber(double x, double v, double K, double scale = 1., double shift = 0.)
{
    return (1. / (sqrt(2. * M_PI) * K) * (abs((x - shift) / scale) <= v ? exp(-pow((x - shift) / scale, 2.) / 2.) : exp(pow(v, 2.) / 2. - v * abs((x - shift) / scale))))  / scale;
}

//функция стандартного нормального распределения
double phi(double x)
{
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

//плотность стандартного нормального распределения
 double phi_lower(double x)
{
     return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}

//дисперсия
double huber_variance(double v, double K)
{
    return 1. + 2. * phi_lower(v) * (pow(v, 2.) + 2.) / (pow(v, 3.) * K);
}

//коэффициент эксцесса
double huber_kurtosis(double v, double K)
{
    return (3. * (2. * phi(v) - 1.) + 2. * phi_lower(v) * (24. / pow(v, 5.) + 24. / pow(v, 3.) + 12. / v + v)) / (pow(huber_variance(v, K), 2.) * K) - 3.;
}

double P(double v, double K)
{
    return (2. * phi(v) - 1.) / K;
}

double K (double v)
{
    return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}

double calculate_x(double v, double K, double scale = 1., double shift = 0.)
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

double empirical_expected_value(int n, vector<double> x_s)
{
    double sum = 0;
    for (double &x : x_s)
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

double empirical_asymetry(int n, vector<double> x_s)
{
    double expected_value = empirical_expected_value(n, x_s);
    double variance = empirical_variance(n, x_s);

    double sum = 0;
    for (double& x : x_s)
        sum += pow(x - expected_value, 3);

    return 1 / (n * pow(variance, 3/2)) * sum;
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

double mixture(double x, double v1, double v2, double K1, double K2, double p, double scale1 = 1., double shift1 = 0., double scale2 = 1., double shift2 = 0.)
{
    return (1 - p) * Huber(x, v1, K1, scale1, shift1) + p * Huber(x, v2, K2, scale2, shift2);
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


TEST_CASE("standard distribution, scale = 1, shift = 0")
{
    double v = 0.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(8.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.94).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.214).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.223).epsilon(0.01));

    v = 0.75;
    REQUIRE(huber_variance(v, K(v)) == Approx(3.71).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.75).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.405).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.296).epsilon(0.01));

    v = 1;
    REQUIRE(huber_variance(v, K(v)) == Approx(2.24).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.37).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.585).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.342).epsilon(0.01));

    v = 1.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.31).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(1.30).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.834).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.384).epsilon(0.01));

    v = 2;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.51).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.946).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.396).epsilon(0.01));

    v = 2.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.02).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.16).epsilon(0.1));
    REQUIRE(P(v, K(v)) == Approx(0.986).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.398).epsilon(0.01));

    v = 3;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.00).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.04).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.997).epsilon(0.01));
    REQUIRE(Huber(0., v, K(v)) == Approx(0.399).epsilon(0.01));
}

TEST_CASE("scaled and shifted distribution, scale = 2, shift = 5")
{
    double v = 0.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(8.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.94).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.214).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.223).epsilon(0.01));

    v = 0.75;
    REQUIRE(huber_variance(v, K(v)) == Approx(3.71).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.75).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.405).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.296).epsilon(0.01));

    v = 1;
    REQUIRE(huber_variance(v, K(v)) == Approx(2.24).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(2.37).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.585).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.342).epsilon(0.01));

    v = 1.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.31).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(1.30).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.834).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.384).epsilon(0.01));

    v = 2;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.51).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.946).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.396).epsilon(0.01));

    v = 2.5;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.02).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.16).epsilon(0.1));
    REQUIRE(P(v, K(v)) == Approx(0.986).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.398).epsilon(0.01));

    v = 3;
    REQUIRE(huber_variance(v, K(v)) == Approx(1.00).epsilon(0.01));
    REQUIRE(huber_kurtosis(v, K(v)) == Approx(0.04).epsilon(0.01));
    REQUIRE(P(v, K(v)) == Approx(0.997).epsilon(0.01));
    REQUIRE(Huber(0. * 2. + 5., v, K(v), 2., 5.) == Approx(0.399).epsilon(0.01));
}

TEST_CASE("mixture of distributions: shift1 = shift2 =/= 0, scale1 = scale2 = 2, v1 = v2, p - random")
{
    double p = 0.6;
    double shift1 = 8;
	double shift2 = 8;
    double scale1 = 2;
	double scale2 = 2;

    double v1 = 0.5;
	double v2 = 0.5;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));


    v1 = 0.75;
    v2 = 0.75;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));

    v1 = 1;
    v2 = 1;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));

    v1 = 1.5;
    v2 = 1.5;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));

    v1 = 2;
    v2 = 2;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));

    v1 = 2.5;
    v2 = 2.5;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));

    v1 = 3;
    v2 = 3;
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v1, K(v1), scale1, shift1)).epsilon(0.01));
    REQUIRE(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) == Approx(Huber(0., v2, K(v2), scale2, shift2)).epsilon(0.01));
}



TEST_CASE("mixture of distributions: shift1 =/= shift2, scale1, scale2, v1, v2 - random, p = 0.5")
{
    double shift1 = 1;
    double shift2 = 2;
    double scale1 = 1;
    double scale2 = 1;
    double p = 0.5;

    double v1 = 1;
    double v2 = 1;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1)) < 0.1);
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v2, K(v2), scale2, shift2)) < 0.1);

	shift1 = 5;
	shift2 = 10;
    scale1 = 2;
    scale2 = 4;

    v1 = 2;
    v2 = 3;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.1));
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v2, K(v2), scale2, shift2) < 0.1));
}

TEST_CASE("mixture of distributions: shift1 = shift2 = 0, scale1 = 1, scale2 = 3, v1 = v2, p = 0.5")
{
    double shift1 = 0;
    double shift2 = 0;
    double scale1 = 1;
    double scale2 = 3;
    double p = 0.5;

    double v1 = 0.5;
    double v2 = 0.5;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 0.75;
    v2 = 0.75;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 1;
    v2 = 1;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 1.5;
    v2 = 1.5;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 2;
    v2 = 2;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 2.5;
    v2 = 2.5;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));

    v1 = 3;
    v2 = 3;
    REQUIRE(abs(mixture(0., v1, v2, K(v1), K(v2), p, scale1, shift1, scale2, shift2) - Huber(0., v1, K(v1), scale1, shift1) < 0.001));
}



TEST_CASE("compare theoretical and empirical values")
{
    vector<double> v_params = { 0.5, 0.75, 1, 1.5, 2, 2.5, 3 };
    vector<double> x_s;
    int n = 10;

    double v = 0.5;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 0.75;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 1;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 1.5;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 2;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 2.5;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);

    v = 3;
    x_s = generate_sequence(n, v, K(v));
    REQUIRE(abs(huber_variance(v, K(v)) - empirical_variance(n, x_s)) < 5);
    REQUIRE(abs(huber_kurtosis(v, K(v)) - empirical_kurtosis(n, x_s)) < 5);
}



TEST_CASE("calculate theoretical and empirical distributions for analysis, n - small")
{
    vector<double> v_params = { 0.5, 0.75, 1, 1.5, 2, 2.5, 3 };
    vector<double> x_s;
    int n = 100;

    ofstream xs;
    ofstream fs_theoretical;
    ofstream fs_empirical;

    for (double& v : v_params)
    {
        xs.open("xs v=" + to_string(v) + ".txt");
        fs_theoretical.open("fs_theoretical v=" + to_string(v) + ".txt");
        fs_empirical.open("fs_empirical v=" + to_string(v) + ".txt");

        x_s = generate_sequence(n, v, K(v));
        sort(x_s.begin(), x_s.end());

        for (double& x : x_s)
        {
            xs << x << endl;

            double f_theoretical = Huber(x, v, K(v));
            fs_theoretical << f_theoretical << endl;

            double f_empirical = empirical_huber(n, x, x_s);
            fs_empirical << f_empirical << endl;
        }

        xs.close();
        fs_theoretical.close();
        fs_empirical.close();
    }
}


TEST_CASE("compare 2 empirical distributions")
{
    int n = 10;
    
    double v = 0.5;
    vector<double> x_s1 = generate_sequence(n, v, K(v));
    vector<double> x_s2 = generate_sequence(n, v, K(v));

    double variance_theoretical = huber_variance(v, K(v));
    double kurtosis_theoretical = huber_kurtosis(v, K(v));
    double variance_empirical1 = empirical_variance(n, x_s1);
    double kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    double variance_empirical2 = empirical_variance(n, x_s2);
    double kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 0.75;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 1;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 1.5;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 2;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 2.5;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

    v = 3;
    x_s1 = generate_sequence(n, v, K(v));
    x_s2 = generate_sequence(n, v, K(v));

    variance_theoretical = huber_variance(v, K(v));
    kurtosis_theoretical = huber_kurtosis(v, K(v));
    variance_empirical1 = empirical_variance(n, x_s1);
    kurtosis_empirical1 = empirical_kurtosis(n, x_s1);
    variance_empirical2 = empirical_variance(n, x_s2);
    kurtosis_empirical2 = empirical_kurtosis(n, x_s2);

    REQUIRE(abs(variance_theoretical - variance_empirical1) < 10);
    REQUIRE(abs(variance_theoretical = variance_empirical2) < 10);
    REQUIRE(abs(variance_empirical1 - variance_empirical2) < 10);
    REQUIRE(abs(kurtosis_theoretical - kurtosis_empirical1) < 10);
    REQUIRE(abs(kurtosis_theoretical = kurtosis_empirical2) < 10);
    REQUIRE(abs(kurtosis_empirical1 - kurtosis_empirical2) < 10);

}



int main(int argc, char** argv)
{
    int result = Catch::Session().run(argc, argv);
    return result;

}