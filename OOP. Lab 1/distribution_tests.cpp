#include "catch.hpp"
#include "huber_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.h"

using namespace std;
/*
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

//--------------------------------------------------------------------------------------------------------
// ֲûגמה ג פאיכ - מענופאךעמנטעü
//--------------------------------------------------------------------------------------------------------
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
*/