#include "catch.hpp"
#include "huber_distribution.h"
#include "empirical_distribution.h"
#include "mixture_distribution.h"

using namespace std;

TEST_CASE("standard distribution, scale = 1, shift = 0")
{
    HuberDistribution* HB = new HuberDistribution();
    HB->v = 0.5;
    HB->K = K(HB->v);
    HB->scale = 1.;
    HB->shift = 0.;
    REQUIRE(huber_variance(HB) == Approx(8.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(2.94).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.214).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.223).epsilon(0.01));

    HB->v = 0.75;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(3.71).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(2.75).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.405).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.296).epsilon(0.01));

    HB->v = 1;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(2.24).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(2.37).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.585).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.342).epsilon(0.01));

    HB->v = 1.5;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(1.31).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(1.30).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.834).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.384).epsilon(0.01));

    HB->v = 2;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(1.08).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(0.51).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.946).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.396).epsilon(0.01));

    HB->v = 2.5;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(1.02).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(0.16).epsilon(0.1));
    REQUIRE(P(HB) == Approx(0.986).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.398).epsilon(0.01));

    HB->v = 3;
    HB->K = K(HB->v);
    REQUIRE(huber_variance(HB) == Approx(1.00).epsilon(0.01));
    REQUIRE(huber_kurtosis(HB) == Approx(0.04).epsilon(0.01));
    REQUIRE(P(HB) == Approx(0.997).epsilon(0.01));
    REQUIRE(Huber(0., HB) == Approx(0.399).epsilon(0.01));
}
