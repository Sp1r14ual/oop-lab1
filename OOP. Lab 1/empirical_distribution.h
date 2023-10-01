#pragma once
#include "lib.h"

using namespace std;

double empirical_expected_value(int n, vector<double> x_s);

double empirical_variance(int n, vector<double> x_s);

double empirical_asymmetry(int n, vector<double> x_s);

double empirical_kurtosis(int n, vector<double> x_s);

double empirical_huber(int n, double x, vector<double> x_s);