#pragma once
#include <tuple>
#include <vector>
using namespace std;

tuple<double, double, vector<int>> mip_config(int m, int n, int T, vector<int> mu, vector<int> sigma, bool is_integer);