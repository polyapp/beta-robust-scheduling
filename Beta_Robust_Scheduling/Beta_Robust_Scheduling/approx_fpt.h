#pragma once
#include <tuple>
#include <vector>
using namespace std;


tuple<double, double, vector<int>> FPT(int m, int n, int T, double eta, vector<int> mu, vector<int> sigma, int max_thread);