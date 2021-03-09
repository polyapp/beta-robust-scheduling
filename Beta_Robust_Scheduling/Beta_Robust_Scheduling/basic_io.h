#pragma once
#include <tuple>
#include <vector>
#include <string>
using namespace std;

tuple<int, int, int, vector<int>, vector<int> > load_data(string file_name);

void write_ans(string file_name, int index, double opt, double t, vector<int> x);