#include "basic_io.h"
#include <iostream>
#include <fstream>
#include <iomanip>

tuple<int, int, int, vector<int>, vector<int> > load_data(string file_name) {
	int n, m, T;
	vector<int> mu;
	vector<int> sigma;

	ifstream fin(file_name);

	fin >> m >> n >> T;

	for (int i = 0; i < n; i++) {
		int m, s;
		fin >> m >> s;
		mu.push_back(m);
		sigma.push_back(s);
	}

	fin.close();

	return make_tuple(m, n, T, mu, sigma);
}

void write_ans(string file_name, int index, double opt, double t, vector<int> x) {

	ofstream fout(file_name, ios::app);
	fout << setiosflags(ios::fixed);
	fout << setprecision(10);

	fout << index << " " << opt << " " << t << " ";
	for (int i = 0; i < x.size(); i++)
		fout << x[i] << " ";
	fout << endl;

	fout.close();
}