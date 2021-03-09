#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <tuple>
#include <string>

#include "basic_io.h"
#include "std_config_ip.h"
#include "approx_fpt.h"

using namespace std;


int main(int argc, char *argv[]) {


	string file_prefix = "test";
	int s_index = 0;
	int tot = 1;

	string alg = "config_ip";
	double approx_fpt_parameter_eta = 10;

	int max_thread = 8;


	if (argc >= 4) {
		file_prefix = argv[1];

		s_index = atoi(argv[2]);

		tot = atoi(argv[3]);
	}


	if (argc >= 5) {
		alg = argv[4];

		if (alg == "approx_fpt" && argc >= 6)
			approx_fpt_parameter_eta = 1 + atof(argv[5]);

		if (alg == "approx_fpt" && argc >= 7)
			max_thread = atoi(argv[6]);
	}


	for (int i = s_index; i < tot; i++) {
		int m, n, T;
		vector<int> mu;
		vector<int> sigma;

		cout << file_prefix + "//" + to_string(i) + ".txt" << endl;

		tuple<int, int, int, vector<int>, vector<int> > data = load_data(file_prefix + "//" + to_string(i) + ".txt");

		m = get<0>(data);
		n = get<1>(data);
		T = get<2>(data);
		mu = get<3>(data);
		sigma = get<4>(data);




		tuple<double, double, vector<int>> ans;
		
		////////////////////////////////////////////////////////////////////////
		

		if (alg == "config_ip") {
			ans = mip_config(m, n, T, mu, sigma, true);
			write_ans(file_prefix + "//opt_std_config.txt", i, get<0>(ans), get<1>(ans), get<2>(ans));
		}

		if (alg == "approx_fpt") {
			ans = FPT(m, n, T, approx_fpt_parameter_eta, mu, sigma, 16);
			write_ans(file_prefix + "//approx_fpt_" + to_string(approx_fpt_parameter_eta - 1) + ".txt", i, get<0>(ans), get<1>(ans), get<2>(ans));
		}


		////////////////////////////////////////////////////////////////////////
	}


	//system("pause");

}