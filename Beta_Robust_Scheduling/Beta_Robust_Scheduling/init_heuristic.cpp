#include <tuple>
#include <vector>
#include <iostream>

#include "distribution.h"

using namespace std;


tuple<double, vector<int>> init_heuristic(int m, int n, int T, vector<int> mu, vector<int> sigma) {

	double ret = 1;
	vector<int> ans(n);

	vector<int> rec(n);

	for (int i = 0; i < n; i++)
		rec[i] = i;

	for (int i = 0; i < mu.size(); i++)
		for (int j = 0; j + 1 < mu.size(); j++)
			if (mu[j] < mu[j + 1]) {
				int tmp = mu[j];
				mu[j] = mu[j + 1];
				mu[j + 1] = tmp;

				tmp = sigma[j];
				sigma[j] = sigma[j + 1];
				sigma[j + 1] = tmp;

				tmp = rec[j];
				rec[j] = rec[j + 1];
				rec[j + 1] = tmp;
			}


	int *sum_mu = new int[m];
	int *sum_sigma = new int[m];


	for (int j = 0; j < m; j++)
		sum_mu[j] = sum_sigma[j] = 0;


	for (int i = 0; i < mu.size(); i++) {

		int index = 0;

		for (int j = 0; j < m; j++)
			if (sum_mu[index] > sum_mu[j])
				index = j;

		sum_mu[index] += mu[i];
		sum_sigma[index] += sigma[i];
		ans[rec[i]] = index;
	}

	for (int j = 0; j < m; j++)
		ret *= phi((T - sum_mu[j]) / sqrt(sum_sigma[j]));

	delete[] sum_mu;
	delete[] sum_sigma;

	cout << "init_heuristic = " << ret << endl;


	return make_tuple(ret, ans);
}

