#include "approx_fpt.h"
#include "distribution.h"

#include "gurobi_c++.h"

#include <map>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "init_heuristic.h"

int sum_vec(vector<int> x) {
	int ret = 0;
	for (int i = 0; i < x.size(); i++)
		ret += x[i];
	return ret;
}

vector<tuple<int, int> > construct_sigma_segment(double eta, vector<int> sigma) {

	int sum_sigma = sum_vec(sigma);


	vector<double> s, f;
	double pts = 1;

	while (pts <= sum_sigma) {
		s.push_back(pts);
		f.push_back(pts*eta);


		pts *= eta;
	}

	vector<tuple<int, int>> seg;

	for (int i = 0; i < s.size(); i++)
		if (ceil(s[i]) <= floor(f[i]))
			seg.push_back(make_tuple(int(ceil(s[i])), int(floor(f[i]))));

	return seg;

}

vector<tuple<int, int>> construct_mu_segment(double eta, int T, int tx, vector<int> mu) {

	int sum_mu = sum_vec(mu);


	vector<double> s, f;
	s.push_back(T);
	f.push_back(T);

	s.push_back(T + 1);
	f.push_back(T + 1);

	double dtx = tx > 1 ? tx : tx*eta;

	double pt = sqrt(dtx);

	s.push_back(T + 1);
	f.push_back(T + pt);

	while (T + pt <= sum_mu) {
		s.push_back(T + pt);
		f.push_back(T + pt*eta);
		pt *= eta;
	}

	pt = sqrt(dtx);

	s.push_back(T - 1);
	f.push_back(T - 1);

	s.push_back(T - pt);
	f.push_back(T - 1);



	while (T - pt >= 0) {
		f.push_back(T - pt);
		s.push_back(T - pt*eta);
		pt *= eta;
	}

	vector<tuple<int, int>> seg;

	vector<int> sigma_l, sigma_u;

	for (int i = 0; i < s.size(); i++)
		if (ceil(s[i]) <= floor(f[i]))
			seg.push_back(make_tuple(int(ceil(s[i])), int(floor(f[i]))));

	for (int i = 0; i < seg.size(); i++)
		for (int j = 0; j < seg.size() - 1; j++)
			if (abs(get<0>(seg[j]) - T) > abs(get<0>(seg[j + 1]) - T)) {

				tuple<int, int> tmp = seg[j];
				seg[j] = seg[j + 1];
				seg[j + 1] = tmp;

			}

	return seg;

}

tuple<vector<tuple<int, int>>, map<int, vector<tuple<int, int>>>> construct_segment(double eta, int T, vector<int> sigma, vector<int> mu) {

	vector<tuple<int, int>> seg_sigma = construct_sigma_segment(eta, sigma);

	map<int, vector<tuple<int, int>>> seg_mu;

	for (int i = 0; i < seg_sigma.size(); i++)
		seg_mu[i] = construct_mu_segment(eta, T, get<1>(seg_sigma[i]), mu);

	return make_tuple(seg_sigma, seg_mu);

}



double max_phi(int T, tuple<int, int> m, tuple<int, int> s) {
	double ret = phi((T - get<0>(m)) / sqrt(get<0>(s)));

	if (ret < phi((T - get<0>(m)) / sqrt(get<1>(s))))
		ret = phi((T - get<0>(m)) / sqrt(get<1>(s)));

	if (ret < phi((T - get<1>(m)) / sqrt(get<0>(s))))
		ret = phi((T - get<1>(m)) / sqrt(get<0>(s)));

	if (ret < phi((T - get<1>(m)) / sqrt(get<1>(s))))
		ret = phi((T - get<1>(m)) / sqrt(get<1>(s)));

	return ret;
}

vector<tuple<int, int, vector<int>>> compress_data(vector<int> mu, vector<int> sigma) {

	vector<tuple<int, int, vector<int>>> ret;

	for (int i = 0; i < mu.size(); i++) {

		int index = 0;

		while (index < ret.size() && (mu[i] != get<0>(ret[index]) || sigma[i] != get<1>(ret[index])))
			index++;

		if (index == ret.size()) {
			vector<int> tmp;
			tmp.push_back(i);
			ret.push_back(make_tuple(mu[i], sigma[i], tmp));
		}
		else {
			vector<int> tmp = get<2>(ret[index]);
			tmp.push_back(i);

			int t0 = get<0>(ret[index]);
			int t1 = get<1>(ret[index]);

			ret[index] = make_tuple(t0, t1, tmp);
		}
	}

	return ret;
}




tuple<double, vector<int> > solve(int machine_num, int n, int T, vector<tuple<int, int>> &guess_mu, vector<tuple<int, int>> &guess_sigma, vector<tuple<int, int, vector<int>>> &item_data, GRBEnv env) {

	//GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_Presolve, 0);
	model.set(GRB_IntParam_OutputFlag, 0);

	double opt = -1;
	vector<int> ans(n);

	GRBVar **x = new GRBVar*[machine_num];

	try {

		for (int i = 0; i < machine_num; i++)
			x[i] = new GRBVar[item_data.size()];

		for (int i = 0; i < machine_num; i++)
			for (int j = 0; j < item_data.size(); j++)
				x[i][j] = model.addVar(0, int(get<2>(item_data[j]).size()), 0.0, GRB_INTEGER);

		GRBLinExpr obj = 0;
		model.setObjective(obj, GRB_MAXIMIZE);


		for (int i = 0; i < machine_num; i++)
			for (int j = 0; j < item_data.size(); j++) {
				GRBLinExpr cons = x[i][j];
				model.addConstr(cons >= 0);
				model.addConstr(cons <= int(get<2>(item_data[j]).size()));
			}



		for (int j = 0; j < item_data.size(); j++) {
			GRBLinExpr cons = 0;
			for (int i = 0; i < machine_num; i++)
				cons += x[i][j];
			model.addConstr(cons == int(get<2>(item_data[j]).size()));
		}

		for (int i = 0; i < machine_num; i++) {
			GRBLinExpr cons_mu = 0;
			for (int j = 0; j < item_data.size(); j++)
				cons_mu += x[i][j] * get<0>(item_data[j]);

			model.addConstr(cons_mu >= get<0>(guess_mu[i]));
			model.addConstr(cons_mu <= get<1>(guess_mu[i]));
		}

		for (int i = 0; i < machine_num; i++) {
			GRBLinExpr cons_sigma = 0;
			for (int j = 0; j < item_data.size(); j++)
				cons_sigma += x[i][j] * get<1>(item_data[j]);

			model.addConstr(cons_sigma >= get<0>(guess_sigma[i]));
			model.addConstr(cons_sigma <= get<1>(guess_sigma[i]));
		}


		model.optimize();



		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {



			int *pt = new int[item_data.size()];

			for (int j = 0; j < item_data.size(); j++)
				pt[j] = 0;

			opt = 1.0;

			for (int i = 0; i < machine_num; i++) {
				int sum_mu = 0;
				int sum_sigma = 0;

				for (int j = 0; j < item_data.size(); j++) {

					int k = int(round(x[i][j].get(GRB_DoubleAttr_X)));

					sum_mu += get<0>(item_data[j])*k;
					sum_sigma += get<1>(item_data[j])*k;

					for (int t = 0; t < k; t++)
						ans[get<2>(item_data[j])[pt[j] + t]] = i;
					pt[j] += k;

				}
				opt *= phi((T - sum_mu) / sqrt(sum_sigma));
			}



			delete[] pt;

		}

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}


	for (int i = 0; i < machine_num; i++)
		delete[] x[i];
	delete[] x;

	return make_tuple(opt, ans);
}


void update(vector<tuple<vector<tuple<int, int>>, vector<tuple<int, int>> >> &to_d_list, tuple<double, vector<int> > &ans, int max_thread, int m, int n, int T, vector<tuple<int, int, vector<int>>> &item_data, GRBEnv &env) {

	tuple<double, vector<int>> *tmp = new tuple<double, vector<int>>[max_thread];

#pragma omp parallel for
	for (int thread = 0; thread < to_d_list.size(); thread++)
		tmp[thread] = solve(m, n, T, get<0>(to_d_list[thread]), get<1>(to_d_list[thread]), item_data, env);

	for (int thread = 0; thread < to_d_list.size(); thread++)
		if (get<0>(ans) < get<0>(tmp[thread])) {
			ans = tmp[thread];
			cout << get<0>(ans) << endl;
		}
	to_d_list.clear();

	delete[] tmp;
}

void dfs(vector<tuple<vector<tuple<int, int>>, vector<tuple<int, int>> >> &to_d_list, int i, int machine_num, int n, int T, double cur, int mu_l, int mu_u, int sigma_l, int sigma_u, int sm, int ss,
	map<tuple<int, int>, double> &m_p, tuple<double, vector<int> > &ans,
	vector<tuple<int, int>> &guess_mu, vector<tuple<int, int>> &guess_sigma, map<int, vector<tuple<int, int>>> &seg_mu, vector<tuple<int, int>> &seg_sigma,
	vector<tuple<int, int, vector<int>>> &item_data, GRBEnv &env, int max_thread, size_t t0) {

	///////////////////////////////////
	if (cur <= get<0>(ans))  return;

	if (sm < mu_l || ss < sigma_l)  return;

	if (i +2 < machine_num )
		if (get<1>(guess_sigma[i+1]) > get<1>(guess_sigma[i + 2]))
			return;
	///////////////////////////////////


	if (i < 0) {

		///////////////////////////////////
		if (sm > mu_u || ss > sigma_u) return;

		if (clock() - t0 > 300000) return;
		///////////////////////////////////

		tuple<vector<tuple<int, int>>, vector<tuple<int, int>> > element = make_tuple(guess_mu, guess_sigma);

		to_d_list.push_back(element);


		if (to_d_list.size() == max_thread)
			update(to_d_list, ans, max_thread, machine_num, n, T, item_data, env);

	}
	else {
		for (int ps = 0; ps < seg_sigma.size(); ps++)
			for (int pm = 0; pm < seg_mu[ps].size(); pm++) {
				guess_mu[i] = seg_mu[ps][pm];
				guess_sigma[i] = seg_sigma[ps];

				dfs(to_d_list, i - 1, machine_num, n, T, cur * m_p[make_tuple(ps, pm)],
					mu_l + get<0>(guess_mu[i]), mu_u + get<1>(guess_mu[i]), sigma_l + get<0>(guess_sigma[i]), sigma_u + get<1>(guess_sigma[i]),
					sm, ss, m_p, ans, guess_mu, guess_sigma, seg_mu, seg_sigma, item_data, env, max_thread, t0);
			}


	}

}






tuple<double, double, vector<int>> FPT(int m, int n, int T, double eta, vector<int> mu, vector<int> sigma, int max_thread) {

	vector<tuple<int, int>> guess_mu(m);
	vector<tuple<int, int>> guess_sigma(m);


	tuple<vector<tuple<int, int>>, map<int, vector<tuple<int, int>>>> seg = construct_segment(eta, T, sigma, mu);

	vector<tuple<int, int>> seg_sigma = get<0>(seg);

	map<int, vector<tuple<int, int>>> seg_mu = get<1>(seg);

	map<tuple<int, int>, double> m_p;

	for (int i = 0; i < seg_sigma.size(); i++)
		for (int j = 0; j < seg_mu[i].size(); j++)
			m_p[make_tuple(i, j)] = max_phi(T, seg_mu[i][j], seg_sigma[i]);


	vector<tuple<int, int, vector<int>>> item_data = compress_data(mu, sigma);

	cout << "sigma_pieces# = " << seg_sigma.size() << " | ";
	cout << "mu_pieces# = " << seg_mu[0].size() << " | ";
	cout << endl;

	//cout << seg_sigma.size() << " " << seg_mu[0].size() << " " << pow(seg_sigma.size()*seg_mu[0].size(), m) << endl;


	size_t t0 = clock();

	tuple<double, vector<int>> ans = init_heuristic(m, n, T, mu, sigma);

	GRBEnv env = GRBEnv();
	vector<tuple<vector<tuple<int, int>>, vector<tuple<int, int>> >> to_d_list;

	cout << setiosflags(ios::fixed);
	cout << setprecision(10);

	dfs(to_d_list, m - 1, m, n, T, 1.0, 0, 0, 0, 0, sum_vec(mu), sum_vec(sigma), m_p,
		ans, guess_mu, guess_sigma, seg_mu, seg_sigma, item_data, env, max_thread, t0);

	update(to_d_list, ans, max_thread, m, n, T, item_data, env);


	return make_tuple(get<0>(ans), (clock() - t0) / 1000.0, get<1>(ans));

}
