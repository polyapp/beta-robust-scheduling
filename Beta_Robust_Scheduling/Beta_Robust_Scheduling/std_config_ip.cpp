#include "std_config_ip.h"
#include "distribution.h"

#include "gurobi_c++.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <tuple>


bool item_inside_config(int i, int j) {
	return (j&(1 << i)) > 0;
}


double log_phi_config(int n, int j, int T, vector<int> mu, vector<int> sigma, bool flag = false) {
	int sum_mu = 0;
	int sum_sigma = 0;
	for (int i = 0; i < n; i++)
		if (item_inside_config(i, j)) {
			sum_mu += mu[i];
			sum_sigma += sigma[i];
		}
	if (flag) {
		cout << j << " " << sum_mu << " " << sum_sigma << endl;
		cout << " phi " << (T - sum_mu) / sqrt(sum_sigma) << " = " << phi((T - sum_mu) / sqrt(sum_sigma)) << endl;
	}

	double ret = log(phi((T - sum_mu) / sqrt(sum_sigma)));

	if (ret < -1e8) ret = -1e8;
	return ret;
}


tuple<double, double, vector<int>> mip_config(int m, int n, int T, vector<int> mu, vector<int> sigma, bool is_integer) {

	cout <<"Standard_Configuration IP" << endl << endl;

	cout <<"m = "<< m << "| n = " << n << "| T = " << T << endl;

	size_t t0 = clock();

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//	model.set(GRB_IntParam_Presolve, 0);
	//	model.set(GRB_IntParam_Threads, 1);

	double opt;
	vector<int> ans(n);
	GRBVar *x = new GRBVar[int(1 << n)];

	double *cof = new double[int(1 << n)];

	cof[0] = 0;

	try {

		for (int j = 0; j < (1 << n); j++) {
			if (is_integer)
				x[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			else
				x[j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
			

		GRBLinExpr obj = 0;

		for (int j = 1; j < (1 << n); j++)
			cof[j] = log_phi_config(n, j, T, mu, sigma);
		obj.addTerms(cof, x, (1 << n));
		//for (int j = 1; j < (1 << n); j++)
		//	obj += log_phi_config(n, j, T, mu, sigma)* x[j];
		model.setObjective(obj, GRB_MAXIMIZE);


		GRBLinExpr cons = 0;
		for (int j = 1; j < (1 << n); j++)
			cof[j] = 1.0;
		cons.addTerms(cof, x, (1 << n));
		//for (int j = 1; j < (1 << n); j++)
		//	cons += x[j];
		model.addConstr(cons == m);

		for (int i = 0; i < n; i++) {
			GRBLinExpr exp = 0;
			for (int j = 1; j < (1 << n); j++)
				if (item_inside_config(i, j))
					cof[j] = 1.0;
				else
					cof[j] = 0;
			exp.addTerms(cof, x, (1 << n));
			//for (int j = 1; j < (1 << n); j++)
			//	if (item_inside_config(i, j))
			//		exp += x[j];
			model.addConstr(exp == 1);
		}



		model.optimize();


		int cnt = 0;
		for (int i = 0; i < n; i++)
			ans[i] = -1;

		for (int j = 1; j < (1 << n); j++)
			if (fabs(x[j].get(GRB_DoubleAttr_X) - 1) < 1e-4) {

				//	cout << j << " " << x[j].get(GRB_DoubleAttr_X) << endl;

				for (int i = 0; i < n; i++)
					if (item_inside_config(i, j))
						ans[i] = cnt;
				cnt++;
			}

		opt = exp(model.get(GRB_DoubleAttr_ObjVal));

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}


	delete[] cof;
	delete[] x;

	return make_tuple(opt, (clock() - t0) / 1000.0, ans);
}
