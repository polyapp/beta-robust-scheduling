#include <vector>
#include <tuple>
#include <cmath>
#include "distribution.h"
#include <iostream>
using namespace std;



double phi(double x) {
	double ret = 0;

	double  a = 7.0 * exp(-0.5*x*x);

	double  b = 16 * exp(-x*x*(2 - sqrt(2)));

	double  c = (7 + 0.25*PI*x*x)*exp(-x*x);

	ret = 0.5 + 0.5*sqrt(1 - 1.0 / 30 * (a + b + c));

	if (x < 0) ret = 1 - ret;

	return ret;
}


double inverse_phi(double y) {
	double l = -5;
	double u = 5;

	while (u - l > 1e-8) {
		if (phi((u + l) / 2) > y)
			u = (u + l) / 2;
		else
			l = (u + l) / 2;
	}

	return l;
}



double inverse_logphi(double y) {
	double l = -5;
	double u = 5;

	while (u - l > 1e-8) {
		if (log(phi((u + l) / 2)) > y)
			u = (u + l) / 2;
		else
			l = (u + l) / 2;
	}

	return l;
}