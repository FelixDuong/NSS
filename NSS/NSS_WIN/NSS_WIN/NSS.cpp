
#include "pch.h"
#include <iostream>

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "optimization.h"
#include "interpolation.h"

using namespace alglib;
using namespace std;

real_1d_array T = "[0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20]";
real_1d_array Y = "[3.12,	3.22,	3.41	,3.77	,4.08,	4.59	,4.98,	5.4,	5.83,	6.18]";

double f_NSS(const real_1d_array &x, double T) {

	double y = x(0) + x(1) * (1 - exp(-T / x(4))) * 1 / (T / x(4))
			+ x(2) * ((1 - exp(-T / x(4))) * 1 / (T / x(4)) - exp(-T / x(4)))
			+ x(3) * ((1 - exp(-T / x(5))) * 1 / (T / x(5)) - exp(-T / x(5)));
	return y;
}

void NSS(const real_1d_array &x, real_1d_array &fi, void *ptr) {

//	double b0 = x(0);
//	double b1 = x(1);
//	double b2 = x(2);
//	double b3 = x(3);
//	double i1 = x(4);
//	double i2 = x(5);
	for (int i = 0; i < 10; ++i) {
		fi[i] = pow(
				x(0) + x(1) * (1 - exp(-T(i) / x(4))) * 1 / (T(i) / x(4))
						+ x(2)
								* ((1 - exp(-T(i) / x(4))) * 1 / (T(i) / x(4))
										- exp(-T(i) / x(4)))
						+ x(3)
								* ((1 - exp(-T(i) / x(5))) * 1 / (T(i) / x(5))
										- exp(-T(i) / x(5))) - Y(i), 2);
	}
}

void NSS_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac,
		void *ptr) {

//	double b0 = x(0);
//	double b1 = x(1);
//	double b2 = x(2);
//	double b3 = x(3);
//	double i1 = x(4);
//	double i2 = x(5);
	for (int i = 0; i < 10; ++i) {

		fi[i] = pow(
				x(0) + x(1) * (1 - exp(-T(i) / x(4))) * 1 / (T(i) / x(4))
						+ x(2)
								* ((1 - exp(-T(i) / x(4))) * 1 / (T(i) / x(4))
										- exp(-T(i) / x(4)))
						+ x(3)
								* ((1 - exp(-T(i) / x(5))) * 1 / (T(i) / x(5))
										- exp(-T(i) / x(5))) - Y(i), 2);
//
		jac[i][0] = 2 * x(0) - 2 * Y(i)
				- 2 * x(2)
						* (exp(-T(i) / x(4))
								+ (x(4) * (exp(-T(i) / x(4)) - 1)) / T(i))
				- 2 * x(3)
						* (exp(-T(i) / x(5))
								+ (x(5) * (exp(-T(i) / x(5)) - 1)) / T(i))
				- (2 * x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i);
		jac[i][1] = (2 * x(4) * (exp(-T(i) / x(4)) - 1)
				* (Y(i) - x(0)
						+ x(2)
								* (exp(-T(i) / x(4))
										+ (x(4) * (exp(-T(i) / x(4)) - 1))
												/ T(i))
						+ x(3)
								* (exp(-T(i) / x(5))
										+ (x(5) * (exp(-T(i) / x(5)) - 1))
												/ T(i))
						+ (x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i)))
				/ T(i);
		jac[i][2] = 2
				* (exp(-T(i) / x(4)) + (x(4) * (exp(-T(i) / x(4)) - 1)) / T(i))
				* (Y(i) - x(0)
						+ x(2)
								* (exp(-T(i) / x(4))
										+ (x(4) * (exp(-T(i) / x(4)) - 1))
												/ T(i))
						+ x(3)
								* (exp(-T(i) / x(5))
										+ (x(5) * (exp(-T(i) / x(5)) - 1))
												/ T(i))
						+ (x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i));
		jac[i][3] = 2
				* (exp(-T(i) / x(5)) + (x(5) * (exp(-T(i) / x(5)) - 1)) / T(i))
				* (Y(i) - x(0)
						+ x(2)
								* (exp(-T(i) / x(4))
										+ (x(4) * (exp(-T(i) / x(4)) - 1))
												/ T(i))
						+ x(3)
								* (exp(-T(i) / x(5))
										+ (x(5) * (exp(-T(i) / x(5)) - 1))
												/ T(i))
						+ (x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i));
		jac[i][4] = 2
				* (x(2)
						* (exp(-T(i) / x(4)) / x(4)
								+ (exp(-T(i) / x(4)) - 1) / T(i)
								+ (T(i) * exp(-T(i) / x(4))) / x(4) * x(4))
						+ (x(1) * exp(-T(i) / x(4))) / x(4)
						+ (x(1) * (exp(-T(i) / x(4)) - 1)) / T(i))
				* (Y(i) - x(0)
						+ x(2)
								* (exp(-T(i) / x(4))
										+ (x(4) * (exp(-T(i) / x(4)) - 1))
												/ T(i))
						+ x(3)
								* (exp(-T(i) / x(5))
										+ (x(5) * (exp(-T(i) / x(5)) - 1))
												/ T(i))
						+ (x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i));
		jac[i][5] = 2 * x(3)
				* (exp(-T(i) / x(5)) / x(5) + (exp(-T(i) / x(5)) - 1) / T(i)
						+ (T(i) * exp(-T(i) / x(5))) / x(5) * x(5))
				* (Y(i) - x(0)
						+ x(2)
								* (exp(-T(i) / x(4))
										+ (x(4) * (exp(-T(i) / x(4)) - 1))
												/ T(i))
						+ x(3)
								* (exp(-T(i) / x(5))
										+ (x(5) * (exp(-T(i) / x(5)) - 1))
												/ T(i))
						+ (x(1) * x(4) * (exp(-T(i) / x(4)) - 1)) / T(i));
	}
}

int main(int argc, char **argv) {
	clock_t tStart = clock();
	real_1d_array x =
			"[8.018777397,	-4.878138602,	0.009703287	,-6.223070427,	2.317207769	,3.956212859]";
	real_1d_array bndl = "[0,-inf,-inf,-inf,0.5,0.5]";

	real_1d_array bndu = "[100,+inf,+inf,+inf,+20,+20]";

	spline1dinterpolant s;

	double epsx = 0.00000001;
	ae_int_t maxits = 0;
	minlmstate state;
	minlmreport rep;

	// Load data and compute
	ifstream ip("Raw.csv");
	if (!ip.is_open())
		std::cout << "ERROR: File Open" << '\n';

	string Date;
	string Y3M;
	string Y6M;
	string Y1Y;
	string Y2Y;
	string Y3Y;
	string Y5Y;
	string Y7Y;
	string Y10Y;
	string Y15Y;
	string Y20Y;
	std::string filename = "NSS_Results.csv";
	ofstream fs;
	fs.open(filename);
	while (ip.good()) {

		getline(ip, Date, ',');
		getline(ip, Y3M, ',');
		getline(ip, Y6M, ',');
		getline(ip, Y1Y, ',');
		getline(ip, Y2Y, ',');
		getline(ip, Y3Y, ',');
		getline(ip, Y5Y, ',');
		getline(ip, Y7Y, ',');
		getline(ip, Y10Y, ',');
		getline(ip, Y15Y, ',');
		getline(ip, Y20Y, '\n');

		Y[0] = std::stod(Y3M);
		Y[1] = std::stod(Y6M);
		Y[2] = std::stod(Y1Y);
		Y[3] = std::stod(Y2Y);
		Y[4] = std::stod(Y3Y);
		Y[5] = std::stod(Y5Y);
		Y[6] = std::stod(Y7Y);
		Y[7] = std::stod(Y10Y);
		Y[8] = std::stod(Y15Y);
		Y[9] = std::stod(Y20Y);

		bndl(0) = Y(0);
		minlmcreatev(10, x, 0.00000001, state);
		//minlmcreatevj(10, x, state);

		minlmsetcond(state, epsx, maxits);

		minlmsetbc(state, bndl, bndu);
//		spline1dbuildcubic(T, Y, s);
//		alglib::minlmoptimize(state, NSS);
		alglib::minlmoptimize(state, NSS, NSS_jac);
		minlmresults(state, x, rep);

		printf("%s  %s\n", Date.c_str(), x.tostring(2).c_str());
//		printf("%d\n", int(rep.terminationtype));

// Data generates via Cubic spline interpolation
		fs << Date.c_str() <<",";
		for (int i = 0; i < 360; ++i) {
			fs << f_NSS(x, double((i+1) * 0.083333)) <<",";
		}
		fs << endl;
		//				<< spline1dcalc(s, double((i + 1) * 0.05)) << endl;

	}
	fs.close();
	ip.close();
	//
	printf("Time completed in : %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	return 0;
}

