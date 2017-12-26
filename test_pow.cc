#include <cmath>
// #include <gsl/gsl_sf_gamma.h>
#include "gamma.cc"
#include <iostream>
#include <ctime>
#include <omp.h>
#include <iomanip>
#include <vector>
using namespace std;



vector< vector<double> > pow_val;
#define POW_BASE_LENGTH 100
#define POW_EXP_LENGTH 110000
void PowInit() {
	pow_val.resize(POW_BASE_LENGTH, vector<double>(POW_EXP_LENGTH));
	for (int i = 0; i < POW_BASE_LENGTH; ++i)
		for (int j = 0; j < POW_EXP_LENGTH; ++j)
			pow_val[i][j] = pow((double)i / POW_BASE_LENGTH, (double)j / POW_EXP_LENGTH * 11 - 1);
}

double GetPow(double base, double exponent) {
	if (exponent >= 10)
		return pow(base, exponent);
	else
		return pow_val[(int)(base * POW_BASE_LENGTH)][(int)((exponent + 1) / 11 * POW_BASE_LENGTH)];
}



int main(int argc, char const *argv[])
{
	double c = 0;
	int begin = clock();
	for(int i = 1; i < 100; ++i)
		for(int j = 0; j < 110000; ++j)
			c += pow(i / 100.0, j / 10000.0 - 1);
	cout << "pow:\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// c = 0;
	// begin = clock();
	// for(int i = 1; i < 100; ++i)
	// 	for(int j = 0; j < 110000; ++j)
	// 		c += __builtin_tgamma(i / 100.0, j / 11000.0 - 1);
	// cout << "__builtin:\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// c = 0;
	// begin = clock();
	// for(int i = 1; i < 100; ++i)
	// 	for(int j = 0; j < 110000; ++j)
	// 		c += Gamma(i / 100.0, j / 11000.0 - 1);
	// cout << "Gamma:\t\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	PowInit();
	c = 0;
	begin = clock();
	for(int i = 1; i < 100; ++i)
		for(int j = 0; j < 110000; ++j)
			c += GetPow(i / 100.0, j / 10000.0 - 1);
	cout << "GetPow:\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// c = 0;
	// begin = clock();
	// for(int i = 1; i < 100000000; ++i)
	// 	c += gsl_sf_gamma(i / 10000000.0);
	// cout << "gsl_sf_gamma: " << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// system("pause");
	return 0;
}