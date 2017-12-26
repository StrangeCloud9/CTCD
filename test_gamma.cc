#include <cmath>
// #include <gsl/gsl_sf_gamma.h>
#include "gamma.cc"
#include <iostream>
#include <ctime>
#include <omp.h>
#include <iomanip>
#include <vector>
using namespace std;


vector<double> gamma_val;
#define GAMMA_LENGTH 100000000
void GammaInit() {
	gamma_val.resize(GAMMA_LENGTH);
	for (int i = 0; i < GAMMA_LENGTH; ++i)
		gamma_val[i] = tgamma(i * 10.0 / GAMMA_LENGTH);
}

double GetGamma(double x) {
	if (x >= 10.0)
		return tgamma(x);
	else
		return gamma_val[(int)(x / 10.0 * GAMMA_LENGTH)];
}


int main(int argc, char const *argv[])
{
	double c = 0;
	int begin = clock();
	for(int i = 1; i < 100000000; ++i)
		c += tgamma(i / 10000000.0);
	cout << "tgamma:\t\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	c = 0;
	begin = clock();
	for(int i = 1; i < 100000000; ++i)
		c += __builtin_tgamma(i / 10000000.0);
	cout << "__builtin:\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	c = 0;
	begin = clock();
	for(int i = 1; i < 100000000; ++i)
		c += Gamma(i / 10000000.0);
	cout << "Gamma:\t\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	GammaInit();
	c = 0;
	begin = clock();
	for(int i = 1; i < 100000000; ++i)
		c += GetGamma(i / 10000000.0);
	cout << "GetGamma:\t" << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// c = 0;
	// begin = clock();
	// for(int i = 1; i < 100000000; ++i)
	// 	c += gsl_sf_gamma(i / 10000000.0);
	// cout << "gsl_sf_gamma: " << setprecision(20) << c << '\t' << ((double)(clock() - begin) / CLOCKS_PER_SEC) << endl;
	// system("pause");
	return 0;
}