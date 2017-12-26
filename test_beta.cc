// #define SELF_GAMMA
// #ifdef SELF_GAMMA
// #include "gamma.cc"
// #else
// #define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
// #endif
// #define OMP_ON
#include <cmath>
double beta(double a, double b) {
	return tgamma(a) * tgamma(b) / tgamma(a + b);
}
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <ctime>
#include <omp.h>
using namespace std;


int main(int argc, char const *argv[])
{
	double c = 0;
	int begin = clock();
#ifdef OMP_ON
	#pragma omp parallel for reduction(+:c)
#endif
	for(int i = 1; i < 10000; ++i)
		for(int j = 1; j < 10000; ++j){
			double a = i / 1000.0;
			double b = j / 1000.0;
			c += gsl_sf_beta(a, b);
		}
	cout << (1 / beta(0, 0)) << endl;
	cout << c << endl;
	cout << ((double)(clock() - begin) / CLOCKS_PER_SEC);
	// system("pause");
	return 0;
}