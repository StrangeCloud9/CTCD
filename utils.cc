#include "utils.h"

string Int2String(int value, int min_len) {
	if (value < 0)
		throw invalid_argument("Negative number is provided to converted to string");
	string result = to_string(value);
	if (result.length() < min_len)
		result = string("00000000000000").substr(0, min_len - result.length()) + result;
	return result;
}

long double RealRand(std::mt19937 &generator, const long double min, const long double max) {
	// static thread_local std::mt19937 generator;
	std::uniform_real_distribution<long double> distribution(min, max);
	return distribution(generator);
}


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


vector< vector<double> > pow_val;
#define POW_BASE_LENGTH 100
#define POW_EXP_LENGTH 110000
void PowInit() {
	pow_val.resize(POW_BASE_LENGTH, vector<double>(POW_EXP_LENGTH));
	for (int i = 0; i < POW_BASE_LENGTH; ++i)
		for (int j = 0; j < POW_EXP_LENGTH; ++j)
			pow_val[i][j] = pow((double)i / POW_BASE_LENGTH, (double)j / POW_EXP_LENGTH * 10 - 1);
}

double GetPow(double base, double exponent) {
	if (exponent >= 10)
		return pow(base, exponent);
	else
		return pow_val[(int)(base * POW_BASE_LENGTH)][(int)((exponent + 1) / 11 * POW_BASE_LENGTH)];
}