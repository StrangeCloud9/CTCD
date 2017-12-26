#pragma once
#ifndef CTCD_UTILS_H
#define CTCD_UTILS_H
#include <typeinfo>
#include <string>
#include <stdexcept>
#include <random>
#include <ctime>
#include <iostream>
#include <cmath>
using namespace std;

inline string TimeString() {
	time_t timep;
	time(&timep);
#ifdef _MSC_VER
	char buffer[100];
	ctime_s(buffer, 100, &timep);
	buffer[strlen(buffer) - 1] = '\0';
	return buffer;
#else
	string s = ctime(&timep);
	return s.substr(0, s.length() - 1);
#endif
}


string Int2String(int value, int min_len = 4);



template<typename T>
struct Trans {
	T operator()(string tag, string s) {
		throw invalid_argument("Unsupported data type provided when parse tag " + tag);
	}
};

template<>
struct Trans<int> {
	int operator()(string tag, string s) {
		return stoi(s);
	}
};

template<>
struct Trans<double> {
	double operator()(string tag, string s) {
		return stod(s);
	}
};

template<>
struct Trans<string> {
	string operator()(string tag, string s) {
		return string(s);
	}
};

template<>
struct Trans<bool> {
	bool operator()(string tag, string s) {
		return true;
	}
};


template<typename T>
T GetArgByTag(int argc, char const *argv[], string tag, T default_value, string introduction) {
	// tag should contain '-'
	int tag_pos = 0;
	for (tag_pos = 0; tag_pos < argc && tag != string(argv[tag_pos]); ++tag_pos);
	T return_value;
	if (tag_pos == argc)
		return_value = default_value;
	else {
		if (tag_pos == argc - 1) {
			if (typeid(T) != typeid(bool))
				throw invalid_argument("There is no value for tag " + tag);
			else
				return_value = Trans<T>()(tag, "");
		} else {
			try {
				return_value = Trans<T>()(tag, argv[tag_pos + 1]);
			} catch (const invalid_argument& e) {
				cerr << TimeString() << '\t' << e.what() << endl;
				throw invalid_argument("Some errors occur when parse tag " + tag);
			}
		}
		// if (typeid(default_value).name() == string("b")) // bool
		// 	return_value = dynamic_cast<T>(true);
		// else {
		// 	if (typeid(default_value).name() == string("i")) // int
		// 		return_value = stoi(argv[tag_pos + 1]);
		// 	else if (typeid(default_value).name() == string("d")) // int
		// 		return_value = stod(argv[tag_pos + 1]);
		// 	else if (typeid(default_value).name() == string("s")) // int
		// 		return_value = string(argv[tag_pos + 1]);
		// 	else
		// 		throw invalid_argument("Unsupported data type provided when parse tag " + tag);
		// }
	}
	cout << introduction << "(default: " << default_value << ")\t=" << return_value << endl;
	return return_value;
}


template<typename T>
inline T Square(T val) {
	return val * val;
}

long double RealRand(std::mt19937 &generator, const long double min, const long double max);

void GammaInit();

double GetGamma(double x);

void PowInit();

double GetPow(double base, double exponent);

template<typename T>
ostream &operator<<(ostream &os, pair<T, T> x){
	os << x.first << ',' << x.second;
	return os;
}

template<typename T>
bool AlmostZero(T x) {
	return abs(x) < 1e-5;
}
#endif