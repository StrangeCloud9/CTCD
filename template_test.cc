#include <iostream>
#include <string>
#include <typeinfo>
using namespace std;

template<typename T>
struct Trans{
	T operator()(string s){
		throw invalid_argument("Unsupported data type provided when parse tag " + tag);
	}
};

template<>
struct Trans<int>{
	int operator()(string s){
		return stoi(s);
	}
};

template<>
struct Trans<double>{
	double operator()(string s){
		return stod(s);
	}
};

template<>
struct Trans<string>{
	string operator()(string s){
		return string(s);
	}
};

template<>
struct Trans<bool>{
	bool operator()(string s){
		return true;
	}
};

int main(int argc, char const *argv[])
{
	int a = 1;
	cout << Trans<double>()("1.2");
	return 0;
}