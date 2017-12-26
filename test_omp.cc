#include <omp.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
using namespace std;
#define SIZE 100000000


double intRand(std::mt19937 &generator, const int & min, const int & max) {
	// static thread_local std::mt19937 generator;
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}


class test
{
public:
	test(){
		cout << "construct";
	}
	~test(){
		cout << "destruct";
	}
	
};


int main(int argc, char const *argv[]) {
	// int *a = new int[SIZE], *b = new int[SIZE];
	// int begin = clock();
	// int sum;
	// // #pragma omp parallel for num_threads(4) private(sum) reduce(sum:+) //schedule(static, 0x1000)
	// for(int i = 0; i < SIZE; ++i)
	// 	for(int j = 0; j < SIZE; ++j)
	// 		sum += a[i] * b[i];
	// cout << clock() - begin << endl;
	// cout << sum << endl;
	cout << omp_get_thread_num();
	#pragma omp parallel
	{
		thread_local std::mt19937 generator;
		generator.seed((unsigned)(time(0) + 10000 * omp_get_thread_num()));
		test tmp = test();
		#pragma omp for
		for (int i = 0; i < 8; ++i) {
			cout << intRand(generator, 0, 10) << '\n';
		}

	}
	// #pragma omp parallel for num_threads(4)
	return 0;
}
