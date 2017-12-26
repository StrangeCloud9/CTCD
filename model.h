#ifndef CTCD_MODEL_H
#define CTCD_MODEL_H
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "dataset.h"
#include "utils.h"
#include <vector>
#include <utility>
#include <random>
#include <cmath>
#ifdef __clang__
#include <libiomp/omp.h>
#else
#include <omp.h>
#endif
#include <functional>
#include <algorithm>
using namespace std;

class Model {
  public:
	Model() = default;
	~Model() = default;

	void Init(int argc, char const *argv[]);
	void Run();
	void Save(string tag = "default");
  private:
	void Load(int argc, char const *argv[]);
	void Allocate();
	template<typename T1, typename T2>
	void RandomFillOne(vector< vector<T1> > &arr, T2 dist, default_random_engine &eng) {
		for(int i = 0; i < arr.size(); ++i)
			for (int j = 0; j < arr[i].size(); ++j) {
				arr[i][j] = dist(eng);
			}
	}
	template<typename T>
	void FillWithZero(vector< vector<T> > &arr) {
		for(int i = 0; i < arr.size(); ++i)
			fill(arr[i].begin(), arr[i].end(), 0);
	}
	template<typename T>
	void FillWithZero(vector<T> &arr) {
		fill(arr.begin(), arr.end(), 0);
	}
	void RandomFill();
	void Iterate();
	void Count();
	void Fitting();
	void Sample();

	Dataset dataset_;

	// Setup Arguments
	string output_file_perfix_, input_folder_path_;
	int thread_number_, community_count_, topic_count_, iteration_count_, save_iteration_count_;

	// Hyper Parameters
	double alpha_, beta_, epsilon_, rho_, sigma_;
	pair<double, double> lambda_;

	// To optimized
	vector< vector<int> > z_, c_, s_, s_apos_;
	vector< vector< pair<double, double> > > delta_, gamma_;
	vector< vector< pair<double, double> > > psi_;
	vector< vector<double> > beta_delta_, beta_gamma_;
	vector< vector<double> > beta_psi_;

	// Can be calculated finally
	vector< vector<double> > theta_, pi_, phi_, omega_;
	vector< vector< pair<double, double> > > eta_;

	// Counter
	// update c_i_j
	vector< vector<int> > n_i_cp_, n_cp_k_;
	vector<int> n_i_all_cp_, n_cp_all_k_;
	//update z_i_j
	vector< vector<int> > n_k_v_;
	vector<int> n_k_all_v_;
	// update s_i_i_apos and s_apos_i_i_apos
	vector< vector<int> > &n_i_ci_ = n_i_cp_;
	vector< vector<int> > n_ci_cp_, n_t_ci_cp_;
	vector<int> &n_i_all_ci_ = n_i_all_cp_;
	vector< vector< vector<double> > > poisson_;
	// update psi_cp_k
	vector< vector<double> > n_t_cp_k_, mean_t_cp_k_, var_t_cp_k_;
	// update delta_i_ci
	vector< vector<double> > n_t_i_ci_, mean_t_i_ci_, var_t_i_ci_;
	// update gamma_i_apos_cp
	vector< vector<double> > n_t_i_apos_cp_, mean_t_i_apos_cp_, var_t_i_apos_cp_;
};
#endif
