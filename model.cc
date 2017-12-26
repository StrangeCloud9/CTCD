#define OMP_ON
// #define DEBUG
// #define POISSON
#include "model.h"
// #include "gamma.h"
inline double beta(pair<double, double> x) {
#ifdef _MSC_VER
	double result = tgamma(x.first) * tgamma(x.second) / tgamma(x.first + x.second);
#else
	double result = beta(x.first, x.second);
#endif
	if (!isfinite(result))
		return 0;
	return result;
}
void Model::Load(int argc, char const *argv[]) {
	input_folder_path_ = GetArgByTag(argc, argv, "-i", string("A:\\Sam\\Lab\\SIGIR2017\\Data\\raw\\dblp_citation\\"), "Input data's folder path");
	dataset_.set_in_dir(input_folder_path_);
	dataset_.LoadData();
	dataset_.NormalizeTime();
	output_file_perfix_ = GetArgByTag(argc, argv, "-o", string("A:\\Sam\\Lab\\SIGIR2017\\Result\\CTCD\\test"), "Output Graph data prefix");
	thread_number_ = GetArgByTag(argc, argv, "-tm", omp_get_num_procs(), "Number of threads for parallelization");
	omp_set_num_threads(thread_number_);
	community_count_ = GetArgByTag(argc, argv, "-cc", 100, "The number of communities to detect");
	topic_count_ = GetArgByTag(argc, argv, "-tc", 100, "The number of topics in one community to detect");
	iteration_count_ = GetArgByTag(argc, argv, "-ic", 5000, "Number of update iterations");
	save_iteration_count_ = GetArgByTag(argc, argv, "-sic", 100, "How many iterations for once save");
	alpha_ = GetArgByTag(argc, argv, "-alpha", 0.1, "Alpha");
	beta_ = GetArgByTag(argc, argv, "-beta", 0.1, "Beta");
	epsilon_ = GetArgByTag(argc, argv, "-epsilon", 0.1, "Epsilon");
	rho_ = GetArgByTag(argc, argv, "-rho", 0.1, "Rho");
	sigma_ = GetArgByTag(argc, argv, "-sigma", 0.1, "Sigma");
	double lambda_weight = GetArgByTag(argc, argv, "-lambdaweight", 1.0, "Lambda's tunable weight");
	lambda_.first = GetArgByTag(argc, argv, "-lambda0", 0.1, "Lambda0");
	double default_lambda1 = lambda_weight * log(
	                             ((double)dataset_.user_count() * (dataset_.user_count() - 1) - dataset_.edge_count())
	                             / (community_count_ * community_count_));
	// double default_lambda1 = dataset_.edge_count() / (community_count_ * community_count_);
	lambda_.second = GetArgByTag(argc, argv, "-lambda1", default_lambda1, "Lambda1");
}


void Model::Allocate() {
	z_.clear();
	z_.reserve(dataset_.user_count());
	c_.clear();
	c_.reserve(dataset_.user_count());
	for (int i = 0; i < dataset_.user_count(); ++i) {
		z_.push_back(vector<int>(dataset_.users_[i].docs.size()));
		c_.push_back(vector<int>(dataset_.users_[i].docs.size()));
	}
	s_.clear();
	s_.reserve(dataset_.user_count());
	s_apos_.clear();
	s_apos_.reserve(dataset_.user_count());
	for (int i = 0; i < dataset_.user_count(); ++i) {
		s_.push_back(vector<int>(dataset_.users_[i].out_edges.size()));
		s_apos_.push_back(vector<int>(dataset_.users_[i].out_edges.size()));
	}
	theta_.assign(community_count_, vector<double>(topic_count_));
	pi_.assign(dataset_.user_count(), vector<double>(community_count_));
	omega_.assign(dataset_.user_count(), vector<double>(community_count_));
	delta_.assign(dataset_.user_count(), vector< pair<double, double> >(community_count_));
	beta_delta_.assign(dataset_.user_count(), vector<double>(community_count_));
	gamma_.assign(dataset_.user_count(), vector< pair<double, double> >(community_count_));
	beta_gamma_.assign(dataset_.user_count(), vector<double>(community_count_));
	eta_.assign(community_count_, vector< pair<double, double> >(community_count_));
	phi_.assign(topic_count_, vector<double>(dataset_.word_count()));
	psi_.assign(community_count_, vector< pair<double, double> >(topic_count_));
	beta_psi_.assign(community_count_, vector<double>(topic_count_));
	//counters
	n_i_cp_.assign(dataset_.user_count(), vector<int>(community_count_));
	n_cp_k_.assign(community_count_, vector<int>(topic_count_));
	n_i_all_cp_.assign(dataset_.user_count(), 0); // used when save
	n_cp_all_k_.assign(community_count_, 0);
	n_k_v_.assign(topic_count_, vector<int>(dataset_.word_count()));
	n_k_all_v_.assign(topic_count_, 0);
	n_i_ci_.assign(dataset_.user_count(), vector<int>(community_count_));
	n_i_all_ci_.assign(dataset_.user_count(), 0); // used when save
	n_ci_cp_.assign(community_count_, vector<int>(community_count_));
	n_t_ci_cp_.assign(community_count_, vector<int>(community_count_));

	n_t_cp_k_.assign(community_count_, vector<double>(topic_count_));
	mean_t_cp_k_.assign(community_count_, vector<double>(topic_count_));
	var_t_cp_k_.assign(community_count_, vector<double>(topic_count_));
	n_t_i_ci_.assign(dataset_.user_count(), vector<double>(community_count_));
	mean_t_i_ci_.assign(dataset_.user_count(), vector<double>(community_count_));
	var_t_i_ci_.assign(dataset_.user_count(), vector<double>(community_count_));
	n_t_i_apos_cp_.assign(dataset_.user_count(), vector<double>(community_count_));
	mean_t_i_apos_cp_.assign(dataset_.user_count(), vector<double>(community_count_));
	var_t_i_apos_cp_.assign(dataset_.user_count(), vector<double>(community_count_));

	poisson_.assign(dataset_.max_pair_count() + 1, vector< vector<double> >(community_count_, vector<double>(community_count_)));
}


void Model::RandomFill() {
	random_device r;
	default_random_engine eng(r());
	RandomFillOne(z_, uniform_int_distribution<int>(0, topic_count_ - 1), eng);
	RandomFillOne(c_, uniform_int_distribution<int>(0, community_count_ - 1), eng);
	RandomFillOne(s_, uniform_int_distribution<int>(0, community_count_ - 1), eng);
	RandomFillOne(s_apos_, uniform_int_distribution<int>(0, community_count_ - 1), eng);
}


void Model::Init(int argc, char const *argv[]) {
	Load(argc, argv);
	cout << TimeString() << " Model loads done" << endl;
	Allocate();
	cout << TimeString() << " Model allocates done" << endl;
	RandomFill();
	cout << TimeString() << " Model initiates done" << endl;
	// GammaInit();
	// cout << TimeString() << " Gamma initiates done" << endl;
	// PowInit();
	// cout << TimeString() << " Pow initiates done" << endl;
}


void Model::Count() {
	FillWithZero(n_i_cp_);
	FillWithZero(n_cp_k_);
	FillWithZero(n_i_all_cp_); // used when save
	FillWithZero(n_cp_all_k_);
	FillWithZero(n_k_v_);
	FillWithZero(n_k_all_v_);
	FillWithZero(n_i_ci_);
	FillWithZero(n_i_all_ci_); // used when save
	FillWithZero(n_ci_cp_);
	FillWithZero(n_t_ci_cp_);

	FillWithZero(n_t_cp_k_);
	FillWithZero(mean_t_cp_k_);
	FillWithZero(var_t_cp_k_);
	FillWithZero(n_t_i_ci_);
	FillWithZero(mean_t_i_ci_);
	FillWithZero(var_t_i_ci_);
	FillWithZero(n_t_i_apos_cp_);
	FillWithZero(mean_t_i_apos_cp_);
	FillWithZero(var_t_i_apos_cp_);
	int cp, k, ci, v, i_apos;
	for (int i = 0; i < dataset_.users_.size(); ++i) {
		for (int j = 0; j < dataset_.users_[i].docs.size(); ++j) {
			cp = c_[i][j];
			k = z_[i][j];
			++n_i_cp_[i][cp];
			++n_cp_k_[cp][k];
			++n_i_all_cp_[i]; // used when save
			++n_cp_all_k_[cp];
			for (int w = 0; w < dataset_.users_[i].docs[j].words.size(); ++w) {
				++n_k_v_[k][dataset_.users_[i].docs[j].words[w]];
			}
			n_k_all_v_[k] += dataset_.users_[i].docs[j].words.size();
		}
		for (int out_i = 0; out_i < dataset_.users_[i].out_edges.size(); ++out_i) {
			ci = s_[i][out_i];
			cp = s_apos_[i][out_i];
			++n_i_ci_[i][ci];
			i_apos = dataset_.users_[i].out_edges[out_i].first;
			++n_i_cp_[i_apos][cp];
			++n_i_all_ci_[i]; // used when save
			++n_i_all_cp_[i_apos]; // used when save
			++n_ci_cp_[ci][cp];
			n_t_ci_cp_[ci][cp] += dataset_.users_[i].out_edges[out_i].second.times.size();
		}
	}

	for (int i = 0; i < dataset_.users_.size(); ++i) {
		for (int j = 0; j < dataset_.users_[i].docs.size(); ++j) {
			cp = c_[i][j];
			k = z_[i][j];
			++n_t_cp_k_[cp][k];
			mean_t_cp_k_[cp][k] += dataset_.users_[i].docs[j].time;
			var_t_cp_k_[cp][k] += Square(dataset_.users_[i].docs[j].time);
		}
		for (int out_i = 0; out_i < dataset_.users_[i].out_edges.size(); ++out_i) {
			ci = s_[i][out_i];
			cp = s_apos_[i][out_i];
			int i_apos = dataset_.users_[i].out_edges[out_i].first;
			n_t_i_ci_[i][ci] += dataset_.users_[i].out_edges[out_i].second.size();
			n_t_i_apos_cp_[i_apos][cp] += dataset_.users_[i].out_edges[out_i].second.size();
			for (int m = 0; m < dataset_.users_[i].out_edges[out_i].second.times.size(); ++m) {
				mean_t_i_ci_[i][ci] += dataset_.users_[i].out_edges[out_i].second.times[m].first;
				var_t_i_ci_[i][ci] += Square(dataset_.users_[i].out_edges[out_i].second.times[m].first);
				mean_t_i_apos_cp_[i_apos][cp] += dataset_.users_[i].out_edges[out_i].second.times[m].second;
				var_t_i_apos_cp_[i_apos][cp] += Square(dataset_.users_[i].out_edges[out_i].second.times[m].second);
			}
		}
	}
	for (int c = 0; c < community_count_; ++c)
		for (int k = 0; k < topic_count_; ++k)
			if (n_t_cp_k_[c][k] > 0) {
				mean_t_cp_k_[c][k] /= n_t_cp_k_[c][k];
				var_t_cp_k_[c][k] = var_t_cp_k_[c][k] / n_t_cp_k_[c][k] - Square(mean_t_cp_k_[c][k]);
			}
	for (int i = 0; i < dataset_.users_.size(); ++i)
		for (int c = 0; c < community_count_; ++c) {
			if (n_t_i_ci_[i][c] > 0) {
				mean_t_i_ci_[i][c] /= n_t_i_ci_[i][c];
				var_t_i_ci_[i][c] = var_t_i_ci_[i][c] / n_t_i_ci_[i][c] - Square(mean_t_i_ci_[i][c]);
			}
			if (n_t_i_apos_cp_[i][c] > 0) {
				mean_t_i_apos_cp_[i][c] /= n_t_i_apos_cp_[i][c];
				var_t_i_apos_cp_[i][c] = var_t_i_apos_cp_[i][c] / n_t_i_apos_cp_[i][c] - Square(mean_t_i_apos_cp_[i][c]);
			}
		}

	// for (int i = 0; i < dataset_.users_.size(); ++i) {
	// 	for (int j = 0; j < dataset_.users_[i].docs.size(); ++j) {
	// 		cp = c_[i][j];
	// 		k = z_[i][j];
	// 		var_t_cp_k_[cp][k] += Square(dataset_.users_[i].docs[j].time - mean_t_cp_k_[cp][k]);
	// 	}
	// 	for (int out_i = 0; out_i < dataset_.users_[i].out_edges.size(); ++out_i) {
	// 		ci = s_[i][out_i];
	// 		cp = s_apos_[i][out_i];
	// 		int i_apos = dataset_.users_[i].out_edges[out_i].first;
	// 		for (int m = 0; m < dataset_.users_[i].out_edges[out_i].second.times.size(); ++m) {
	// 			var_t_i_ci_[i][ci] += Square(dataset_.users_[i].out_edges[out_i].second.times[m].first - mean_t_i_ci_[i][ci]);
	// 			var_t_i_apos_cp_[i_apos][cp] += Square(dataset_.users_[i].out_edges[out_i].second.times[m].second - mean_t_i_apos_cp_[i_apos][cp]);
	// 		}
	// 	}
	// }
	// for (int c = 0; c < community_count_; ++c)
	// 	for (int k = 0; k < topic_count_; ++k)
	// 		if (n_t_cp_k_[c][k] > 1)
	// 			var_t_cp_k_[c][k] /= n_t_cp_k_[c][k];
	// for (int i = 0; i < dataset_.users_.size(); ++i)
	// 	for (int c = 0; c < community_count_; ++c) {
	// 		if (n_t_i_ci_[i][c] > 1)
	// 			var_t_i_ci_[i][c] /= n_t_i_ci_[i][c];
	// 		if (n_t_i_apos_cp_[i][c] > 1)
	// 			var_t_i_apos_cp_[i][c] /= n_t_i_apos_cp_[i][c];
	// 	}
}


void Model::Fitting() {
#ifdef OMP_ON
	#pragma omp parallel
#endif
	{
		// psi_cp_k
		double mean, var;
		// int not_zero_count = 0;
#ifdef OMP_ON
		#pragma omp for schedule(dynamic, 32)
#endif
		for (int c = 0; c < community_count_; ++c)
			for (int k = 0; k < topic_count_; ++k) {
				mean = mean_t_cp_k_[c][k];
				var = var_t_cp_k_[c][k];
				if (!AlmostZero(var)) {
					// ++not_zero_count;
					psi_[c][k].first = mean * (mean * (1 - mean) / var - 1);
					psi_[c][k].second = (1 - mean) * (mean * (1 - mean) / var - 1);
				} else {
					psi_[c][k] = pair<double, double>(1, 1);
				}
				beta_psi_[c][k] = beta(psi_[c][k]);
			}
		// delta_i_ci, gamma_i_apos_cp
#ifdef OMP_ON
		#pragma omp for schedule(dynamic, 32)
#endif
		for (int i = 0; i < dataset_.user_count(); ++i) {
			for (int c = 0; c < community_count_; ++c) {
				mean = mean_t_i_ci_[i][c];
				var = var_t_i_ci_[i][c];
				if (!AlmostZero(var)) {
					// ++not_zero_count;
					delta_[i][c].first = mean * (mean * (1 - mean) / var - 1);
					delta_[i][c].second = (1 - mean) * (mean * (1 - mean) / var - 1);
				} else {
					delta_[i][c] = pair<double, double>(1, 1);
				}
				beta_delta_[i][c] = beta(delta_[i][c]);
				mean = mean_t_i_apos_cp_[i][c];
				var = var_t_i_apos_cp_[i][c];
				if (!AlmostZero(var)) {
					// ++not_zero_count;
					gamma_[i][c].first = mean * (mean * (1 - mean) / var - 1);
					gamma_[i][c].second = (1 - mean) * (mean * (1 - mean) / var - 1);

				} else {
					gamma_[i][c] = pair<double, double>(1, 1);
				}
				beta_gamma_[i][c] = beta(gamma_[i][c]);
			}
		}
	}
}


void Model::Sample() {
#ifdef OMP_ON
	#pragma omp parallel
#endif
	{
		std::mt19937 generator(time(0) + 10000 * omp_get_thread_num());
		// cout << "not zero count: " << not_zero_count << endl;
		// c_i_j, z_i_j
		vector<double> pc(community_count_), pk(topic_count_);
		int k, c, new_c, v, new_z;
		double t;
		long double tmp, tmp2;
#ifdef OMP_ON
		#pragma omp for schedule(dynamic, 32)
#endif
		for (int i = 0; i < dataset_.users_.size(); ++i) {
			for (int j = 0; j < dataset_.users_[i].docs.size(); ++j) {
				c = c_[i][j];
				k = z_[i][j];
				t = dataset_.users_[i].docs[j].time;

				for (int c = 0; c < community_count_; ++c) {
					pc[c] = (n_i_cp_[i][c] + rho_) * (n_cp_k_[c][k] + alpha_) * (pow(t, psi_[c][k].first - 1) * pow(1 - t, psi_[c][k].second - 1))
					/ ((n_cp_all_k_[c] + topic_count_ * alpha_) * beta_psi_[c][k]);
					// if (isinf(pc[c])) {
					// 	cout << (n_i_cp_[i][c] + rho_) << " " << (n_cp_k_[c][k] + alpha_) << " " << pow(1 - t, psi_[c][k].first - 1) << " " << t << " " << psi_[c][k].second - 1
					// 		 << " " << (n_cp_all_k_[c] + topic_count_ * alpha_) << " " << beta_psi_[c][k];
					// 	cout << endl;
					// }
				}
				pc[c] = (n_i_cp_[i][c] - 1 + rho_) * (n_cp_k_[c][k] - 1 + alpha_) * (pow(t, psi_[c][k].first - 1) * pow(1 - t, psi_[c][k].second - 1))
				        / ((n_cp_all_k_[c] - 1 + topic_count_ * alpha_) * beta_psi_[c][k]);
				for (int c = 1; c < community_count_; ++c)
					pc[c] += pc[c - 1];
				tmp = RealRand(generator, 0, pc[community_count_ - 1]);
				for (new_c = 0; pc[new_c] <= tmp && new_c < community_count_ - 1; ++new_c);
				// if (pc[new_c] <= tmp) {
				// 	cout << pc[new_c] << ' ' << tmp << '\n';
				// }
				c_[i][j] = new_c;

				for (int k = 0; k < topic_count_; ++k) {
					tmp = tmp2 = 1;
					for (int w = 0; w < dataset_.users_[i].docs[j].words.size(); ++w) {
						v = dataset_.users_[i].docs[j].words[w];
						tmp *= n_k_v_[k][v] + dataset_.users_[i].docs[j].counter[w] + beta_;
						tmp2 *= n_k_all_v_[k] + w + dataset_.word_count() * beta_;
					}
					pk[k] = (n_cp_k_[c][k] + alpha_) * (pow(t, psi_[c][k].first - 1) * pow(1 - t, psi_[c][k].second - 1)) * tmp
					        / (beta_psi_[c][k] * tmp2);
				}
				tmp = tmp2 = 1;
				for (int w = 0; w < dataset_.users_[i].docs[j].words.size(); ++w) {
					v = dataset_.users_[i].docs[j].words[w];
					tmp *= n_k_v_[k][v] - dataset_.users_[i].docs[j].counter[w] - 1 + beta_;
					tmp2 *= n_k_all_v_[k] - w - 1 + dataset_.word_count() * beta_;
				}
				pk[k] = (n_cp_k_[c][k] - 1 + alpha_) * (pow(t, psi_[c][k].first - 1) * pow(1 - t, psi_[c][k].second - 1)) * tmp
				        / (beta_psi_[c][k] * tmp2);
				for (int k = 1; k < topic_count_; ++k)
					pk[k] += pk[k - 1];
				tmp = RealRand(generator, 0, pk[topic_count_ - 1]);
				for (new_z = 0; pk[new_z] <= tmp && new_z < topic_count_ - 1; ++new_z);
				z_[i][j] = new_z;
			}
		}
		// s_i_i_apos, s_apos_i_i_apos
		int i_apos, ci, cp, n_i_i_apos;
		long double tmp0;
		vector<long double> pcc(community_count_ * community_count_);
		vector< vector< pair<double, double> > > pow_t_delta, pow_t_gamma;
		pow_t_delta.assign(community_count_, vector< pair<double, double> >(dataset_.max_pair_count()));
		pow_t_gamma.assign(community_count_, vector< pair<double, double> >(dataset_.max_pair_count()));
#ifdef POISSON
#ifdef OMP_ON
		#pragma omp for schedule(dynamic, 32)
#endif
		// PreCalculation for Poison Distribution
		for (int n_i_i_apos = 0; n_i_i_apos <= dataset_.max_pair_count(); ++n_i_i_apos) {
			if (dataset_.has_pair_count(n_i_i_apos))
				for (int ci = 0; ci < community_count_; ++ci)
					for (int cp = 0; cp < community_count_; ++cp) {
						tmp0 = 1;
						for (int m = 0; m < n_i_i_apos; ++m)
							tmp0 *= lambda_.first + n_t_ci_cp_[ci][cp] + m;
						poisson_[n_i_i_apos][ci][cp] = pow((long double)(lambda_.second + n_ci_cp_[ci][cp]) / (lambda_.second + n_ci_cp_[ci][cp] + 1), (long double)lambda_.first + n_t_ci_cp_[ci][cp]) * tmp0
						                               / (long double)pow(lambda_.second + n_ci_cp_[ci][cp] + 1, n_i_i_apos);
					}
		}
#endif
#ifdef OMP_ON
		#pragma omp for schedule(dynamic, 32)
#endif
		for (int i = 0; i < dataset_.users_.size(); ++i) {
			for (int out_i = 0; out_i < dataset_.users_[i].out_edges.size(); ++out_i) {
				// ++edge_count;
				// if (edge_count % 10000 == 0) {
				// 	cout << edge_count << endl;
				// }
				ci = s_[i][out_i];
				cp = s_apos_[i][out_i];
				i_apos = dataset_.users_[i].out_edges[out_i].first;
				n_i_i_apos = dataset_.users_[i].out_edges[out_i].second.times.size();
				for (int c = 0; c < community_count_; ++c) {
					if (!AlmostZero(beta_delta_[i][c]))
						for (int m = 0; m < n_i_i_apos; ++m) {
							t = dataset_.users_[i].out_edges[out_i].second.times[m].first;
							// cout << 't' << t << '\n';
							// pow_t_delta[c][m].first = pow(t, delta_[i][c].first - 1);
							// pow_t_delta[c][m].second = pow(1 - t, delta_[i][c].second - 1);
							pow_t_delta[c][m] = pair<double, double>(pow(t, delta_[i][c].first - 1), pow(1 - t, delta_[i][c].second - 1));
							// cout << pair<double, double>(pow(t, delta_[i][c].first - 1), pow(1 - t, delta_[i][c].second - 1)) << '\n';
						}
					if (!AlmostZero(beta_gamma_[i_apos][c]))
						for (int m = 0; m < n_i_i_apos; ++m) {
							t = dataset_.users_[i].out_edges[out_i].second.times[m].second;
							// pow_t_gamma[c][m].first = pow(t, gamma_[i_apos][c].first - 1);
							// pow_t_gamma[c][m].second = pow(1 - t, gamma_[i_apos][c].second - 1);
							pow_t_gamma[c][m] = pair<double, double>(pow(t, gamma_[i_apos][c].first - 1), pow(1 - t, gamma_[i_apos][c].second - 1));
							// cout << pair<double, double>(pow(t, gamma_[i_apos][c].first - 1), pow(1 - t, gamma_[i_apos][c].second - 1)) << '\n';
						}
				}
				for (int ci = 0; ci < community_count_; ++ci) {
					for (int cp = 0; cp < community_count_; ++cp) {
						tmp0 = tmp = tmp2 = 1;
						// for (int m = 0; m < n_i_i_apos; ++m) {
						// 	tmp0 *= lambda_.first + n_t_ci_cp_[ci][cp] + m;
						// }
						if (!AlmostZero(beta_delta_[i][ci]))
							for (int m = 0; m < n_i_i_apos; ++m) {
								t = dataset_.users_[i].out_edges[out_i].second.times[m].first;
								tmp *= pow_t_delta[ci][m].first * pow_t_delta[ci][m].second / beta_delta_[i][ci];
								// tmp *= pow(t, delta_[i][ci].first - 1) * pow(1 - t, delta_[i][ci].second - 1);
								// tmp2 *= beta_delta_[i][ci];
							}
						if (!AlmostZero(beta_gamma_[i_apos][cp]))
							for (int m = 0; m < n_i_i_apos; ++m) {
								t = dataset_.users_[i].out_edges[out_i].second.times[m].second;
								tmp *= pow_t_gamma[cp][m].first * pow_t_gamma[cp][m].second / beta_gamma_[i_apos][cp];
								// tmp *= pow(t, gamma_[i_apos][cp].first - 1) * pow(1 - t, gamma_[i_apos][cp].second - 1);
								// tmp2 *= beta_gamma_[i_apos][cp];
							}
						// cout << tmp << ' ' << tmp2 << '\n';
#ifdef POISSON
						pcc[ci * community_count_ + cp] = (n_i_ci_[i][ci] + sigma_) * (n_i_cp_[i][cp] + rho_) * poisson_[n_i_i_apos][ci][cp] * tmp
						                                  / (tmp2);
#else
						pcc[ci * community_count_ + cp] = (n_i_ci_[i][ci] + sigma_) * (n_i_cp_[i][cp] + rho_) * (lambda_.first + n_ci_cp_[ci][cp]) * tmp
						                                  / ((lambda_.first + lambda_.second + n_ci_cp_[ci][cp]) * tmp2);
#ifdef DEBUG
						if (!isfinite(pcc[ci * community_count_ + cp])) {
							cout << i << '\t' << i_apos << '\t' << ci << '\t' << cp << endl;
							cout << tmp << '\t' << tmp2 << endl;
							cout << pcc[ci * community_count_ + cp] << '\t' << beta_delta_[i][ci] << '\t' << beta_gamma_[i_apos][cp];
							cout << endl;
							exit(0);
						}
#endif
#endif
					}
				}
				tmp0 = tmp = tmp2 = 1;
				for (int m = 0; m < n_i_i_apos; ++m) {
					tmp0 *= lambda_.first + n_t_ci_cp_[ci][cp] - n_i_i_apos + m;
				}
				if (!AlmostZero(beta_delta_[i][ci]))
					for (int m = 0; m < n_i_i_apos; ++m) {
						t = dataset_.users_[i].out_edges[out_i].second.times[m].first;
						tmp *= pow(t, delta_[i][ci].first - 1) * pow(1 - t, delta_[i][ci].second - 1) / beta_delta_[i][ci];
						// tmp2 *= beta_delta_[i][ci];
					}
				if (!AlmostZero(beta_gamma_[i_apos][cp]))
					for (int m = 0; m < n_i_i_apos; ++m) {
						t = dataset_.users_[i].out_edges[out_i].second.times[m].second;
						tmp *= pow(t, gamma_[i_apos][cp].first - 1) * pow(1 - t, gamma_[i_apos][cp].second - 1) / beta_gamma_[i_apos][cp];
						// tmp2 *= beta_gamma_[i_apos][cp];
					}
				// cout << tmp << ' ' << tmp2 << '\n';
#ifdef POISSON
				pcc[ci * community_count_ + cp] = (n_i_ci_[i][ci] - 1 + sigma_) * (n_i_cp_[i][cp] - 1 + rho_) * (pow((lambda_.second + n_ci_cp_[ci][cp] - 1) / (lambda_.second + n_ci_cp_[ci][cp]), lambda_.first + n_t_ci_cp_[ci][cp] - n_i_i_apos) * tmp0) * tmp
				                                  / (pow(lambda_.second + n_ci_cp_[ci][cp], n_i_i_apos) * tmp2);
#else
				pcc[ci * community_count_ + cp] = (n_i_ci_[i][ci] - 1 + sigma_) * (n_i_cp_[i][cp] - 1 + rho_) * (lambda_.first + n_ci_cp_[ci][cp] - 1) * tmp
				                                  / ((lambda_.first + lambda_.second + n_ci_cp_[ci][cp] - 1) * tmp2);
#endif
#ifdef DEBUG
				if (!isfinite(pcc[ci * community_count_ + cp])) {
					cout << i << '\t' << i_apos << '\t' << c << endl;
					cout << pcc[ci * community_count_ + cp] << '\t' << beta_delta_[i][ci] << '\t' << beta_gamma_[i_apos][cp];
					cout << endl;
					exit(0);
				}
				// for (auto t : pcc)
				// 	cout << t << ' ';
				for (int c = 0; c < community_count_ * community_count_; ++c) {
					if (!isfinite(pcc[c])) {
						cout << i << '\t' << i_apos << '\t' << c << endl;
						cout << pcc[c] << '\t' << beta_delta_[i][c / community_count_] << '\t' << beta_gamma_[i_apos][c % community_count_];
						cout << endl;
						exit(0);
					}
				}
#endif
				for (int c = 1; c < community_count_ * community_count_; ++c)
					pcc[c] += pcc[c - 1];
				tmp = RealRand(generator, 0, pcc.back());
				// cout << pcc.back() << endl;
				for (new_c = 0; pcc[new_c] <= tmp && new_c < pcc.size() - 1; ++new_c);
				s_[i][out_i] = new_c / community_count_;
				s_apos_[i][out_i] = new_c % community_count_;
			}
		}
		// cout << "edge_count: " << edge_count << endl;
	}
}


void Model::Iterate() {
	Count();
	Fitting();
	Sample();
}



void Model::Save(string tag) {
	// topic_word
	Count();
	Fitting();
	string filename = output_file_perfix_ + "_" + tag + ".topic_word.txt";
	ofstream fout(filename);
	for (int k = 0; k < topic_count_; ++k) {
		vector< pair<int, int> > freq_words;
		for (int w = 0; w < dataset_.word_count(); ++w)
			if (n_k_v_[k][w] > 0)
				freq_words.push_back(pair<int, int>(n_k_v_[k][w], w));
		partial_sort(freq_words.begin(), freq_words.begin() + min(50, (int)freq_words.size()), freq_words.end(), greater< pair<int, int> >());
		for (int i = 0; i < min(50, (int)freq_words.size()); ++i)
			fout << freq_words[i].first << ',' << dataset_.wid2word(freq_words[i].second) << '\t';
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".pi.txt";
	fout.open(filename);
	for (int i = 0; i < dataset_.user_count(); ++i) {
		for (int cp = 0; cp < community_count_; ++cp) {
			pi_[i][cp] = (n_i_cp_[i][cp] + rho_) / (n_i_all_cp_[i] + community_count_ * rho_);
			fout << pi_[i][cp] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".omega.txt";
	fout.open(filename);
	for (int i = 0; i < dataset_.user_count(); ++i) {
		for (int ci = 0; ci < community_count_; ++ci) {
			omega_[i][ci] = (n_i_ci_[i][ci] + sigma_) / (n_i_all_ci_[i] + community_count_ * sigma_);
			fout << omega_[i][ci] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".theta.txt";
	fout.open(filename);
	for (int cp = 0; cp < community_count_; ++cp) {
		for (int k = 0; k < topic_count_; ++k) {
			theta_[cp][k] = (n_cp_k_[cp][k] + sigma_) / (n_cp_all_k_[cp] + topic_count_ * sigma_);
			fout << theta_[cp][k] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".phi.txt";
	fout.open(filename);
	for (int k = 0; k < topic_count_; ++k) {
		for (int v = 0; v < dataset_.word_count(); ++v) {
			phi_[k][v] = (n_k_v_[k][v] + beta_) / (n_k_all_v_[k] + dataset_.word_count() * beta_);
			fout << phi_[k][v] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".eta.txt";
	fout.open(filename);
	for (int ci = 0; ci < community_count_; ++ci) {
		for (int cp = 0; cp < community_count_; ++cp) {
#ifdef POISSON
			eta_[ci][cp] = pair<double, double>(lambda_.first + n_t_ci_cp_[ci][cp], lambda_.second + n_ci_cp_[ci][cp]);
			fout << eta_[ci][cp] << '\t';
#else
			fout << ((lambda_.first + n_ci_cp_[ci][cp]) / (lambda_.first + lambda_.second + n_ci_cp_[ci][cp])) << '\t';
#endif
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".ncc.txt";
	fout.open(filename);
	for (int ci = 0; ci < community_count_; ++ci) {
		for (int cp = 0; cp < community_count_; ++cp) {
			fout << n_ci_cp_[ci][cp] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".psi.txt";
	fout.open(filename);
	for (int cp = 0; cp < community_count_; ++cp) {
		for (int k = 0; k < topic_count_; ++k) {
			fout << psi_[cp][k] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".delta.txt";
	fout.open(filename);
	for (int i = 0; i < dataset_.user_count(); ++i) {
		for (int ci = 0; ci < community_count_; ++ci) {
			fout << delta_[i][ci] << '\t';
		}
		fout << '\n';
	}
	fout.close();

	filename = output_file_perfix_ + "_" + tag + ".gamma.txt";
	fout.open(filename);
	for (int i = 0; i < dataset_.user_count(); ++i) {
		for (int cp = 0; cp < community_count_; ++cp) {
			fout << gamma_[i][cp] << '\t';
		}
		fout << '\n';
	}
	fout.close();
}


void Model::Run() {
	for (int i = 0; i < iteration_count_; ++i) {
		cout << TimeString() << " iteration " << i << endl;
		Iterate();
		if (i % save_iteration_count_ == 0) {
			Save(Int2String(i));
			cout << TimeString() << '\t' << "Save done, " << i << " iterations done." << endl;
		}
	}
}