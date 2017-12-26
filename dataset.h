#ifndef CTCD_DATASET_H
#define CTCD_DATASET_H
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
using namespace std;


struct Document {
	inline int size() const {
		return words.size();
	}
	inline void clear() {
		words.clear();
	}
	inline void add_word(int wid) {
		words.push_back(wid);
	}
	void calc_counter() {
		counter.resize(words.size());
		int count;
		for (int i = 0; i < words.size(); ++i) {
			count = 0;
			for (int j = 0; j < i; ++j) {
				if (words[j] == words[i])
					++count;
			}
			counter[i] = count;
		}
	}
	vector<int> words;
	vector<int> counter;
	double time;
};


struct Edge {
	inline int size() const {
		return times.size();
	}
	inline void clear() {
		times.clear();
	}
	inline void add_time_pair(pair<double, double> &&time) {
		times.push_back(time);
	}
	vector< pair<double, double> > times;
};


struct User {
	inline void add_doc(Document &doc) {
		docs.push_back(doc);
	}
	inline void add_edge(int user, Edge &edge) {
		out_edges.push_back(pair<int, Edge>(user, edge));
	}
	vector<Document> docs;
	vector< pair<int, Edge> > out_edges;
};

class Dataset {
  public:
	Dataset() : min_time_(0x7fffffff), max_time_(0) {}
	~Dataset() = default;

	inline string in_dir() {
		return in_dir_;
	}
	inline void set_in_dir(string in_dir) {
		in_dir_ = in_dir;
	}
	inline double normalize_one_time(double time) const {
		return (time - min_time_) / (max_time_ - min_time_);
	}
	inline double unnormalize_one_time(double time) const {
		return time * (max_time_ - min_time_) + min_time_;
	}
	inline int user_count() const {
		return user2uid_.size();
	}
	inline int word_count() const {
		return word2wid_.size();
	}
	inline int has_pair_count(int pair_count) const {
		return pair_counts_.find(pair_count) != pair_counts_.end();
	}
	inline int max_pair_count() const {
		return max_pair_count_;
	}
	inline int edge_count() const {
		return edge_count_;
	}
	inline string wid2word(int wid) const {
		return wid2word_[wid];
	}
	void LoadData();
	void NormalizeTime();
	void UnNormalizeTime();

	vector<User> users_;
  private:
	inline void update_time_range(double time) {
		if (time - 1 < min_time_)
			min_time_ = time - 1;
		if (time + 1 > max_time_)
			max_time_ = time + 1;
	}
	inline void init_max_pair_count() {
		max_pair_count_ = 0;
		edge_count_ = 0;
		for (int i = 0; i < users_.size(); ++i) {
			for (int out_i = 0; out_i < users_[i].out_edges.size(); ++out_i) {
				pair_counts_.insert((int)users_[i].out_edges[out_i].second.times.size());
				max_pair_count_ = max(max_pair_count_, (int)users_[i].out_edges[out_i].second.times.size());
			}
			edge_count_ += users_[i].out_edges.size();
		}
	}
	string in_dir_;
	map<string, int> user2uid_;
	map<string, int> word2wid_;
	vector<string> wid2word_;
	double min_time_;
	double max_time_;
	int max_pair_count_;
	int edge_count_;
	set<int> pair_counts_;
};
#endif