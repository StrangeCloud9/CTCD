#include "dataset.h"

void Dataset::LoadData() {
	ifstream fin;
	string line, tmps;
	int uid, uid1, uid2, wid;
	double time1, time2;
	Document doc;
	Edge edge;
	cout << "start to load data in " << in_dir_ << endl;
	fin.open(in_dir_ + "/docs.txt");
	int doc_count = 0;
	while(getline(fin, line, '\n')) {
		++doc_count;
		stringstream ss(line);
		ss >> tmps;
		if(user2uid_.find(tmps) == user2uid_.end()){
			user2uid_.insert(pair<string, int>(tmps, user2uid_.size()));
			users_.push_back(User());
		}
		uid = user2uid_[tmps];
		doc.clear();
		ss >> doc.time;
		update_time_range(doc.time);
		while(ss >> tmps){
			if(word2wid_.find(tmps) == word2wid_.end()){
				word2wid_.insert(pair<string, int>(tmps, word2wid_.size()));
				wid2word_.push_back(tmps);
			}
			wid = word2wid_[tmps];
			doc.add_word(wid);
		}
		doc.calc_counter();
		users_[uid].add_doc(doc);
	}
	fin.close();
	cout << doc_count << " documents loaded" << endl;

	int edge_count = 0;
	fin.open(in_dir_ + "/links.txt");
	while(getline(fin, line, '\n')) {
		++edge_count;
		stringstream ss(line);
		ss >> tmps;
		if(user2uid_.find(tmps) == user2uid_.end()){
			user2uid_.insert(pair<string, int>(tmps, user2uid_.size()));
			users_.push_back(User());
		}
		uid1 = user2uid_[tmps];
		ss >> tmps;
		if(user2uid_.find(tmps) == user2uid_.end()){
			user2uid_.insert(pair<string, int>(tmps, user2uid_.size()));
			users_.push_back(User());
		}
		uid2 = user2uid_[tmps];
		edge.clear();
		while(ss >> time1 >> time2){
			update_time_range(time1);
			update_time_range(time2);
			edge.add_time_pair(pair<double, double>(time1, time2));
		}
		users_[uid1].add_edge(uid2, edge);
	}
	init_max_pair_count();
	cout << edge_count << " edges loaded" << endl;
}

void Dataset::NormalizeTime() {
	for(int i = 0; i < users_.size(); ++i) {
		for(int j = 0; j < users_[i].docs.size(); ++j){
			users_[i].docs[j].time = normalize_one_time(users_[i].docs[j].time);
		}
		for(int j = 0; j < users_[i].out_edges.size(); ++j){
			for(int k = 0; k < users_[i].out_edges[j].second.times.size(); ++k){
				users_[i].out_edges[j].second.times[k].first = normalize_one_time(users_[i].out_edges[j].second.times[k].first);
				users_[i].out_edges[j].second.times[k].second = normalize_one_time(users_[i].out_edges[j].second.times[k].second);
			}
		}
	}
}

void Dataset::UnNormalizeTime() {
	for(int i = 0; i < users_.size(); ++i) {
		for(int j = 0; j < users_[i].docs.size(); ++j){
			users_[i].docs[j].time = unnormalize_one_time(users_[i].docs[j].time);
		}
		for(int j = 0; j < users_[i].out_edges.size(); ++j){
			for(int k = 0; k < users_[i].out_edges[j].second.times.size(); ++k){
				users_[i].out_edges[j].second.times[k].first = unnormalize_one_time(users_[i].out_edges[j].second.times[k].first);
				users_[i].out_edges[j].second.times[k].second = unnormalize_one_time(users_[i].out_edges[j].second.times[k].second);
			}
		}
	}
}