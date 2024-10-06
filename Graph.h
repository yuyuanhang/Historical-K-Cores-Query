//
// Created by DW on 7/29/2020.
//

#ifndef SPAN_CORE_GRAPH_H
#define SPAN_CORE_GRAPH_H


#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstring>
#include <bitset>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <set>
#include <map>
#include <list>
#include <limits>
#include <tuple>
#include <chrono>

#define _LINUX_

#ifdef _LINUX_
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#endif

using namespace std;

bool cmp(const pair<int,int> &a, const pair<int,int> &b);
bool cmp_nbr(const pair<int,int> &a, const pair<int,int> &b);
// bool cmp_ts(const tuple<int,int,int> &a, const tuple<int,int,int> &b);

struct less_core {
    bool operator()(const pair<int,int> &a, const pair<int,int> &b) {
        if (a.second != b.second) {
            return cmp_nbr(a, b);
        } else return cmp(a, b);
    }
};

struct pair_order {
    bool operator()(const pair<int, int>& a, const pair<int, int>& b) const {
        if (a.second != b.second) {
            return a.second > b.second;
        }
        return a.first < b.first;
    }
};

class OrderedPairs {
public:
    map<int, int> m;
    set<pair<int, int>, pair_order> s;

    void insert(int key, int value) {
        m[key] = value;
        s.insert({key, value});
    }

    void modify(int key, int newValue) {
        auto it = m.find(key);
        // always get a minimum value, otherwise we cannot get a min(t_edge, t_valid)
        if (it != m.end()) {
            if (it->second <= newValue) return;
            s.erase({key, it->second});
            it->second = newValue;
            s.insert({key, newValue});
        }
    }

    bool key_exist(int key) {
        if (m.find(key) != m.end()) return true;
        else return false;
    }

    int kvalue(int k) {
        if (s.size() < k) {
            return -1;
        }
        auto it = s.begin();
        advance(it, k - 1);
        return it->second;
    }

    int size() { return s.size(); }
};

class Graph {

    unsigned int n_{};
    unsigned int m_{};
    unsigned int max_deg_{};
    unsigned int effective_m_{};
    int min_k_;
    long long idx_size_;

    unsigned int max_effective_deg_{};

    FILE* log_f_;

    unordered_map<long, int> t_new_to_old_;
    vector<pair<int,int>> edges_;
    vector<int> edges_idx_;
    vector<vector<pair<int,int>>> nbr_;
//    vector<vector<int>> reverse_nbr_idx_;

    vector<unordered_map<int,int>> nbr_cnt_;

    vector<int> core_{};

    vector<int>* query_v_{};
    int* query_deg_{};

    vector<unordered_map<int,int>> cd_;


    vector<bool> v_a_;
    vector<bool> v_b_;


    vector<unordered_map<int,int>>& ct_cnt_ = nbr_cnt_;

    vector<int> t_offset_;

    // new parameters
    vector<vector<pair<int, int>>> t_edges_;
    unordered_map<int, int> sc_verts_;
    vector<vector<int>> sc_edges_;
    int inf_ = numeric_limits<int>::max();
    string g_path_;
    bool ts_third_;
    vector<unordered_map<int,int>> rct_cnt_;
    vector<bool> cache_flag_;
    vector<int> cache_cv_;
    vector<bool> cache_flag_b_;
    vector<int> cache_mcd_;
    vector<int> ub_;
    vector<int> rct_;
    vector<int> pre_rct_;
    vector<bool> update_f_;
    vector<OrderedPairs> ordered_nbr_;

    void init_nbr_cnt();
    void init_core_time();
    void init_ct_cnt(int k);

    void del_nbr(int u, int v);
    bool invalid(int u, int k);

//    functions for baseline index construction
    void decremental_core_bl(const int &t_s);
    void compute_core_time_bl(const int &t_s);
    void compute_core_deg(const int &t_s);

    void test_core_decomposition(const int &t_s);
    void print_idx_size();
    void print_graph_size();

//  functions for updating
    void update_core_t();
    int corevalue(int t_s, int u);
    void find_subcore(int r, int k, int t_s);
    set<int> refine_cand(int k);
    int rct(int u, int k);
    void expand(int u, int v, int t, bool &is_expand, bool &is_duplicate);
    unordered_map<int, int> compute_rct(int old_inv, int k, set<int> k_verts, int src);
    void update_core(int u, int v, int t);
    bool quick_check(int u, int v, int t_s);
    void init_rct_cnt(int k, int t, set<int> k_verts);
    void del_rct(int u, int v);
    int mcd(int t_s, int u, int k);
    int pcd(int t_s, int u, int k);
    void update_nbr(int u, int v, int k, vector<int> &q);
    void init_rct(int k);
    vector<int> init_queue(vector<int> cand_v, int k, int t);
    vector<int> collect(vector<int> cand_v);

public:
    unsigned int t_{};
    unsigned int t_max_{};
    float s_f_ = 1;
    int k_max_{};
    vector<vector<vector<pair<int,int>>>> core_t_{};

    Graph();
    ~Graph();
    void load(const string &path, bool timestamp_third);
    void load(const string &path, const int &total_edge_num, bool timestamp_third);

    bool query(int u, int t_s, int t_e, int k);
    int query_all(int t_s, int t_e, int k);
    void query_subgraph(int u, int t_s, int t_e, int k, vector<int>& r, vector<pair<int,int>>& r_edges);

    void query_init();
    void core_decomposition();

    void hcindex();
    void index();
    void index_baseline();
    void update_dec();
    void update_inc();
    void update_baseline();
    void remove_expired_t(int t);

    void load_idx(const string &path);
    void write_idx(const string &path);
    void init_log(const string &log_path);

    void naive_index();
    void naive_index_size();

    void test();
    void test_u();


//    void online_span_core(const int &u, const int &t_s, const int &t_e);
    void online_core_decomposition(const int &t_s, const int &t_e);
    int online_k_core(const int &t_s, const int &t_e, const int& k);
    int online_query(const int &t_s, const int &t_e, const int& k, int& snapshot_m);
    int online_query(const int &t_s, const int &t_e, const int& k);

//    for case study
    int online_span_core(const int &t_s, const int &t_e, const int& k);
    void edge_intersection(vector<vector<pair<int,int>>>& edges, vector<pair<int,int>>& result);
    int index_span_core(const int &t_s, const int &t_e, const int& k);

    void count_graph_size(long long &snapshot_size, long long &window_size, int t_s, int t_e);
    void print_idx();
};




inline void Graph::del_nbr(int u, int v) {
//    if (core_[u] < k) return;
    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) return;
    --ct_cnt_[u][v];
    if (ct_cnt_[u][v]==0) ct_cnt_[u].erase(v);
}

inline bool Graph::invalid(int u, int k) {
    if (core_[u] < k) return true;
    return (!core_t_[u][k].empty()) && core_t_[u][k].back().second == t_;
}

inline int Graph::rct(int u, int k) {
    if (k > core_[u]) { return -1; }
    if (k == 1) {
        if ( nbr_[u].empty() ) return -1;
        else return nbr_[u].back().second+1;
    }
    return core_t_[u][k].back().first;
}

inline int Graph::corevalue(int t_s, int u) {
    int k = 1; // we do not access a point that is not in the graph
    for (int i = min_k_; i <= core_[u]; i++) {
        if (core_t_[u][i].back().first > t_s) k = i;
        else break;
    }
    return k;
}

inline int Graph::mcd(int t_s, int u, int k) {
    int mcd = 0;
    for (int i = nbr_[u].size()-1; i >= 0; i--) {
        if (nbr_[u][i].second < t_s) break;
        int crt_cv;
        int w = nbr_[u][i].first;
        if (cache_flag_[w]) { crt_cv = cache_cv_[w]; }
        else {
            crt_cv = corevalue(t_s, w);
            cache_flag_[w] = true;
            cache_cv_[w] = crt_cv;
        }
        if (crt_cv >= k) mcd++;
    }
    return mcd;
}

inline int Graph::pcd(int t_s, int u, int k) {
    int pcd = 0;
    for (int i = nbr_[u].size()-1; i >= 0; i--) {
        if (nbr_[u][i].second < t_s) break;
        int crt_cv;
        if (cache_flag_[u]) { crt_cv = cache_cv_[u]; }
        else {
            crt_cv = corevalue(t_s, u);
            cache_flag_[u] = true;
            cache_cv_[u] = crt_cv;
        }
        int crt_mcd;
        int w = nbr_[u][i].first;
        if (cache_flag_b_[w]) { crt_mcd = cache_mcd_[w]; }
        else {
            crt_mcd = mcd(t_s, w, k);
            cache_flag_b_[w] = true;
            cache_mcd_[w] = crt_mcd;
        }
        if (crt_cv > k||(crt_cv==k&&crt_mcd>k)) pcd++;
    }
    return pcd;
}

inline void Graph::expand(int u, int v, int t, bool &is_expand, bool &is_duplicate) {
    if (u + 1 > nbr_.size()) {
        n_ = u+1;
        nbr_.resize(n_);
        nbr_cnt_.resize(n_);
        core_.resize(n_, 0);
        core_[u] = 1;
        v_a_.resize(n_, false);
        v_b_.resize(n_, false);
        sc_edges_.resize(n_);
        core_t_.resize(n_);
        core_t_[u].resize(2);
        rct_cnt_.resize(n_);
        cache_flag_.resize(n_, false);
        cache_cv_.resize(n_);
        cache_flag_b_.resize(n_, false);
        cache_mcd_.resize(n_);
        ub_.resize(n_);
        is_expand = true;
    }
    if (nbr_[u].empty()) core_[u] = 1;
    nbr_[u].emplace_back(make_pair(v, t));
    if (nbr_cnt_[u].find(v) != nbr_cnt_[u].end()) { // (u, v) is not a new edge
        nbr_cnt_[u].find(v)->second++;
        is_duplicate = true;
    } else {
        nbr_cnt_[u][v] = 1;
        if (nbr_cnt_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_[u].size();
    }
}

inline bool Graph::quick_check(int u, int v, int t_s) {
    // ignore the last neighbor because it must be v
    for (int i = nbr_[u].size()-2; i >= 0; i--) {
        if (nbr_[u][i].second < t_s) return true;
        if (nbr_[u][i].first == v) return false;
    }
    return true;
}

inline void Graph::del_rct(int u, int v) {
    if (rct_cnt_[u].find(v) == rct_cnt_[u].end()) return;
    --rct_cnt_[u][v];
    if (rct_cnt_[u][v]==0) rct_cnt_[u].erase(v);
}

#endif //SPAN_CORE_GRAPH_H
