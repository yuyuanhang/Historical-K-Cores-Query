#include "Graph.h"


Graph::Graph() {
    min_k_ = 2;
    idx_size_ = 0;
    query_v_ = nullptr;
    query_deg_ = nullptr;
    // core_t_ = nullptr;
    // v_a_ = nullptr;
    // v_b_ = nullptr;
    // core_ = nullptr;
    // t_offset_ = nullptr;
    // cd_ = nullptr;
    log_f_ = nullptr;
}

Graph::~Graph() {
    // delete[] nbr_cnt_;
    // delete[] core_t_;
    // delete[] core_;
    // delete[] ct_cnt_;
    // delete[] t_offset_;
    // delete[] v_a_;
    // delete[] v_b_;
    // delete[] cd_;
    delete[] query_v_;
    delete[] query_deg_;
    if(log_f_!= nullptr) fclose(log_f_);
}

void Graph::load(const string &path, const int &total_edge_num, bool timestamp_third) {
    g_path_ = path;
    ts_third_ = timestamp_third;
    if(log_f_ != nullptr) fprintf(log_f_,"Graph path: %s, total edge count: %d\n",path.c_str(),total_edge_num);
    //printf("Graph path: %s\n",path.c_str());
    //printf("Loading Graph\n");

    ifstream ifs(path);
    if(!ifs.is_open()){
        cerr << "open file failed!" << endl;
        exit(-1);
    }


    n_ = 0;
    m_ = 0;
    max_deg_ = 0;
    t_ = 0;


    int u,v;
    long ts;
    while(ifs.good() && !ifs.eof()) {
        char line[200];
        ifs.getline(line, 200);
        if (line[0] < '0' || line[0] > '9') continue;
        if (timestamp_third) {
            sscanf(line, "%d %d %ld", &u, &v, &ts);
        } else {
            int weight;
            sscanf(line, "%d %d %d %ld", &u, &v, &weight, &ts);
        }

        if(u == v) continue;

        auto it_ = t_new_to_old_.find(ts);
        if (it_ != t_new_to_old_.end()) {
            t_edges_[it_->second].emplace_back(make_pair(u, v));
        } else {
            t_new_to_old_[ts] = t_new_to_old_.size();
            t_edges_.resize(t_new_to_old_.size()+1);
            t_edges_[t_new_to_old_[ts]].emplace_back(make_pair(u, v));
        }
        // if(max(u, v)+1 > n_) n_ = max(u, v)+1;
    }
    ifs.close();

    t_max_ = t_new_to_old_.size();
    t_ = t_max_ * s_f_;

    edges_idx_.emplace_back(0);
    for (auto i = 0; i < t_; i++) {
        for (auto edge : t_edges_[i]) {
            u = edge.first;
            v = edge.second;

            edges_.emplace_back(make_pair(u,v));

            //adjust size of neighbor list if necessary.
            if(u+1 > nbr_.size()){
                nbr_.resize(u+1);
            }
            if(v+1 > nbr_.size()){
                nbr_.resize(v+1);
            }

//          construct neighbor list
            nbr_[u].emplace_back(make_pair(v,i));
            if(nbr_[u].size() > max_deg_) max_deg_ = nbr_[u].size();

            nbr_[v].emplace_back(make_pair(u,i));
            if(nbr_[v].size() > max_deg_) max_deg_ = nbr_[v].size();

            ++m_;
            if (total_edge_num!=-1 && m_ >= total_edge_num) break;
        }
        edges_idx_.emplace_back(edges_.size());
    }

    n_ = nbr_.size();

//    cout<<"Vertex Num: "<<n_<<", Edge Num: "<<m_<<"; Time Span: "<<t_<<endl;

    init_nbr_cnt();

    //printf("n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
    //printf("span = %d.\n",t_);

    if(log_f_!= nullptr){
        fprintf(log_f_, "n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
        fprintf(log_f_, "span = %d.\n",t_);
    }

    //print_graph_size();

    v_a_.resize(n_);
    v_b_.resize(n_);
    for (int i = 0; i < n_; ++i) {
        v_a_[i] = false;
        v_b_[i] = false;
    }
    query_v_ = new vector<int>[n_];
    query_deg_ = new int[n_];
    t_offset_.resize(n_, 0);
}

void Graph::hcindex() {
    core_t_.resize(n_);
    v_a_.resize(n_);
    v_b_.resize(n_);

    core_decomposition();

    size_t hc_idx_size = 0;
    for (int i = 0; i < n_; i++) {
        for (int j = 2; j < core_t_[i].size(); j++) {
            auto idx = core_t_[i][j];
            for (int k = 0; k < idx.size()-1; k++) {
                auto num_ts = idx[k+1].first - idx[k].first;
                auto num_te = t_max_ - idx[k].second;
                hc_idx_size += num_ts * num_te * sizeof(int);
            }
        }
    }

    if(log_f_!= nullptr){
        fprintf(log_f_, "hc_idx size: %lldMB\n", hc_idx_size/1024/1024);
    }
}

//initialize effective_deg_[], forward_nbr_cnt_[], and backward_nbr_cnt_[]
void Graph::init_nbr_cnt() {

    effective_m_ = 0;

    nbr_cnt_.resize(n_);


    max_effective_deg_ = 0;

    for (int u = 0; u < n_; ++u) {

        for(auto &i : nbr_[u]){
            if (nbr_cnt_[u].find(i.first) != nbr_cnt_[u].end()){
                ++nbr_cnt_[u][i.first];
            } else{
                nbr_cnt_[u].insert(make_pair(i.first,1));
            }
        }

        if(nbr_cnt_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_[u].size();
        effective_m_ += nbr_cnt_[u].size();
    }

    effective_m_ /= 2;
}

void Graph::index() {

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

//  CT
    core_t_.resize(n_);
    v_a_.resize(n_);
    v_b_.resize(n_);


    //printf("starting core decomposition...\n");
//    compute core number for all edges
    core_decomposition();
    for (int u = 0; u < n_; ++u) {
        v_a_[u] = false;
        v_b_[u] = false;
        core_t_[u].resize(core_[u]+1);
    }
    printf("k_max = %d\n",k_max_);

    // printf("initialize core time.\n");
    compute_core_deg(0);
    init_core_time();


//    init ct_deg_ and ct_cnt
    queue<int> q;
    for (int k = min_k_; k <= k_max_; ++k) {
        // printf("Iteration k = %d.\n",k);
        init_ct_cnt(k);

        for (int t_s = 1; t_s < t_; ++t_s) {
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;

                if (invalid(u,k) || invalid(v,k)) continue;

//                process u
                if (!v_a_[u]){
                    del_nbr(u,v);
                    if (ct_cnt_[u].size()<k){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }


//                process v
                if (!v_a_[v]) {
                    del_nbr(v, u);
                    if (ct_cnt_[v].size() < k) {
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            while (!q.empty()){
                int u = q.front();
                q.pop();
                v_a_[u] = false;

                ct_cnt_[u].clear();
                vector<int> nbr_t;
                vector<int> bm_history;


                int ct = 0;

//                  compute new core time of u
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {

                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k && t > ct) break;
                    if (t < t_s){
                        t_offset_[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
//                    v is not valid or has already been visited
                    if (invalid(v,k) || v_b_[v]) continue;

//                    mark v has been visited
                    v_b_[v] = true;
                    int v_t = core_t_[v][k].back().second;
                    nbr_t.emplace_back(max(t,v_t));
                    bm_history.emplace_back(v);

                    if (nbr_t.size() <= k) ct = max(ct,v_t);

                }
                for (auto &v:bm_history) v_b_[v] = false;

                int new_t = t_;
                if (nbr_t.size() >= k){
                    nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
//                    sort(nbr_t.begin(),nbr_t.end());
                    new_t = nbr_t[k-1];
                }

                int old_t = core_t_[u][k].back().second;
                if (core_t_[u][k].back().first == t_s){
                    core_t_[u][k].back().second = new_t;
                }else{
                    core_t_[u][k].emplace_back(make_pair(t_s,new_t));
                }

//                compute ct_cnt_[u] and add neighbor to queue if necessary
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
//                    compute ct_cnt_[u]
                    int t = nbr_[u][i].second;
                    int v = nbr_[u][i].first;
                    if (t > new_t) break;
                    if (invalid(v,k) || core_t_[v][k].back().second > new_t) continue;

                    if (new_t != t_){
                        if (ct_cnt_[u].find(v)==ct_cnt_[u].end() ){
                            ct_cnt_[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt_[u][v];
                        }
                    }

//                    add neighbor to queue if necessary
                    if (v_a_[v]) continue;
                    if (core_t_[v][k].back().second < old_t || new_t <= core_t_[v][k].back().second) continue;
//                    del_nbr(v,u);
                    ct_cnt_[v].erase(u);
                    if (ct_cnt_[v].size() < k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

        }
    }





#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    //printf("Running time: %lld s, %lld mins\n", t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time: %lld s\n",t_msec/1000);


    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    if(log_f_ != nullptr) fprintf(log_f_,"Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
    //printf("Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
#else
    clock_t end = clock();
    //printf("Running time: %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);

#endif
    if(log_f_ != nullptr) fprintf(log_f_,"kmax = %d\n",k_max_);
    //print_idx_size();
}



//core decomposition for all edges
void Graph::core_decomposition() {
    core_.resize(n_);

    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= max_effective_deg_ ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n_; ++u) {
        t_offset_[u] = bin[core_[u]];
        vert[t_offset_[u]] = u;
        bin[core_[u]]++;
    }

    for (int i = max_effective_deg_; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    k_max_ = 0;

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_[u]){
            if (v_a_[item.first]) continue;
            v_a_[item.first] = true;
            if (core_[item.first] > core_[u]){
                int dv = core_[item.first], pv = t_offset_[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    t_offset_[item.first] = pw, vert[pv] = w;
                    t_offset_[w] = pv, vert[pw] = item.first;
                }
                ++bin[dv];
                --core_[item.first];
            }

        }

        for (auto& item : nbr_[u]){
            v_a_[item.first] = false;
        }

        if (core_[u] > k_max_) k_max_ = core_[u];
    }

    delete[] bin;
    delete[] vert;
}

void Graph::init_core_time() {
//    vector<int> core_history;
    int* core = new int[n_];
    for (int u = 0; u < n_; ++u) {
        core[u] = core_[u];
    }

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = t_-1; t_e >= 0 ; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {

            int u = edges_[i].first;
            int v = edges_[i].second;

    //        process u
            if (core[u] <= core[v]){
                --cd_[u][v];
                if (cd_[u][v] == 0){
                    cd_[u].erase(v);
                    if (cd_[u].size() < core[u] && !v_a_[u]){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
            }

//            process v
            if (core[v] <= core[u]){
                --cd_[v][u];
                if (cd_[v][u] == 0){
                    cd_[v].erase(u);
                    if (cd_[v].size() < core[v] && !v_a_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }



        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a_[u] = false;

            int oc = core[u];

            memset(cnt,0,sizeof(int)*(oc+1));

            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (v_b_[v]) continue;
                v_b_[v] = true;

                ++cnt[core[v] < core[u] ? core[v]:core[u]];
            }

            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            int cd = 0;
            for (int k = oc; k >= 0 ; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core[u] = k;
                    break;
                }
            }

//            update cd_;
            cd_[u].clear();
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (core[v] < core[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }


    //        add influenced neighbor to the queue
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                int v = nbr_[u][i].first;
                if (core[u] < core[v] && core[v] <= oc && !v_b_[v]){
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < core[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            for (int k = oc; k > core[u]; --k) {
                core_t_[u][k].emplace_back(make_pair(0, t_e));
            }

        }

    }



    delete[] cnt;
    delete[] core;
}

void Graph::init_ct_cnt(int k) {

//    reuse nbr_cnt_ and avoid applying space
//    ct_cnt_ = nbr_cnt_;

    for (int u = 0; u < n_; ++u) {
        if (invalid(u,k)) continue;
        t_offset_[u] = 0;
        ct_cnt_[u].clear();

        int t = core_t_[u][k].front().second;
        for (auto &i : nbr_[u]){
            if (i.second > t) break;
            int v = i.first;
            if (core_[v]<k || t < core_t_[v][k].front().second) continue;
            if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                ct_cnt_[u].insert(make_pair(v,1));
            }else{
                ++ ct_cnt_[u][v];
            }
        }
    }

}

void Graph::update_dec() {
    update_core_t(); // convert t_{max} to inf_
    init_nbr_cnt(); // recover nbr_cnt_
    rct_cnt_.resize(n_);
    cache_flag_.resize(n_);
    cache_cv_.resize(n_);
    cache_flag_b_.resize(n_);
    cache_mcd_.resize(n_);

    int num_upt = t_max_ - t_;
    auto start = chrono::high_resolution_clock::now();

    int cnt = 0;
    for (auto i = t_; i < t_max_; i++) {
        for (auto edge : t_edges_[i]) {
            cnt++;
            edges_.emplace_back(edge);

            auto u = edge.first;
            auto v = edge.second;

            bool is_expand = false;
            bool is_duplicate = false;

            expand(u, v, i, is_expand, is_duplicate);
            expand(v, u, i, is_expand, is_duplicate);

            if (is_expand) continue;

            for (int k = min_k_; k <= min(core_[u], core_[v]); k++) {
                int t_u = rct(u, k);
                int t_v = rct(v, k);
                if (t_u > t_v) {
                    swap(u, v);
                    swap(t_u, t_v);
                }
                if (t_u == t_+1 || t_u == t_max_-1) continue;
                int k_u = corevalue(t_u, u);
                // can not reach k for both u and v because one edge can only inc core by 1
                if (k_u < k-1) continue;
                // at t_u, core(v)=
                if (quick_check(u, v, t_u)) {
                    find_subcore(u, k-1, t_u);
                    auto k_verts = refine_cand(k-1);
                    if (k_verts.empty()) continue;
                    init_rct_cnt(k, t_u, k_verts);
                    auto km_verts = compute_rct(t_u, k, k_verts, u);
                    for (auto km_vert : km_verts) {
                        int sz = core_t_[km_vert.first][k].size();
                        if (core_t_[km_vert.first][k][sz-2].second == i) { // not a new t_e
                            core_t_[km_vert.first][k].back().first = km_vert.second;
                        } else {
                            core_t_[km_vert.first][k].back().second = i;
                            core_t_[km_vert.first][k].emplace_back(make_pair(km_vert.second, inf_));
                        }
                    }
                }
            }
            if (!is_duplicate) {
                update_core(u, v, i);
            }
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> diff = end - start;
            cout << cnt << endl;
            if (diff.count() >= 43200) { //12 hrs
                // cout << "updating running time per timestamp: " << (double)diff.count()/num_upt << " seconds" << endl;
                fprintf(log_f_, "%d ", cnt);
                fprintf(log_f_,"updating running time: %.2f seconds\n", (double)diff.count());
                exit(0);
            }
        }
        edges_idx_.emplace_back(edges_.size());
        t_ = i+1;
        //if (i % 1000 == 0) test_u();
    }
    fprintf(log_f_, "%d ", cnt);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    fprintf(log_f_,"updating running time: %.2f seconds\n", (double)diff.count());
}

unordered_map<int, int> Graph::compute_rct(int old_inv, int k, set<int> k_verts, int src) {
    vector<int> ch_flag(n_, false);
    set<int> poten_verts = k_verts;
    vector<int> new_inv(n_);
    bool is_finished = false;
    v_b_.assign(n_, false);
    set<int> changed_verts;
    queue<int> q;
    for (auto i = old_inv; i < t_; i++) {
        for (auto edge : t_edges_[i]) {
            int u = edge.first;
            int v = edge.second;

            if (!v_a_[u] && !v_a_[v]) continue;

            if (v_a_[u] && !v_b_[u]) { // the computation of u has not been finished yet
                del_rct(u ,v);
                if (rct_cnt_[u].size() < k) {
                    q.push(u);
                    v_b_[u] = true;
                }
            }

            if (v_a_[v] && !v_b_[v]) {
                del_rct(v ,u);
                if (rct_cnt_[v].size() < k) {
                    q.push(v);
                    v_b_[v] = true;
                }
            }
        }

        for (auto &u : poten_verts) {
            auto rct_copy = rct_cnt_[u];
            for (auto &v : rct_copy) {
                if (!v_a_[v.first] && rct(v.first, k) == i + 1) {
                    rct_cnt_[u].erase(v.first);
                    if (rct_cnt_[u].size() < k && !v_b_[u]) {
                        q.push(u);
                        v_b_[u] = true;
                    }
                }
            }
        }

        while (!q.empty()) {
            int u = q.front(); q.pop();
            poten_verts.erase(u);
            int old_rct = rct(u, k);
            if (old_rct < i+1) {
                ch_flag[u] = true;
                new_inv[u] = i+1;
                changed_verts.insert(u);
            }
            if (old_rct > i+1) {
                cout << "unexcepted new rct" << endl;
                exit(0);
            }
            for (int j = nbr_[u].size()-1; j >= 0; j--) {
                int v = nbr_[u][j].first;
                int t = nbr_[u][j].second;
                if (t <= i) break;
                if (!v_b_[v]&&v_a_[v]) {
                    rct_cnt_[v].erase(u);
                    if (rct_cnt_[v].size() < k) {
                        q.push(v);
                        v_b_[v] = true;
                    }
                    if (ch_flag[u]) {
                        poten_verts.insert(v);
                    }
                }
            }
        }

        if (poten_verts.find(src) == poten_verts.end()) {
            is_finished = true;
            break;
        }
    }

    if (!is_finished) {
        for (int i = edges_idx_[t_]; i < edges_.size(); i++) {
            auto edge = edges_[i];
            int u = edge.first;
            int v = edge.second;

            if (!v_a_[u] || !v_a_[v]) continue;

            if (!v_b_[u]) { // the computation of u has not been finished yet
                del_rct(u ,v);
                if (rct_cnt_[u].size() < k) {
                    q.push(u);
                    v_b_[u] = true;
                }
            }

            if (!v_b_[v]) {
                del_rct(v ,u);
                if (rct_cnt_[v].size() < k) {
                    q.push(v);
                    v_b_[v] = true;
                }
            }
        }
        while (!q.empty()) {
            int u = q.front(); q.pop();
            int old_rct = rct(u, k);
            if (old_rct < t_+1) {
                ch_flag[u] = true;
                new_inv[u] = t_+1;
                changed_verts.insert(u);
            }
            if (old_rct > t_+1) {
                cout << "unexcepted new rct" << endl;
                exit(0);
            }
        }
    }

    unordered_map<int, int> km_verts;
    for (auto &u : changed_verts) {
        km_verts.insert(make_pair(u, new_inv[u]));
    }

    return km_verts;
}

void Graph::update_core(int u, int v, int t) {
    int r = u;
    if (core_[v] < core_[u]) r = v;
    int k = core_[r];

    find_subcore(r, core_[r], 0);
    auto k_verts = refine_cand(core_[r]);
    if (k_verts.empty()) return;
    for (auto k_vert : k_verts) {
        core_[k_vert] = k+1;
        core_t_[k_vert].resize(k+2);
        core_t_[k_vert][k+1].emplace_back(make_pair(0, inf_));
    }
    init_rct_cnt(k+1, 0, k_verts);
    auto km_verts = compute_rct(0, k+1, k_verts, r);
    for (auto km_vert : km_verts) {
        int sz = core_t_[km_vert.first][k+1].size();
        if (sz >= 2) {
            if (core_t_[km_vert.first][k+1][sz-2].second == t) { // not a new t_e
                core_t_[km_vert.first][k+1].back().first = km_vert.second;
            } else {
                core_t_[km_vert.first][k+1].back().second = t;
                core_t_[km_vert.first][k+1].emplace_back(km_vert.second, inf_);
            }
        } else {
            core_t_[km_vert.first][k+1].back().second = t;
            core_t_[km_vert.first][k+1].emplace_back(km_vert.second, inf_);
        }
    }

    if (!k_verts.empty()&&k+1 > k_max_) {
        k_max_ = k+1;
    }
}

void Graph::init_rct_cnt(int k, int t, set<int> k_verts) {
    v_a_.assign(n_, false);
    for (auto &u : k_verts) {
        v_a_[u] = true;
        rct_cnt_[u].clear();

        for (int i = nbr_[u].size()-1; i >= 0; i--){
            int v = nbr_[u][i].first;
            int ts = nbr_[u][i].second;
            if (ts < t) break;
            if (rct(v, k) > t || k_verts.find(v) != k_verts.end()) {
                if (rct_cnt_[u].find(v) == rct_cnt_[u].end()){
                    rct_cnt_[u][v] = 1;
                }else{
                    rct_cnt_[u][v]++;
                }
            }
        }
    }
}

void Graph::update_core_t() {
    for (int i = 0; i < n_; i++) {
        // if (core_[i] == 0) core_[i] = 1;
        for (int j = min_k_; j <= core_[i]; j++) {
            if (core_t_[i][j].back().second == t_) {
                core_t_[i][j].back().second = inf_;
            } else {
                core_t_[i][j].emplace_back(t_, inf_);
            }
            // last item is (t_-1, t_-1)
            // core_t_[i][j].emplace_back(std::make_pair(t_, inf_));
        }
    }
}

void Graph::find_subcore(int r, int k, int t_s) {
    sc_verts_.clear();
    sc_edges_.clear();
    sc_edges_.resize(n_);
    fill(v_a_.begin(), v_a_.end(), false);
    fill(v_b_.begin(), v_b_.end(), false);
    fill(cache_flag_.begin(), cache_flag_.end(), false);
    fill(cache_flag_b_.begin(), cache_flag_b_.end(), false);

    queue<int> q;
    q.push(r);
    v_a_[r] = true;
    while (!q.empty()) {
        int v = q.front(); q.pop();
        sc_verts_[v] = 0;
        for (int i = nbr_[v].size()-1; i >= 0; i--) {
            auto n = nbr_[v][i];
            if (n.second < t_s) break;
            int w = n.first;
            if (!v_b_[w]) {
                int crt_cv;
                if (cache_flag_[w]) { crt_cv = cache_cv_[w]; }
                else {
                    crt_cv = corevalue(t_s, w);
                    cache_flag_[w] = true;
                    cache_cv_[w] = crt_cv;
                }
                //int crt_pcd;
                //if (crt_cv == k) { crt_pcd = pcd(t_s, w, k); }
                //int crt_mcd;
                //if (crt_cv == k) { crt_mcd = mcd(t_s, w, k); }
                if (crt_cv >= k) {
                //if (crt_cv > k||(crt_cv == k&&crt_mcd > k)) {
                //if (crt_cv > k||(crt_cv == k&&crt_pcd > k)) {
                        sc_verts_[v]++;
                        if (crt_cv == k) {
                            sc_edges_[v].emplace_back(w);
                            if (!v_a_[w]) {
                                q.push(w);
                                v_a_[w] = true;
                            }
                        }
                }
                v_b_[w] = true;
            }
        }
        // reset
        fill(v_b_.begin(), v_b_.end(), false);
    }
}

set<int> Graph::refine_cand(int k) {
    set<pair<int, int>, less_core> sorted_verts(sc_verts_.begin(), sc_verts_.end());

    while (!sorted_verts.empty()) {
        pair<int, int> p = *sorted_verts.begin();
        if (p.second <= k) {
            sorted_verts.erase(sorted_verts.begin());
            for (auto n : sc_edges_[p.first]) {
                if (sc_verts_[n] > p.second) {
                    sorted_verts.erase(make_pair(n, sc_verts_[n]));
                    sorted_verts.insert(make_pair(n, sc_verts_[n]-1));
                    auto it = find(sc_edges_[n].begin(), sc_edges_[n].end(), p.first);
                    if (it != sc_edges_[n].end()) {
                        sc_edges_[n].erase(it);
                    }
                    sc_verts_[n]--;
                }
            }
        } else break;
    }

    set<int> k_verts;
    for (auto p : sorted_verts) {
        k_verts.insert(p.first);
    }
    return k_verts;
}

void Graph::update_nbr(int u, int v, int k, vector<int> &q) {
    //if (u == 1688&&v==277&&k==3) {
    if (v == 3&&k == 3) {
        int a = inf_ + 1;
    }
    if (update_f_[v]) {
        if (ordered_nbr_[v].size() < k) return;
        ordered_nbr_[v].modify(u, ub_[u]);
        if (ub_[v] == ordered_nbr_[v].kvalue(k)) return;
    } else {
        for (int i = nbr_[v].size()-1; i >= 0; i--) {
            int n = nbr_[v][i].first;
            int tb = nbr_[v][i].second;
            if (tb < rct_[v]-1) break;
            if (v_a_[n]) continue;
            v_a_[n] = true;
            int tv;
            if (u == n) tv = ub_[n];
            else tv = pre_rct_[n]-1;
            if (tv < rct_[v]-1) continue;
            tv = min(tb, tv);
            tv = max(tv, -1);
            ordered_nbr_[v].insert(n, tv);
        }
        fill(v_a_.begin(), v_a_.end(), false);
    }
    if (!update_f_[v]) {
        update_f_[v] = true;
    }
    ub_[v] = ordered_nbr_[v].kvalue(k);
    if (ub_[v] < rct_[v]-1) {
        cout << "smaller rct" << endl;
    }
    if (ub_[v] >= rct(v, k-1)&&ub_[v] != -1) {
        cout << "beyond ub" << endl;
    }
    if (find(q.begin(), q.end(), v) == q.end()) q.push_back(v);
};

vector<int> Graph::init_queue(vector<int> cand_v, int k, int t) {
    ub_.assign(n_, -1);
    update_f_.assign(n_, false);
    ordered_nbr_.clear();
    ordered_nbr_.resize(n_);

    vector<int> q;
    /*for (auto u : cand_v) {
        for (int i = nbr_[u].size()-1; i >= 0; i--) {
            int v = nbr_[u][i].first;
            int tb = nbr_[u][i].second;
            if (tb < rct_[u]-1) break;
            if (v_a_[v]) continue;
            v_a_[v] = true;
            int tv;
            tv = pre_rct_[v]-1;
            if (tv < rct_[u]-1) continue;
            tv = min(tb, tv);
            tv = max(tv, -1);
            ordered_nbr_[u].insert(v, tv);
        }
        fill(v_a_.begin(), v_a_.end(), false);
        ub_[u] = ordered_nbr_[u].kvalue(k);
        q.push_back(u);
        update_f_[u] = true;
    }
    while (!q.empty()) {
        int u = q.front(); q.erase(q.begin());
        for (int i = nbr_[u].size()-1; i >= 0; i--) {
            int v = nbr_[u][i].first;
            int tb = nbr_[u][i].second;
            if (tb < t) break;
            if (v_a_[v]) continue;
            v_a_[v] = true;
            update_nbr(u, v, k, q);
        }
        fill(v_a_.begin(), v_a_.end(), false);
    }*/
    for (auto u : cand_v) {
        q.push_back(u);
    }
    return q;
}

void Graph::init_rct(int k) {
    rct_.resize(n_);
    for (int i = 0; i < n_; i++) {
        rct_[i] = rct(i, k);
    }
}

vector<int> Graph::collect(vector<int> cand_v) {
    vector<bool> v_c(n_, false);
    vector<int> c;
    queue<int> q;
    for (auto u : cand_v) {
        v_c[u] = true;
        if (ub_[u] >= max(0, rct_[u])) {
            q.push(u); c.push_back(u);
        }
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto p : ordered_nbr_[u].m) {
            int v = p.first;
            if (v_a_[v]) continue;
            v_a_[v] = true;
            if (v_c[v]) continue;
            v_c[v] = true;
            if (ub_[v] >= max(0, rct_[v])) {
                q.push(v); c.push_back(v);
            }
        }
        fill(v_a_.begin(), v_a_.end(), false);
    }
    return c;
}

void Graph::update_inc() {
    update_core_t(); // convert t_{max} to inf_
    init_nbr_cnt(); // recover nbr_cnt_
    fill(v_a_.begin(), v_a_.end(), false);
    fill(v_b_.begin(), v_b_.end(), false);
    // rct_cnt_.resize(n_);

    int num_upt = t_max_ - t_;
    auto start = chrono::high_resolution_clock::now();
    int num_visited = 0;
    int num_updated = 0;
    for (int i = t_; i < t_max_; i++) {
        cout << "updating " << i << endl;
        vector<int> cand_v;

        for (auto edge: t_edges_[i]) {
            edges_.emplace_back(edge);

            auto u = edge.first;
            auto v = edge.second;

            bool is_expand = false;
            bool is_duplicate = false;

            expand(u, v, i, is_expand, is_duplicate);
            expand(v, u, i, is_expand, is_duplicate);
            if (!is_expand) {
                if (find(cand_v.begin(), cand_v.end(), u) == cand_v.end()) cand_v.push_back(u);
                if (find(cand_v.begin(), cand_v.end(), v) == cand_v.end()) cand_v.push_back(v);
            }
        }
        if (cand_v.empty()) continue;

        int k_bound = 0;
        for (auto u : cand_v) {
            int core_ub = core_[u];
            for (int j = nbr_[u].size()-1; j >= 0; j--) {
                int v = nbr_[u][j].first;
                int t = nbr_[u][j].second;
                if (t < i) break;
                if (core_[v] >= core_[u]) {
                    core_ub++;
                }
            }
            k_bound = max(k_bound, core_ub);
        }

        pre_rct_.resize(n_);
        for (int j = 0; j < n_; j++) { pre_rct_[j] = rct(j, 1); }

        ub_.resize(n_);
        int k;
        for (k = min_k_; k <= k_bound; k++) {
            init_rct(k);
            // auto q = init_queue(cand_v, k, i);
            for (int j = 0; j < n_; j++) {
                ub_[j] = pre_rct_[j] - 1;
                ub_[j] = max(-1, ub_[j]);
            }
            update_f_.assign(n_, false);
            queue<int> q;
            for (auto u : cand_v) {
                q.push(u);
                v_a_[u] = true;
            }

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                v_a_[u] = false;
                vector<int> nbr_t;
                vector<int> bm_history;
                int ct = inf_;

                // compute new core time of u
                for (int j = nbr_[u].size()-1; j >= 0; j--) {
                    int t = nbr_[u][j].second;
                    int v = nbr_[u][j].first;
                    if (nbr_t.size() >= k && t < ct) break;
                    // mark v has been visited
                    if (v_b_[v]) continue;

                    v_b_[v] = true;
                    int v_t = ub_[v];
                    nbr_t.emplace_back(min(t, v_t));
                    bm_history.emplace_back(v);
                    if (nbr_t.size() <= k) ct = min(ct, v_t);
                }
                for (auto &v:bm_history) v_b_[v] = false;

                int new_t = -1;
                if (nbr_t.size() >= k){
                    nth_element(nbr_t.begin(), nbr_t.begin()+nbr_t.size()-k, nbr_t.end());
                    new_t = nbr_t[nbr_t.size()-k];
                }

                if (new_t == ub_[u]&&update_f_[u]) {
                    continue;
                } else {
                    update_f_[u] = true;
                    ub_[u] = new_t;
                }

                // add neighbor to queue if necessary
                for (int j = nbr_[u].size()-1; j >= 0; j--) {
                    int t = nbr_[u][j].second;
                    int v = nbr_[u][j].first;
                    if (!v_a_[v]) {
                        v_a_[v] = true;
                        q.push(v);
                    }
                }
            }
            for (int j = 0; j < n_; j++) {
                if (update_f_[j]) { num_visited++; }
                if (ub_[j] >= max(0, rct_[j])&&update_f_[j]) {
                    if (core_[j] < k) {
                        core_[j] = k;
                        core_t_[j].resize(k+1);
                        core_t_[j][k].push_back(make_pair(0, inf_));
                    }
                    core_t_[j][k].back().second = i;
                    core_t_[j][k].push_back(make_pair(ub_[j]+1, inf_));
                    rct_[j] = ub_[j]+1;
                    num_updated++;
                }
            }
            pre_rct_ = rct_;
            /*while (!q.empty()) {
                int u = q.front(); q.erase(q.begin());
                for (int j = nbr_[u].size() - 1; j >= 0; j--) {
                    int v = nbr_[u][j].first;
                    int t = nbr_[u][j].second;
                    if (v_b_[v]) continue;
                    v_b_[v] = true;
                    if (t >= rct_[v]) {
                        update_nbr(u, v, k, q);
                    }
                }
                fill(v_b_.begin(), v_b_.end(), false);
            }
            auto c = collect(cand_v);
            for (auto u : c) {
                if (core_[u] < k) {
                    core_[u] = k;
                    core_t_[u].resize(k+1);
                    core_t_[u][k].push_back(make_pair(0, inf_));
                }
                core_t_[u][k].back().second = i;
                core_t_[u][k].push_back(make_pair(ub_[u]+1, inf_));
                rct_[u] = ub_[u]+1;
            }
            pre_rct_ = rct_;*/
        }

        edges_idx_.emplace_back(edges_.size());
        t_ = i+1;
        //if (i % 10 == 0) { test_u(); }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    // cout << "updating running time per timestamp: " << (double)diff.count()/num_upt << " seconds" << endl;
    cout << "updating running time: " << (double)diff.count() << " seconds" << endl;
    cout << "the number of visited vertices: " << num_visited << endl;
    cout << "the number of changed vertices: " << num_updated << endl;
}

void Graph::update_baseline() {
    update_core_t();
    init_nbr_cnt();

    auto start = chrono::high_resolution_clock::now();

    for (auto i = t_; i < t_max_; i++) {
        // G = G \cup E_{i}
        t_++;
        edges_idx_.resize(t_ + 1);
        edges_idx_[t_] = edges_idx_[t_ - 1];

        for (auto edge: t_edges_[i]) {
            edges_.emplace_back(edge);

            auto src = edge.first;
            auto dest = edge.second;

            auto expand = [&](int u, int v) {
                if(u + 1 > nbr_.size()){
                    n_ = u + 1;
                    nbr_.resize(n_);
                    nbr_cnt_.resize(n_);
                }

                nbr_[u].emplace_back(make_pair(v,i));

                if (nbr_cnt_[u].find(v) != nbr_cnt_[u].end()) { // (u, v) is not a new edge
                    nbr_cnt_[u].find(v)->second++;
                } else {
                    nbr_cnt_[u][v] = 1;
                    if(nbr_cnt_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_[u].size();
                }
            };

            expand(src, dest);
            expand(dest, src);
            edges_idx_[t_]++;

            // prepare parameters for core decompose
            t_offset_.assign(n_, 0);
            core_.resize(n_, 0);
            v_a_.assign(n_, false);

            // core decompose
            vector<int> old_core(core_);
            core_decomposition();

            // initialize index for new vertices/core values
            core_t_.resize(n_);
            for (auto u = 0; u < n_; u++) {
                if (old_core[u] == core_[u]) continue;
                core_t_[u].resize(core_[u] + 1);
                for (auto k = max(min_k_, old_core[u] + 1); k <= core_[u]; k++) core_t_[u][k].emplace_back(0, inf_);
            }
            // prepare parameters
            compute_core_deg(0);
            v_a_.assign(n_, false);
            v_b_.assign(n_, false);
            // compute rct
            vector<int> core(core_);
            queue<int> q;
            int* cnt = new int[k_max_+1];
            for (auto t_s = 0; t_s < t_; t_s++) {
                for (auto t = edges_idx_[t_s]; t < edges_idx_[t_s+1]; t++) {
                    int u = edges_[t].first;
                    int v = edges_[t].second;

                    auto process = [&](int u, int v) {
                        if (core[u] <= core[v]){
                            --cd_[u][v];
                            if (cd_[u][v] == 0){
                                cd_[u].erase(v);
                                if (cd_[u].size() < core[u] && !v_a_[u]){
                                    q.push(u);
                                    v_a_[u] = true;
                                }
                            }
                        }
                    };

                    process(u, v);
                    process(v, u);
                }

                while (!q.empty()){
                    int u = q.front();
                    q.pop();
                    v_a_[u] = false;
                    int oc = core[u];
                    memset(cnt, 0,sizeof(int) * (oc + 1));
                    // compute local core
                    for (int j = nbr_[u].size() - 1; j >= 0; j--) {
                        int v = nbr_[u][j].first;
                        int t = nbr_[u][j].second;
                        if (t <= t_s) break;
                        if (v_b_[v]) continue;
                        v_b_[v] = true;
                        ++cnt[core[v] < core[u] ? core[v]:core[u]];
                    }
                    int cd = 0;
                    for (int k = oc; k >= 0 ; --k) {
                        cd += cnt[k];
                        if(cd >= k){
                            core[u] = k;
                            break;
                        }
                    }
                    // recover flag
                    for (int j = nbr_[u].size()-1; j >= 0; j--) {
                        if (nbr_[u][j].second <= t_s) break;
                        v_b_[nbr_[u][j].first] = false;
                    }
                    // update cd_;
                    cd_[u].clear();
                    for (int j = nbr_[u].size()-1; j >= 0; j--) {
                        int v = nbr_[u][j].first;
                        int t = nbr_[u][j].second;
                        if (t <= t_s) break;
                        if (core[v] < core[u]) continue;
                        if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                        else ++cd_[u][v];
                    }
                    // add influenced neighbor to the queue
                    for (int j = nbr_[u].size()-1; j >= 0; j--) {
                        if (nbr_[u][j].second <= t_s) break;
                        int v = nbr_[u][j].first;
                        if (core[u] < core[v] && core[v] <= oc && !v_b_[v]){
                            v_b_[v] = true;
                            cd_[v].erase(u);
                            if (!v_a_[v] && cd_[v].size() < core[v]){
                                q.push(v);
                                v_a_[v] = true;
                            }
                        }
                    }
                    // recover flag
                    for (int j = nbr_[u].size()-1; j >= 0; j--) {
                        if (nbr_[u][j].second <= t_s) break;
                        v_b_[nbr_[u][j].first] = false;
                    }
                    // update index
                    if (oc < min_k_) continue;
                    for (int k = oc; k > max(core[u], min_k_-1); k--) {
                        auto old_rct = rct(u, k);
                        if (old_rct != t_s+1) {
                            int sz = core_t_[u][k].size();
                            if (sz >= 2) {
                                if (core_t_[u][k][sz-2].second == t_ - 1) { // not a new t_e
                                    core_t_[u][k].back().first = t_s+1;
                                } else {
                                    core_t_[u][k].back().second = t_ - 1;
                                    core_t_[u][k].emplace_back(t_s+1, inf_);
                                }
                            } else {
                                core_t_[u][k].back().second = t_ - 1;
                                core_t_[u][k].emplace_back(t_s+1, inf_);
                            }
                        }
                    }
                }
            }
            delete[] cnt;
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> diff = end - start;
            if (diff.count() >= 43200) { //12 hrs
                fprintf(log_f_,"updating running time: %.2f seconds\n", diff.count());
                exit(0);
            }
        }
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    if(log_f_ != nullptr) fprintf(log_f_,"updating running time: %.2f seconds\n", (double)diff.count());
}

void Graph::test() {

    for (int u = 0; u < 20; ++u) {
        printf("vertex %d:\n",u);

        for (int k = min_k_; k < core_t_[u].size(); ++k) {
            printf("k=%d:\n",k);
            for (auto &i:core_t_[u][k]) printf("[%d,%d], ",i.first,i.second);
            printf("\n");
        }
        printf("\n");
    }

//    auto fp = fopen(R"(C:\Users\DW\Desktop\idx\email-b)","rb");
//    int c,d;
//    fread(&c,sizeof(unsigned int),1,fp);
//
//
//    for (int u = 0; u < n_; ++u) {
//        int cs;
//        fread(&cs,sizeof(int),1,fp);
//        if (cs != core_t_[u].size()) printf("neq %d\n",u);
//
//        for (int k = min_k_; k < cs; ++k) {
//            int ccs;
//            fread(&ccs,sizeof(int),1,fp);
//            if (ccs != core_t_[u][k].size()) printf("neq ccs %d[%d]\n",u,k);
//            for (int i = 0; i < ccs; ++i) {
//                int a,b;
//                fread(&a,sizeof(int),1,fp);
//                fread(&b,sizeof(int),1,fp);
//
//            }
//        }
//    }
//    fclose(fp);

}

void Graph::remove_expired_t(int t) {
    update_core_t();
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < 1; i++) {
        auto core_t_back = core_t_;
        for (int u = 0; u < n_; ++u) {
            for (int k = min_k_; k < core_t_back[u].size(); k++) {
                auto it = std::upper_bound(core_t_back[u][k].begin(), core_t_back[u][k].end(), std::make_pair(t, std::numeric_limits<int>::min()),
                                           [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                                               return a.first < b.first;
                                           });

                if (it != core_t_back[u][k].begin()) {
                    auto target = std::prev(it);
                    if (target->first <= t) {
                        target->first = t;
                        core_t_back[u][k].erase(core_t_back[u][k].begin(), target);
                    } else {
                        std::cout << "removal error.\n";
                        exit(0);
                    }
                } // else no need to change index
                if (core_t_back[u][k].begin()->second == inf_) {
                    core_t_back[u].erase(core_t_back[u].begin()+k, core_t_back[u].end());
                    break;
                }
            }
        }
    }
    auto end = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if(log_f_ != nullptr) fprintf(log_f_,"removing expired timestamps: %.2f microseconds\n", (double)duration.count());
}

void Graph::print_idx_size() {
    if (idx_size_ != 0) printf("Index size: %lld MB.",idx_size_/1024/1024);

    idx_size_ += sizeof(int);
    idx_size_ += sizeof(int)*n_;

    double average_t = 0;
    int average_t_d = 0;
    int max_t = 0;

    for (int u = 0; u < n_; ++u) {
        if (core_t_[u].size() <= min_k_) continue;
        for (int k = min_k_; k < core_t_[u].size(); ++k){
            idx_size_ += core_t_[u][k].size()*2*sizeof(int);
            average_t += core_t_[u][k].size();
            if (core_t_[u][k].size() > max_t) max_t = core_t_[u][k].size();
        }
        average_t_d += core_t_[u].size()-min_k_;
    }

    printf("Index size: %.2f MB.\n",(float)idx_size_/1024/1024);
    printf("Average T = %.2f, max T = %d.\n",average_t/average_t_d,max_t);

    if(log_f_ != nullptr) fprintf(log_f_,"Index size: %.2f MB\n",(float)idx_size_/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Average T = %.2f, max T = %d.\n",average_t/average_t_d,max_t);
}

void Graph::print_graph_size() {
    printf("Graph size: %.2f MB.\n",(float)edges_.size()*3*sizeof(int)/1024/1024);
    fprintf(log_f_,"Graph size: %.2f MB.\n",(float)edges_.size()*3*sizeof(int)/1024/1024);
}

bool Graph::query(int u, int t_s, int t_e, int k) {
    if (core_t_[u].size() < k+1) return false;
    if (k < 2){
        printf("Enter a parameter larger than 1\n");
        return false;
    }
    auto it = upper_bound(core_t_[u][k].begin(),core_t_[u][k].end(),make_pair(t_s,t_e),cmp);
    --it;

    return it->second <= t_e;
}

void Graph::write_idx(const string &path) {
    auto fp = fopen(path.c_str(),"wb");
    fwrite(&n_,sizeof(unsigned int),1,fp);

    // debug
    for (int u = 0; u < n_; ++u) {
        if (log_f_ != nullptr) fprintf(log_f_, "vertex %d\n", u);
        int cs = core_t_[u].size();
        if (log_f_ != nullptr) fprintf(log_f_, "core_value %d\n", cs);

        for (int k = min_k_; k < cs; ++k) {
            for (auto &i:core_t_[u][k]){
                if (log_f_ != nullptr) fprintf(log_f_, "[%d, %d]", i.first, i.second);
            }
            if (log_f_ != nullptr) fprintf(log_f_, "\n");
        }
    }

    for (int u = 0; u < n_; ++u) {
        int cs = core_t_[u].size();
        fwrite(&cs,sizeof(int),1,fp);
        for (int k = min_k_; k < cs; ++k) {
            int ccs = core_t_[u][k].size();
            fwrite(&ccs,sizeof(int),1,fp);
            for (auto &i:core_t_[u][k]){
                fwrite(&i.first,sizeof(int),1,fp);
                fwrite(&i.second,sizeof(int),1,fp);
            }
        }
    }
    fclose(fp);
}

void Graph::load_idx(const string &path) {
    auto fp = fopen(path.c_str(),"rb");
    fread(&n_,sizeof(unsigned int),1,fp);
    k_max_ = 0;
    //t_ = 0;

    core_t_.resize(n_);
    core_.resize(n_);

    for (int u = 0; u < n_; ++u) {
        int cs;
        fread(&cs,sizeof(int),1,fp);
        core_t_[u].resize(cs);
        core_[u] = cs-1;
        if(k_max_ < cs-1) k_max_ = cs-1;

        for (int k = min_k_; k < cs; ++k) {
            int ccs;
            fread(&ccs,sizeof(int),1,fp);

            for (int i = 0; i < ccs; ++i) {
                int a,b;
                fread(&a,sizeof(int),1,fp);
                fread(&b,sizeof(int),1,fp);
                core_t_[u][k].emplace_back(make_pair(a,b));
            }
            //if (t_ < core_t_[u][k].back().second) t_ = core_t_[u][k].back().second;
        }
    }
    fclose(fp);
    print_idx_size();
}

void Graph::index_baseline() {
#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    v_a_.resize(n_);
    v_b_.resize(n_);
    if (core_t_.empty()) core_t_.resize(n_);
    else{
        printf("Index structure exists!\n");
        exit(1);
    }
    if (t_offset_.empty()) t_offset_.resize(n_);

    core_decomposition();
    for (int u = 0; u < n_; ++u) {
        v_a_[u] = false;
        v_b_[u] = false;
        core_t_[u].resize(core_[u]+1);
    }
    t_offset_.assign(n_, 0);
    compute_core_deg(0);


    for (int t_s = 0; t_s < t_; ++t_s) {

        if (t_s%100 == 0) printf("t = %d.\n",t_s);

//        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
//            int u = edges_[i].first;
//            int v = edges_[i].second;
//
//            --nbr_cnt_[u][v];
//            if (nbr_cnt_[u][v] == 0) nbr_cnt_[u].erase(v);
//
//            --nbr_cnt_[v][u];
//            if (nbr_cnt_[v][u] == 0) nbr_cnt_[v].erase(u);
//        }
//        test_core_decomposition(t_s);

        compute_core_time_bl(t_s);
        if (t_s == t_-1) break;
        compute_core_deg(t_s);
        decremental_core_bl(t_s);
    }


#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time (baseline): %lld s, %lld mins\n", t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time (Baseline): %lld s\n",t_msec/1000);

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);

#else
    clock_t end = clock();
    printf("Running time (baseline): %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);
#endif

    print_idx_size();
}



void Graph::compute_core_time_bl(const int &t_s) {
    vector<int> core_history;
    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = t_-1; t_e >= t_s ; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {

            int u = edges_[i].first;
            int v = edges_[i].second;

            //        process u
            if (core_[u] <= core_[v]){
                --cd_[u][v];
                if (cd_[u][v] == 0){
                    cd_[u].erase(v);
                    if (cd_[u].size() < core_[u] && !v_a_[u]){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
            }

//            process v
            if (core_[v] <= core_[u]){
                --cd_[v][u];
                if (cd_[v][u] == 0){
                    cd_[v].erase(u);
                    if (cd_[v].size() < core_[v] && !v_a_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }



        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a_[u] = false;

            int oc = core_[u];

            memset(cnt,0,sizeof(int)*(oc+1));

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].second;
                if (t < t_s){
                    t_offset_[u] = i+1;
                    continue;
                }
                int v = nbr_[u][i].first;
                if (t >= t_e) break;
                if (v_b_[v]) continue;
                v_b_[v] = true;
                ++cnt[core_[v] < core_[u] ? core_[v]:core_[u]];
            }

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            int cd = 0;
            for (int k = oc; k >= 0 ; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core_[u] = k;
                    break;
                }
            }

//            update cd_;
            cd_[u].clear();
            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (core_[v] < core_[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }


            //        add influenced neighbor to the queue
            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                int v = nbr_[u][i].first;
                if (core_[u] < core_[v] && core_[v] <= oc && !v_b_[v]){
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < core_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            for (int k = oc; k > core_[u]; --k) {

                core_history.emplace_back(u);
                if (t_s == 0 || core_t_[u][k].empty() || core_t_[u][k].back().second < t_e) {
                    core_t_[u][k].emplace_back(make_pair(t_s, t_e));
                }
            }

        }

    }

    delete[] cnt;
    for (auto &i:core_history) ++core_[i];

}

void Graph::decremental_core_bl(const int &t_s) {

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int i = edges_idx_[t_s]; i < edges_idx_[t_s + 1]; ++i) {

        int u = edges_[i].first;
        int v = edges_[i].second;

        //        process u
        if (core_[u] <= core_[v]) {
            --cd_[u][v];
            if (cd_[u][v] == 0) {
                cd_[u].erase(v);
                if (cd_[u].size() < core_[u] && !v_a_[u]) {
                    q.push(u);
                    v_a_[u] = true;
                }
            }
        }

//            process v
        if (core_[v] <= core_[u]) {
            --cd_[v][u];
            if (cd_[v][u] == 0) {
                cd_[v].erase(u);
                if (cd_[v].size() < core_[v] && !v_a_[v]) {
                    q.push(v);
                    v_a_[v] = true;
                }
            }
        }

    }
    while (!q.empty()){
        int u = q.front();
        q.pop();
        v_a_[u] = false;

        int oc = core_[u];

        memset(cnt,0,sizeof(int)*(oc+1));

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            int t = nbr_[u][i].second;
            if (t <= t_s) break;

            int v = nbr_[u][i].first;
            if (v_b_[v]) continue;
            v_b_[v] = true;
            ++cnt[core_[v] < core_[u] ? core_[v]:core_[u]];
        }

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            v_b_[nbr_[u][i].first] = false;
        }

        int cd = 0;
        for (int k = oc; k >= 0 ; --k) {
            cd += cnt[k];
            if(cd >= k){
                core_[u] = k;
                break;
            }
        }

//            update cd_;
        cd_[u].clear();
        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            int v = nbr_[u][i].first;
            int t = nbr_[u][i].second;
            if (t <= t_s) break;
            if (core_[v] < core_[u]) continue;
            if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
            else ++cd_[u][v];
        }


        //        add influenced neighbor to the queue
        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            int v = nbr_[u][i].first;
            if (core_[u] < core_[v] && core_[v] <= oc && !v_b_[v]){
                v_b_[v] = true;
                cd_[v].erase(u);
                if (!v_a_[v] && cd_[v].size() < core_[v]){
                    q.push(v);
                    v_a_[v] = true;
                }
            }
        }

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            v_b_[nbr_[u][i].first] = false;
        }

        for (int k = oc; k > core_[u]; --k) {
            core_t_[u][k].emplace_back(make_pair(t_s + 1, t_));
        }

    }

    delete[] cnt;

}

void Graph::online_core_decomposition(const int &t_s, const int &t_e) {

    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;


    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (v < u) swap(u,v);
        if (edge_mp.find(u) != edge_mp.end() && edge_mp[u].find(v) != edge_mp[u].end()) continue;

//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int span_max_deg = 0;
    for (auto &i:span_nbr){
        if (span_max_deg < i.size()) span_max_deg = i.size();
    }


    int n = new_to_old.size();
    int* span_core = new int[n];
    int* pos = new int[n];
    auto vert = new int[n];
    auto bin = new int[span_max_deg+1];
    memset(bin,0,sizeof(int)*(span_max_deg+1));

    for (int u = 0; u < n; ++u) {
        int d = span_nbr[u].size();
        span_core[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= span_max_deg ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n; ++u) {
        pos[u] = bin[span_core[u]];
        vert[pos[u]] = u;
        bin[span_core[u]]++;
    }

    for (int i = span_max_deg; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;


    for (int i = 0; i < n; ++i) {
        int u = vert[i];

        for (auto& item : span_nbr[u]){

            if (span_core[item] > span_core[u]){
                int dv = span_core[item], pv = pos[item];
                int pw = bin[dv], w = vert[pw];
                if (item != w){
                    pos[item] = pw, vert[pv] = w;
                    pos[w] = pv, vert[pw] = item;
                }
                ++bin[dv];
                --span_core[item];
            }

        }
    }


//    span_core stores the core number in the current span

    delete[] bin;
    delete[] vert;
    delete[] pos;
    delete[] span_core;
}

void Graph::compute_core_deg(const int &t_s) {
    cd_.resize(n_);

    for (int u = 0; u < n_; ++u) {
        cd_[u].clear();
        for (int i = nbr_[u].size()-1;i>=0;--i){
            int v = nbr_[u][i].first;
            int t = nbr_[u][i].second;
            if (t < t_s) break;
            if (core_[v] < core_[u]) continue;

            if (cd_[u].find(v) == cd_[u].end()){
                cd_[u].insert(make_pair(v,1));
            }else{
                ++cd_[u][v];
            }
        }
    }
}

void Graph::test_core_decomposition(const int &t_s) {

    if (core_.empty()) core_.resize(n_);

    vector<int> old_core(n_);
    for (int u = 0; u < n_; ++u) {
        old_core[u] = core_[u];
    }

    auto pos = new int[n_];
    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= max_effective_deg_ ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n_; ++u) {
        pos[u] = bin[core_[u]];
        vert[pos[u]] = u;
        bin[core_[u]]++;
    }

    for (int i = max_effective_deg_; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    k_max_ = 0;

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_[u]){
            if (item.second < t_s) continue;
            if(v_a_[item.first]) continue;
            v_a_[item.first] = true;

            if (core_[item.first] > core_[u]){
                int dv = core_[item.first], pv = pos[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    pos[item.first] = pw, vert[pv] = w;
                    pos[w] = pv, vert[pw] = item.first;
                }
                ++bin[dv];
                --core_[item.first];
            }

        }

        for (auto& item : nbr_[u]){
            v_a_[item.first] = false;
        }

        if (core_[u] > k_max_) k_max_ = core_[u];
    }

    delete[] bin;
    delete[] vert;
    delete[] pos;

    for (int u = 0; u < n_; ++u) {
        for (int k = old_core[u]; k > core_[u]; --k) {
            core_t_[u][k].emplace_back(make_pair(t_s,t_));
        }
    }
}

void Graph::naive_index() {

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif
    long long idx_size = sizeof(int);

    for (int t_s = 0; t_s < t_; ++t_s) {
        if (t_s % 100 == 0) printf("t_s = %d\n",t_s);
        for (int t_e = t_s; t_e < t_; ++t_e) {
            online_core_decomposition(t_s,t_e);
            idx_size += n_*sizeof(int);
        }
    }

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time (Naive Index): %lld s, %lld mins\n", t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time (Naive Index): %lld s\n",t_msec/1000);
#else
    clock_t end = clock();
    printf("Running time (naive index): %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);
#endif

    printf("Index size: %.2f MB.\n",(float)idx_size/1024/1024);

}

void Graph::query_init() {
    if (v_a_.empty()) v_a_.resize(core_t_.size());
    if (v_b_.empty()) v_b_.resize(core_t_.size());
}

void Graph::query_subgraph(int u, int t_s, int t_e, int k, vector<int>& r, vector<pair<int,int>>& r_edges) {

    if (v_a_.empty()) query_init();

//    vector<int> r;
//    vector<pair<int,int>> r_edges;

    queue<int> q;
    if (query(u,t_s,t_e,k)){
        q.push(u);
        v_a_[u] = true;
        r.emplace_back(u);
    }
    vector<int> bm_history;

    while (!q.empty()){
        int v = q.front();
        q.pop();
        auto nbr_it = lower_bound(nbr_[v].begin(),nbr_[v].end(),make_pair(0,t_s),cmp_nbr);
        while (nbr_it != nbr_[v].end()){
            if (nbr_it->second > t_e) break;
            int w = nbr_it->first;
            ++nbr_it;
            if (v_b_[w]) continue;
            if (v_a_[w]){
                r_edges.emplace_back(make_pair(v,w));
                continue;
            }

            if (query(w,t_s,t_e,k)){
                q.push(w);
                v_a_[w] = true;
                r.emplace_back(w);
                r_edges.emplace_back(make_pair(v,w));
            }else{
                v_b_[w] = true;
                bm_history.emplace_back(w);
            }
        }

    }
    for (auto &i:bm_history) v_b_[i] = false;
    for (auto &i:r) v_a_[i] = false;

//    r is the result array



}

void Graph::init_log(const string &log_path) {
    log_f_ = fopen(log_path.c_str(),"a");
    fprintf(log_f_,"\n\n==================\n");
    time_t now = time(0);
    fprintf(log_f_,"%s\n",ctime(&now));
}

int Graph::query_all(int t_s, int t_e, int k) {
    int r = 0;
    for (int i = 0; i < n_; ++i) {
        if(query(i,t_s,t_e,k)) ++r;
    }
    return r;
}

void Graph::load(const string &path, bool timestamp_third) {
    load(path,-1,timestamp_third);
}

void Graph::naive_index_size() {

    long long idx_size = sizeof(int);

    for (int t_s = 0; t_s < t_; ++t_s) {
//        if (t_s % 100 == 0) printf("t_s = %d\n",t_s);
        for (int t_e = t_s; t_e < t_; ++t_e) {
//            online_core_decomposition(t_s,t_e);
            idx_size += n_*sizeof(int);
        }
    }

    printf("Index size (Naive Index): %.2f MB.\n",(float)idx_size/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Index size (Naive Index): %.2f MB.\n",(float)idx_size/1024/1024);
}

int Graph::online_k_core(const int &t_s, const int &t_e, const int& k) {
    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;


    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (v < u) swap(u,v);
        if (edge_mp.find(u) != edge_mp.end() && edge_mp[u].find(v) != edge_mp[u].end()) continue;

//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}

int Graph::online_query(const int &t_s, const int &t_e, const int &k) {
    int sm;
    return online_query(t_s,t_e,k,sm);
}

int Graph::online_query(const int &t_s, const int &t_e, const int &k, int& snapshot_m) {
    vector<int> vertices;

//    reset snapshot edge number
    snapshot_m = 0;

    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (!v_a_[u]){
            vertices.emplace_back(u);
            v_a_[u] = true;
        }
        if (!v_a_[v]){
            vertices.emplace_back(v);
            v_a_[v] = true;
        }
        query_v_[u].emplace_back(v);
        query_v_[v].emplace_back(u);
    }
    int r = vertices.size();

    queue<int> q;
    for(auto &i : vertices){
        query_deg_[i] = 0;
        vector<int> vb_history;
        for (auto &nbr: query_v_[i]){
            if (v_b_[nbr]) continue;
            v_b_[nbr] = true;
            vb_history.emplace_back(nbr);
            ++ query_deg_[i];


        }
        for(auto &nbr:vb_history) v_b_[nbr] = false;

        snapshot_m += query_deg_[i];

        if (query_deg_[i] < k){
            q.push(i);
            v_a_[i] = false;
            --r;
        }

    }

    snapshot_m /= 2;

    while (!q.empty()){
        int u = q.front();
        q.pop();
        vector<int> vb_history;
        for(auto &i:query_v_[u]){
            if (!v_a_[i] || v_b_[i]) continue;
            v_b_[i] = true;
            vb_history.emplace_back(i);
            -- query_deg_[i];
            if (query_deg_[i] < k){
                v_a_[i] = false;
                --r;
                q.push(i);
            }
        }
        for(auto &nbr:vb_history) v_b_[nbr] = false;
    }

    for (auto &u:vertices){
        query_v_[u].clear();
        v_a_[u] = false;
    }

    return r;
}

int Graph::online_span_core(const int &t_s, const int &t_e, const int &k) {
    vector<vector<pair<int,int>>> interval_edges;
    vector<pair<int,int>> intersection;


    for (int t = t_s; t <= t_e; ++t) {
        interval_edges.emplace_back(vector<pair<int,int>>());
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; ++i){
            int u = edges_[i].first;
            int v = edges_[i].second;
            if (v < u) swap(u,v);
            interval_edges.back().emplace_back(make_pair(u,v));
        }

        if (interval_edges.back().size() < k) return 0;
    }
    if (t_s == t_e){
        intersection = interval_edges.back();
    }else{
        edge_intersection(interval_edges,intersection);
        if (intersection.size() < k) return 0;
    }





    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;

    for (auto &edge:intersection){
        int u = edge.first;
        int v = edge.second;
//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}

void Graph::edge_intersection(vector<vector<pair<int, int>>> &edges, vector<pair<int, int>> &result) {
    vector<int> pos(edges.size(),0);
    pair<int,int> cur = make_pair(-1,-1);


    for (int i = 0; i < edges.size(); ++i){
        if (pos[i] == edges[i].size()) break;
        if (cur.first == -1){
            cur.first = edges[i][pos[i]].first;
            cur.second = edges[i][pos[i]].second;
            continue;
        }
        if (edges[i][pos[i]] < cur) {
            pos[i]++;
            i--;
            continue;
        }

        if (edges[i][pos[i]] > cur) {
            cur.first = edges[i][pos[i]].first;
            cur.second = edges[i][pos[i]].second;
            i = -1;
            continue;
        }

        if (i == edges.size()-1){
            result.emplace_back(cur);
            cur.first = -1;
            pos[0] ++;
            i = -1;
        }

    }


}

int Graph::index_span_core(const int &t_s, const int &t_e, const int& k) {
    vector<vector<pair<int,int>>> interval_edges;
    vector<pair<int,int>> intersection;

    vector<int> bm_history;

    for (int t = t_s; t <= t_e; ++t) {
        interval_edges.emplace_back(vector<pair<int,int>>());
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; ++i){
            int u = edges_[i].first;
            int v = edges_[i].second;
            if (v < u) swap(u,v);

            bool edge_available = true;
            if(v_a_[u]){
                edge_available = false;
            }else if (!query(u,t,t,k)){
                v_a_[u] = true;
                bm_history.emplace_back(u);
                edge_available = false;
            }

            if(v_a_[v]){
                edge_available = false;
            }else if (!query(v,t,t,k)){
                v_a_[v] = true;
                bm_history.emplace_back(v);
                edge_available = false;
            }

            if(edge_available) interval_edges.back().emplace_back(make_pair(u,v));
        }

        if (interval_edges.back().size() < k){
            for(auto &item:bm_history) v_a_[item] = false;
            return 0;
        }
    }
    for(auto &item:bm_history) v_a_[item] = false;

    if (t_s == t_e){
        intersection = interval_edges.back();
    }else{
        edge_intersection(interval_edges,intersection);
        if (intersection.size() < k) return 0;
    }



    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;

    for (auto &edge:intersection){
        int u = edge.first;
        int v = edge.second;
//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}

void Graph::test_u() {
    cout << "testing" << endl;
    auto *g = new Graph();
    g->s_f_ = (float)t_ / t_max_;
    g->load(g_path_, ts_third_);

    g->index();
    cout << "groundtruth done" << endl;

    g->update_core_t();

    bool is_correct = true;
    if (g->n_ != n_) {
        is_correct = false;
        cout << "the number of vertices is wrong" << endl;
    }
    for (int i = 0; i < n_; i++) {
        if (g->core_[i] != core_[i]) {
            is_correct = false;
            cout << "core value error at " << i << " " << g->core_[i] << " " << core_[i] << endl;
        }
    }
    for (int i = 0; i < n_; i++) {
        for (int j = min_k_; j < min(core_t_[i].size(), g->core_t_[i].size()); j++) {
            if (g->core_t_[i][j].size() != core_t_[i][j].size()) {
                is_correct = false;
                cout << "size error at " << i << " " << j << " " << g->core_t_[i][j].size() << " " << core_t_[i][j].size() << endl;
                cout << "index error at " << i << " " << j << " " << " (" <<
                     g->core_t_[i][j].back().first << ", " << g->core_t_[i][j].back().second << ") (" <<
                     core_t_[i][j].back().first << ", " << core_t_[i][j].back().second << ")" << endl;
            } else {
                for (int k = 0; k < core_t_[i][j].size(); k++) {
                    if (g->core_t_[i][j][k].first != core_t_[i][j][k].first ||
                        g->core_t_[i][j][k].second != core_t_[i][j][k].second) {
                        is_correct = false;
                        cout << "index error at " << i << " " << j << " " << k << " (" <<
                             g->core_t_[i][j][k].first << ", " << g->core_t_[i][j][k].second << ") (" <<
                             core_t_[i][j][k].first << ", " << core_t_[i][j][k].second << ")" << endl;
                    }
                }
            }
        }
    }

    delete g;
    cout << "testing done" << endl;
    if (!is_correct) exit(0);
}

void Graph::count_graph_size(long long &snapshot_size, long long &window_size, int t_s, int t_e) {
    snapshot_size = 0;
    window_size = 0;
    nbr_cnt_.clear();
    nbr_cnt_.resize(n_);
    for (auto i = t_s; i <= t_e; i++) {
        window_size += t_edges_[i].size();
        for (auto &e : t_edges_[i]) {
            int u = e.first;
            int v = e.second;

            if (nbr_cnt_[u].find(v) != nbr_cnt_[u].end()){
                ++nbr_cnt_[u][v];
                ++nbr_cnt_[v][u];
            } else{
                nbr_cnt_[u].insert(make_pair(v,1));
                nbr_cnt_[v].insert(make_pair(u,1));
            }
        }
    }
    for (int u = 0; u < n_; u++) {
        snapshot_size += nbr_cnt_[u].size();
    }
    snapshot_size /= 2;
}

void Graph::print_idx() {
    for (int i = 0; i < n_; i++) {
        cout << "vertex " << i << endl;
        for (int j = min_k_; j < core_t_[i].size(); j++) {
            cout << "k " << j << endl;
            for (int k = 0; k < core_t_[i][j].size(); k++) {
                cout << "(" << core_t_[i][j][k].first << ", " << core_t_[i][j][k].second << ")" << endl;
            }
        }
    }
}

bool cmp(const pair<int,int> &a, const pair<int,int> &b){
    return a.first < b.first;
}

bool cmp_nbr(const pair<int,int> &a, const pair<int,int> &b){
    return a.second < b.second;
}

/*bool cmp_ts(const tuple<int,int,int> &a, const tuple<int,int,int> &b) {
    return get<2>(a) < get<2>(b);
}*/
