#include "Graph.h"

#include <iostream>
#include <random>
#include <cmath>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {

    if(strcmp(argv[1],"-idx") == 0 || strcmp(argv[1],"-idx-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);

        auto *g = new Graph();
        g->s_f_ = 1.0;
        g->init_log(log_path);


        if (argc > 5){
            g->load(graph_path,atoi(argv[5]),strcmp(argv[1],"-idx") == 0);
        }else{
            g->load(graph_path,strcmp(argv[1],"-idx") == 0);
        }

        g->index();
        g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-idx-bl") == 0 || strcmp(argv[1],"-idx-bl-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        auto *g = new Graph();
        g->init_log(log_path);

        if (argc > 5){
            g->load(graph_path,atoi(argv[5]),strcmp(argv[1],"-idx-bl") == 0);
        }else{
            g->load(graph_path,strcmp(argv[1],"-idx-bl") == 0);
        }

        g->index_baseline();
        g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-idx-size") == 0 || strcmp(argv[1],"-idx-size-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        float ratio = 1;
        if (argc > 5) {
            ratio = atof(argv[5]);
        }
        auto *g = new Graph();
        g->s_f_ = ratio;
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-idx-size") == 0);

        g->index();
        g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-hc-idx-size") == 0 || strcmp(argv[1],"-hc-idx-size-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        auto *g = new Graph();
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-hc-idx-size") == 0);
        g->load_idx(idx_path);

        g->hcindex();
        delete g;
    }
//---------------------------------------------------------------------------
// query-original
//---------------------------------------------------------------------------

    if(strcmp(argv[1],"-q") == 0 || strcmp(argv[1],"-q-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        int t_range = atoi(argv[5]);
//        int k_range = atoi(argv[6]);
        int cnt = atoi(argv[6]);

        auto *g = new Graph();
//        g->init_log(log_path);
        g->load(graph_path,strcmp(argv[1],"-q") == 0);
        g->load_idx(idx_path);


        int t_gap = g->t_*t_range/100-1;

        printf("Prepare %d queries\n",cnt);

        vector<pair<int,int>> queries;
        int t_s_max = g->t_ - t_gap - 1;

        random_device rd;
        uniform_int_distribution<int> dist(0,t_s_max-1);

        for (int i = 0; i < cnt; ++i) {
            int t_s = 0;
            if (t_s_max != 0) t_s = dist(rd);
            int t_e = t_s+t_gap;
            queries.emplace_back(make_pair(t_s,t_e));
        }

        FILE* log_f = fopen(log_path.c_str(),"a");
        fprintf(log_f,"\n\n:====Query Processing====\n");
        time_t now = time(0);
        fprintf(log_f,"%s\n",ctime(&now));
        fprintf(log_f,"Index path:%s\n",idx_path.c_str());
        fprintf(log_f,"t_range: %d, t_max = %d, k_max = %d,\n",t_range,g->t_,g->k_max_);
	    printf("t_range: %d, t_max = %d, k_max = %d,\n",t_range,g->t_,g->k_max_);

        for (int k_range = 2; k_range < 10 ; k_range+=2){
            int k = g->k_max_*k_range/10;
            fprintf(log_f,"-------\n");
            fprintf(log_f,"k_range = %d,k = %d\n",k_range,k);
            printf("k_range = %d,k = %d\n",k_range,k);

            int empty_q = 0;
            int result_size = 0;

    #ifdef _LINUX_
            struct timeval t_start,t_end;
            gettimeofday(&t_start, NULL);
    #endif

            for (auto &ti:queries){
                int s = g->query_all(ti.first,ti.second,k);
                result_size += s;
                if(s <= 0) ++empty_q;
    //            printf("cnt:%d\n",g->query_all(ti.first,ti.second,k));
            }

    #ifdef _LINUX_
            gettimeofday(&t_end, NULL);
            long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);

            float average_result_size = float(result_size)/cnt;
            printf("Average query time: %lld *e-6 s\n", t_msec/cnt);
            printf("Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);
            fprintf(log_f,"Average query time: %lld *e-6 s\n", t_msec/cnt);
            fprintf(log_f,"Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);
    #endif

            fprintf(log_f,"--------\n");
            empty_q = 0;
            result_size = 0;

            long long average_snapshot_m = 0;

    #ifdef _LINUX_
            gettimeofday(&t_start, NULL);
    #endif

            for (auto &ti:queries){
                int snapshot_edges;
                int s = g->online_query(ti.first,ti.second,k,snapshot_edges);
                average_snapshot_m += snapshot_edges;
                result_size += s;
                if(s <= 0) ++empty_q;
    //            printf("cnt:%d\n",g->online_query(ti.first,ti.second,k));

            }


    #ifdef _LINUX_
            gettimeofday(&t_end, NULL);
            t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);

            average_result_size = float(result_size)/cnt;
            printf("Average online query time: %lld *e-6 s\n", t_msec/cnt);
            printf("Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);
            fprintf(log_f,"Average online query time: %lld *e-6 s\n", t_msec/cnt);
            fprintf(log_f,"Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);

            long long average_snapshot_m2 = average_snapshot_m/(long long)cnt;
            printf("Average snapshot edges number: %lld\n", average_snapshot_m2);
            fprintf(log_f,"Average snapshot edges number: %lld\n", average_snapshot_m2);
    #endif


        }
        delete g;
        fclose(log_f);
    }


//---------------------------------------------------------------------------
// query - compute standard devision for running time
//---------------------------------------------------------------------------
  
    if(strcmp(argv[1],"-q-sd") == 0 || strcmp(argv[1],"-q-sd-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        int t_range = atoi(argv[5]);
        int k_range = atoi(argv[6]);
        int cnt = atoi(argv[7]);

        auto *g = new Graph();
//        g->init_log(log_path);
        g->load(graph_path,strcmp(argv[1],"-q-sd") == 0);
        g->load_idx(idx_path);


        int t_gap = g->t_*t_range/100-1;

        printf("Prepare %d queries\n",cnt);

        vector<pair<int,int>> queries;
        vector<long long> times;
        long long snapshot_size = 0;
        long long window_size = 0;

        int t_s_max = g->t_ - t_gap - 1;

        random_device rd;
        uniform_int_distribution<int> dist(0,t_s_max-1);

        for (int i = 0; i < cnt; ++i) {
            int t_s = 0;
            if (t_s_max != 0) t_s = dist(rd);
            int t_e = t_s+t_gap;
            cout << t_s << " " << t_e << endl;
            queries.emplace_back(make_pair(t_s,t_e));
        }

        FILE* log_f = fopen(log_path.c_str(),"a");
        fprintf(log_f,"\n\n:====Query Processing Standard Deviation====\n");
        time_t now = time(0);
        fprintf(log_f,"%s\n",ctime(&now));
        fprintf(log_f,"Index path:%s\n",idx_path.c_str());
        fprintf(log_f,"t_range: %d, t_max = %d, k_max = %d,\n",t_range,g->t_,g->k_max_);

        //for (int k_range = 6; k_range <= 6; k_range+=2){
            int k = g->k_max_*k_range/10;
            fprintf(log_f,"-------\n");
            fprintf(log_f,"k_range = %d,k = %d\n",k_range,k);
            printf("k_range = %d,k = %d\n",k_range,k);



            for (auto &ti:queries){
#ifdef _LINUX_
                struct timeval t_start,t_end;
                gettimeofday(&t_start, NULL);
#endif
                g->query_all(ti.first,ti.second,k);

#ifdef _LINUX_
                gettimeofday(&t_end, NULL);
                long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);
                times.emplace_back(t_msec);
#endif
            }

//            compute average time
            float average_time = 0;
            for (long long &i:times) average_time += i;
            average_time /= times.size();

//            compute standard deviation
            float sd = 0;
            for (long long &i:times){
                sd += pow(average_time-i,2);
            }
            sd = sqrt(sd/times.size());

//            float average_result_size = float(result_size)/cnt;

            printf("Average query time: %.2f *e-6 s\n", average_time);
            printf("Standard deviation: %.2f *e-6 s\n", sd);
//                printf("Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);
            fprintf(log_f,"Average query time: %.2f *e-6 s\n", average_time);
            fprintf(log_f,"Standard deviation: %.2f *e-6 s\n", sd);



            fprintf(log_f,"--------\n");


            times.clear();


            for (auto &ti:queries){
#ifdef _LINUX_
                struct timeval t_start,t_end;
                gettimeofday(&t_start, NULL);
#endif

                g->online_query(ti.first,ti.second,k);

#ifdef _LINUX_
                gettimeofday(&t_end, NULL);
                long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);
                times.push_back(t_msec);
#endif

            }

            //            compute average time
            average_time = 0;
            for (long long &i:times) average_time += i;
            average_time /= times.size();

//            compute standard deviation
            sd = 0;
            for (long long &i:times){
                sd += pow(average_time-i,2);
            }
            sd = sqrt(sd/times.size());

//            float average_result_size = float(result_size)/cnt;

            printf("Average online query time: %.2f *e-6 s\n", average_time);
            printf("Standard deviation: %.2f *e-6 s\n", sd);
//                printf("Empty result queries: %d, average result size: %.2f\n",empty_q,average_result_size);
            fprintf(log_f,"Average online query time: %.2f *e-6 s\n", average_time);
            fprintf(log_f,"Standard deviation: %.2f *e-6 s\n", sd);

            times.clear();

        //}
        delete g;
        fclose(log_f);
    }

    if(strcmp(argv[1],"-snap-size") == 0 || strcmp(argv[1],"-snap-size-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        int t_range = atoi(argv[5]);
//        int k_range = atoi(argv[6]);
        int cnt = atoi(argv[6]);

        auto *g = new Graph();
//        g->init_log(log_path);
        g->load(graph_path,strcmp(argv[1],"-q-sd") == 0);
        g->load_idx(idx_path);


        int t_gap = g->t_*t_range/100-1;

        printf("Prepare %d queries\n",cnt);

        vector<pair<int,int>> queries;
        vector<long long> times;
        long long snapshot_size = 0;
        long long window_size = 0;

        int t_s_max = g->t_ - t_gap - 1;

        random_device rd;
        uniform_int_distribution<int> dist(0,t_s_max-1);

        for (int i = 0; i < cnt; ++i) {
            int t_s = 0;
            if (t_s_max != 0) t_s = dist(rd);
            int t_e = t_s+t_gap;
            queries.emplace_back(make_pair(t_s,t_e));
        }

        FILE* log_f = fopen(log_path.c_str(),"a");
        fprintf(log_f,"\n\n:====Snapshot Size====\n");
        for (auto &ti:queries){
            long long crt_window = 0;
            long long crt_snapshot = 0;
            g->count_graph_size(crt_snapshot, crt_window, ti.first, ti.second);
            snapshot_size += crt_snapshot;
            window_size += crt_window;
        }
        printf("Average window size: %d \n", window_size/cnt);
        printf("Average snapshot size: %d \n", snapshot_size/cnt);
        fprintf(log_f,"Average window size: %d \n", window_size/cnt);
        fprintf(log_f,"Average snapshot size: %d \n", snapshot_size/cnt);
    }
//---------------------------------------------------------------------------
// query - what is q search?
//---------------------------------------------------------------------------

    if(strcmp(argv[1],"-q-search") == 0 || strcmp(argv[1],"-q-search-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        int t_s = atoi(argv[5]);
        int t_e = atoi(argv[6]);
        int k = atoi(argv[7]);
        int query_id = atoi(argv[8]);


        auto *g = new Graph();

        g->load(graph_path,strcmp(argv[1],"-q-search") == 0);
        g->load_idx(idx_path);

        FILE* log_f = fopen(log_path.c_str(),"a");
        fprintf(log_f,"\n\n:====Query Search for Case Study====\n");
        time_t now = time(0);
        fprintf(log_f,"%s\n",ctime(&now));
        fprintf(log_f,"Index path:%s\n",idx_path.c_str());
        fprintf(log_f,"Interval=[%d,%d], k = %d, query_id = %d\n",t_s,t_e,k,query_id);

        vector<int> subgraph_vertices;
        vector<pair<int,int>> subgraph_edges;
        g->query_subgraph(query_id, t_s, t_e, k, subgraph_vertices, subgraph_edges);

        if (subgraph_vertices.empty()){
            printf("No subgraph found\n");
            fprintf(log_f,"result size: 0\n");
        } else{
            printf("Subgraph size: %lu\n",subgraph_vertices.size());
            fprintf(log_f,"Subgraph size: %lu\n",subgraph_vertices.size());
            fprintf(log_f,"Vertices:\n");
            for (int &i:subgraph_vertices) fprintf(log_f,"%d, ",i);
            fprintf(log_f,"\n\nEdges:\n");
            for (auto &i:subgraph_edges) fprintf(log_f,"(%d,%d), ",i.first,i.second);
            fprintf(log_f,"\n");
        }

        delete g;
        fclose(log_f);
    }


//---------------------------------------------------------------------------
// query - span core case
//---------------------------------------------------------------------------

    if(strcmp(argv[1],"-case-span-core") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        int t_s = atoi(argv[5]);
        int t_e = atoi(argv[6]);
        int k = atoi(argv[7]);

        auto *g = new Graph();
        g->load(graph_path,true);
        g->load_idx(idx_path);

        if (t_s >= g->t_ || t_e >= g->t_){
            printf("time is too large!\n");
            delete g;
            return 0;
        }

        FILE* log_f = fopen(log_path.c_str(),"a");
        fprintf(log_f,"\n\n:====Case Study Span Core====\n");
        time_t now = time(0);
        fprintf(log_f,"%s\n",ctime(&now));
        fprintf(log_f,"query interval = [%d,%d], k = %d\n",t_s,t_e,k);


#ifdef _LINUX_
        struct timeval t_start,t_end;
        gettimeofday(&t_start, NULL);
#endif

        int result_size = g->index_span_core(t_s,t_e,k);

#ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);
        printf("Index-based span core query time: %lld *e-6 s\n", t_msec);
        printf("Index-based Result size: %d\n",result_size);
//        fprintf(log_f,"Index-based span core query time: %lld *e-6 s\n", t_msec);
//        fprintf(log_f,"Result size: %d\n",result_size);
#endif


#ifdef _LINUX_
        gettimeofday(&t_start, NULL);
#endif

        result_size = g->online_span_core(t_s,t_e,k);

#ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        t_msec = (t_end.tv_sec - t_start.tv_sec)*1000*1000 + (t_end.tv_usec - t_start.tv_usec);
        printf("Online span core query time: %lld *e-6 s\n", t_msec);
        printf("Online result size: %d\n",result_size);
//        fprintf(log_f,"Index-based span core query time: %lld *e-6 s\n", t_msec);
//        fprintf(log_f,"Result size: %d\n",result_size);
#endif



        delete g;
        fclose(log_f);
    }

    if(strcmp(argv[1],"-uidx") == 0 || strcmp(argv[1],"-uidx-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        float ratio = 1;
        if (argc > 5) {
            ratio = atof(argv[5]);
        }

        auto *g = new Graph();
        g->s_f_ = ratio;
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-uidx") == 0);

        //g->load_idx(idx_path);
        g->index();
        g->update_dec();
        g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-uidx+") == 0 || strcmp(argv[1],"-uidx+-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);

        auto *g = new Graph();
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-uidx+") == 0);

        g->index();
        g->update_inc();
        // g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-uidx-bl") == 0 || strcmp(argv[1],"-uidx-bl-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        float ratio = 1;
        if (argc > 5) {
            ratio = atof(argv[5]);
        }

        auto *g = new Graph();
        g->s_f_ = ratio;
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-uidx-bl") == 0);

        g->index();
        //g->load_idx(idx_path);
        g->update_baseline();
        // g->write_idx(idx_path);
        delete g;
    }

    if(strcmp(argv[1],"-rm") == 0 || strcmp(argv[1],"-rm-4col") == 0){
        string graph_path(argv[2]);
        string idx_path(argv[3]);
        string log_path(argv[4]);
        float ratio = 1;
        if (argc > 5) {
            ratio = atof(argv[5]);
        }

        auto *g = new Graph();
        g->init_log(log_path);

        g->load(graph_path,strcmp(argv[1],"-rm") == 0);

        g->load_idx(idx_path);
        g->remove_expired_t(g->t_max_*ratio);

        // g->write_idx(idx_path);
        delete g;
    }

    return 0;
}
