//Copyright [2014] [Wei Zhang]

//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//http://www.apache.org/licenses/LICENSE-2.0
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

///////////////////////////////////////////////////////////////////
// Date: 2015/4/9                                                //
// Transition Geo-Temporal Topic Model for Successive Location   //
//   Recommendation.                                             //
///////////////////////////////////////////////////////////////////


#pragma once

#include "../utils.hpp"

using namespace std;
using namespace __gnu_cxx;

struct VISIT_PAIR{
    int uidx;
    int pidx;
    int pidx1;
    int timeidx;
};
typedef struct VISIT_PAIR VisitPair;


class TGTTM{
    // Specified parameters
    const static int sampling_iters = 100;
    const static int T = 40;
    const static double gamma_t = 0.1;
    const static double gamma_ut = 0.1;
    const static double gamma_lt = 0.1;
    const static double gamma_tt = 0.1;
    const static double gamma_tl = 0.1;

    const static int factor_type = 2;
    const static double size_lat=0.05;
    const static double size_lng=0.05;
    int index_extent[4];
    int ndimx;
    int ndimy;
    
    // Storing variables
    Grid ** grids_pois;
    map<string, Poi*>* pois_latlng;
    vector<VisitPair*>* tr_pairs;
    vector<VisitPair*>* te_pairs;
    vector<VisitPair*>* va_pairs;
    
    map<string, int> user_ids;
    map<int, string> ruser_ids;
    map<string, int> poi_ids;
    map<int, string> rpoi_ids;

    int n_users;
    int n_pois;
    int n_docs;

    double alpha;
    int * pair_y;
    int * sum_pair_y;
    int * pair_loc_topic;
    int * pair_trans_topic;
    int * pair_user_topic;
    int ** topic_topic;
    int * sum_topic_topic;
    int ** location_topic;
    int * sum_location_topic;
    int ** user_topic;
    int * sum_user_topic;
    int ** topic_location;
    int * sum_topic_location;
    double ** theta_user;
    double ** theta_pois;
    double ** theta_poie;
    double ** theta_tt;

    int re_topk;
    bool restart_tag;
    char* trdata_path;
    char* vadata_path;
    char* tedata_path;

public:
    TGTTM(char *trdata_path, char *vadata_path, char *tedata_path,
            char *model_path, char * poi_path, char * grid_path, int data_num,
            bool restart_tag, int re_topk) {
        this->trdata_path = trdata_path;
        this->vadata_path = vadata_path;
        this->tedata_path = tedata_path;
        index_extent[0] = -90, index_extent[1] = -180, index_extent[2] = 90, index_extent[3] = 180;
        ndimx = int((index_extent[3]-index_extent[1])/size_lng);
        ndimy = int((index_extent[2]-index_extent[0])/size_lat);
        printf("Model init.\n");
        utils::getIdMapRelation(trdata_path, &user_ids, &ruser_ids,
                &poi_ids, &rpoi_ids, &n_users, &n_pois);       
        printf("User: %ld, Poi: %ld\n", user_ids.size(), poi_ids.size());
        grids_pois = utils::loadGridInfo(grid_path, ndimx, ndimy);
        pois_latlng = utils::loadPoiInfo(&poi_ids, poi_path, data_num);
        tr_pairs = genVisitPairs(trdata_path);
        n_docs = tr_pairs->size();
        printf("Number of pairs: %d\n", n_docs);
        te_pairs = genVisitPairs(tedata_path);
        va_pairs = genVisitPairs(vadata_path);
        
        if (restart_tag == true) {
            modelInit();
        } else {
            loadModel();
        }
        printf("Finishing initing.\n");
        
    }
   
    ~TGTTM() {
        user_ids.clear();
        map<string, int>(user_ids).swap(user_ids);
        ruser_ids.clear();
        map<int, string>(ruser_ids).swap(ruser_ids);
        poi_ids.clear();
        map<string, int>(poi_ids).swap(poi_ids);
        rpoi_ids.clear();
        map<int, string>(rpoi_ids).swap(rpoi_ids);
        
        delete tr_pairs;
        delete te_pairs;
        delete va_pairs;
    }

    vector<VisitPair*>* genVisitPairs(char * data_path) {
        string line;
        int timeidx;
        string uid, pid, pid1, pid2, day, hour;
        vector<string> parts;
        vector<VisitPair*>* visit_pairs = new vector<VisitPair*>; 

        printf("Start generating visiting pairs.\n");
        //timeval t1, t2;
        int finished_num = 0;
        ifstream* in = utils::ifstream_(data_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid = parts[1];
            day = parts[2];
            hour = parts[3];
            pid1 = parts[4];
            
            timeidx = utils::getTimeId(atoi(day.c_str()), atoi(hour.c_str()), 6);
            VisitPair * pair = new VisitPair();
            pair->uidx = user_ids[uid];
            pair->pidx = poi_ids[pid];
            pair->pidx1 = poi_ids[pid1];
            pair->timeidx = timeidx;
            visit_pairs->push_back(pair);
            finished_num++;
            
            if (finished_num % 100000 == 0)
                printf("\rGenerated training pairs: %d", finished_num);

            //if (finished_num > 88962*0.2)
            //    break;
            //if (finished_num > 798280*0.8)
            //    break;
        }
        printf("\nGenerated Visiting Pair Num: %ld\n", visit_pairs->size());
        return visit_pairs;
    }
    
    void modelInit() {
        int uidx, pidx, pidx1;
        int y, z_ut, z_lt, z_tt;

        // Factor init
        theta_user = new double*[n_users];
        for (int u=0; u<n_users; u++)
            theta_user[u] = new double[n_users];
        theta_pois = new double*[n_pois];
        theta_poie = new double*[n_pois];
        for (int i=0; i<n_pois; i++) {
            theta_pois[i] = new double[T];
            theta_poie[i] = new double[T];
        }
        theta_tt = new double*[T];
        for (int t=0; t<T; t++) {
            theta_tt[t] = new double[T];
        }

        pair_y = new int[n_docs];
        sum_pair_y = new int[factor_type];
        memset(sum_pair_y, 0, sizeof(int)*factor_type);
        pair_loc_topic = new int[n_docs];
        pair_trans_topic = new int[n_docs];
        pair_user_topic = new int[n_docs];
        topic_topic = new int*[T];
        for (int t=0; t<T; t++) {
            topic_topic[t] = new int[T];
            memset(topic_topic[t], 0, sizeof(T));
        }
        sum_topic_topic = new int[T];
        memset(sum_topic_topic, 0, sizeof(int)*T);
        location_topic = new int*[n_pois];
        for (int i=0; i<n_pois; i++) {
            location_topic[i] = new int[T];
            memset(location_topic[i], 0, sizeof(int)*T);
        }
        sum_location_topic = new int[n_pois];
        memset(sum_location_topic, 0, sizeof(int)*n_pois);
        user_topic = new int*[n_users];
        for (int u=0; u<n_users; u++) {
            user_topic[u] = new int[T];
            memset(user_topic[u], 0, sizeof(int)*T);
        }
        sum_user_topic = new int[n_users];
        memset(sum_user_topic, 0, sizeof(int)*n_users);
        topic_location = new int*[T];
        for (int t=0; t<T; t++) {
            topic_location[t] = new int[n_pois];
            memset(topic_location[t], 0, sizeof(int)*n_pois);
        }
        sum_topic_location = new int[T];
        memset(sum_topic_location, 0, sizeof(T));

        // Sampling init
        int idx = 0;
        for (vector<VisitPair*>::iterator it=tr_pairs->begin();
                it!=tr_pairs->end(); it++) {
            uidx = (*it)->uidx;
            pidx = (*it)->pidx;
            pidx1 = (*it)->pidx1;
            y = (int)(((double)random()/RAND_MAX)*factor_type);
            z_ut = (int)(((double)random()/RAND_MAX)*T);
            z_lt = (int)(((double)random()/RAND_MAX)*T);
            z_tt = (int)(((double)random()/RAND_MAX)*T);

            pair_y[idx] = y;
            sum_pair_y[y]++;
            pair_loc_topic[idx] = z_lt;
            pair_trans_topic[idx] = z_tt;
            pair_user_topic[idx] = z_ut;
        
            if (y == 0) {
                topic_topic[z_lt][z_tt]++;
                sum_topic_topic[z_lt]++;
                location_topic[pidx][z_lt]++;
                sum_location_topic[pidx]++;
                topic_location[z_tt][pidx1]++;
                sum_topic_location[z_tt]++;
            } else if (y == 1) {
                user_topic[uidx][z_ut]++;
                sum_user_topic[uidx]++;
                topic_location[z_ut][pidx1]++;
                sum_topic_location[z_ut]++;
            }
            idx++;
        }
    }

    void train() {
        int uidx, pidx, pidx1, idx, sidx;
        int y, z_ut, z_lt, z_tt;
        double * P = new double[T+T*T], sval;
        double gamma_T = gamma_t*factor_type;
        double gamma_UT = gamma_ut*T;
        double gamma_LT = gamma_lt*T;
        double gamma_TT = gamma_tt*T;
        double gamma_TL = gamma_tl*n_pois;

        timeval start_t, end_t;
        for (int i=0; i<sampling_iters; i++) {
            utils::tic(start_t);
            idx = 0;
            for (vector<VisitPair*>::iterator it=tr_pairs->begin();
                    it!=tr_pairs->end(); it++) {
                uidx = (*it)->uidx;
                pidx = (*it)->pidx;
                pidx1 = (*it)->pidx1;
                memset(P, 0, sizeof(double)*(T+T*T));
                
                y = pair_y[idx];
                sum_pair_y[y]--;
                z_lt = pair_loc_topic[idx];
                z_tt = pair_trans_topic[idx];
                z_ut = pair_user_topic[idx];
                if (y == 0) {
                    topic_topic[z_lt][z_tt]--;
                    sum_topic_topic[z_lt]--;
                    location_topic[pidx][z_lt]--;
                    sum_location_topic[pidx]--;
                    topic_location[z_tt][pidx1]--;
                    sum_topic_location[z_tt]--;
                } else if (y == 1) {
                    user_topic[uidx][z_ut]--;
                    sum_user_topic[uidx]--;
                    topic_location[z_ut][pidx1]--;
                    sum_topic_location[z_ut]--;
                }

                // from user topic
                P[0] = (sum_pair_y[0]+gamma_t)/(n_docs+gamma_T)*(user_topic[uidx][0]+gamma_ut)/(sum_user_topic[uidx]+gamma_UT)*(topic_location[0][pidx1]+gamma_tl)/(sum_topic_location[0]+gamma_TL);
                for (int t=1; t<T; t++) {
                    P[t] = P[t-1] + (sum_pair_y[0]+gamma_t)/(n_docs+gamma_T)*(user_topic[uidx][t]+gamma_ut)/(sum_user_topic[uidx]+gamma_UT)*(topic_location[t][pidx1]+gamma_tl)/(sum_topic_location[t]+gamma_TL);
                }
                for (int t=0; t<T; t++)
                    for (int t1=0; t1<T; t1++)
                        P[t*T+t1+T] = P[t*T+t1+T-1] + (sum_pair_y[1]+gamma_t)/(n_docs+gamma_T)*(location_topic[pidx][t]+gamma_lt)/(sum_location_topic[pidx]+gamma_LT)*(topic_topic[t][t1]+gamma_tt)/(sum_topic_topic[t]+gamma_TT)*(topic_location[t1][pidx1]+gamma_tl)/(sum_topic_location[t1]+gamma_TL);
                
                sval = (double)random()/RAND_MAX*P[T*T+T-1];
                for (sidx=0; sidx<T*T+T-1; sidx++) {
                    if (sval < P[sidx])
                        break;
                }
                if (sidx < T) {
                    pair_y[sidx]=0;
                    sum_pair_y[0]++;
                    user_topic[uidx][sidx]++;
                    sum_user_topic[uidx]++;
                    topic_location[sidx][pidx1]++;
                    sum_topic_location[sidx]++;
                    
                } else {
                    pair_y[sidx]=1;
                    sum_pair_y[1]++;
                    sidx = sidx-T;
                    z_lt = (int)sidx/T;
                    z_tt = (int)sidx%T;
                    topic_topic[z_lt][z_tt]++;
                    sum_topic_topic[z_lt]++;
                    location_topic[pidx][z_lt]++;
                    sum_location_topic[pidx]++;
                    topic_location[z_tt][pidx1]++;
                    sum_topic_location[z_tt]++;
                }
                idx++;
                if ((idx+1) % 10000 == 0) {
                    printf("\rFinished scanning pairs: %d", idx);
                    utils::toc(start_t, end_t, true);
                }
            }
            printf("Finished Iteration: %d\n", i+1);
        }
        computeTheta();
    }
   
    void computeTheta() {
        double gamma_UT = gamma_ut*T;
        double gamma_LT = gamma_lt*T;
        double gamma_TT = gamma_tt*T;
        double gamma_TL = gamma_lt*n_pois;
        for (int u=0; u<n_users; u++)
            for (int t=0; t<T; t++)
                theta_user[u][t] = (user_topic[u][t]+gamma_ut)/(sum_user_topic[u]+gamma_UT);
        for (int i=0; i<n_pois; i++)
            for (int t=0; t<T; t++)
                theta_pois[i][t] = (location_topic[i][t]+gamma_lt)/(sum_location_topic[i]+gamma_LT);
        for (int t=0; t<T; t++)
            for (int t1=0; t1<T; t1++)
                theta_tt[t][t1] = (topic_topic[t][t1]+gamma_tt)/(sum_topic_topic[t]+gamma_TT);
        for (int t=0; t<T; t++)
            for (int i=0; i<n_pois; i++)
                theta_poie[t][i] = (topic_location[t][i]+gamma_tl)/(sum_topic_location[t]+gamma_TL);
        
        alpha = (sum_pair_y[0]+gamma_t)/(n_docs+factor_type*gamma_t);
    }

    double predRating(int uid, int cpid, int npid, int tid) {
        double rating1 = 0.0;
        double rating2 = 0.0;
        for (int t=0; t<T; t++)
            rating1 += theta_user[uid][t]*theta_poie[npid][t];
        for (int t=0; t<T; t++)
            for (int t1=0; t1<T; t1++)
                rating2 += theta_pois[cpid][t]*theta_tt[t][t1]*theta_poie[npid][t1];
        return alpha*rating1 + (1-alpha)*rating2;
    }

    void recommendationNewPOI(char* submission_path) {
        int re_num = 0, finished_num = 0, timeidx;
        double score;
        string uid, pid1, pid2, c_pid, line, hour, day;
        vector<string> parts;
        vector<Coordinate*>* near_grids = NULL;
        vector<string>* candidate_pois = NULL;
        vector<string>* onere_result = NULL;
        vector<Rateval>* inter_result = NULL;
        ifstream* in = NULL;
        Rateval* rateval = NULL; 
        hash_set<string>* tmp_set = NULL;
        map<string, hash_set<string>*>* user_visit = NULL;
        vector<vector<string> >* recommendation_result = NULL;

        printf("Start recommendation!\n");
        user_visit = new map<string, hash_set<string>*>();
        for (map<string, int>::iterator it=user_ids.begin();
                it!=user_ids.end(); it++) {
            tmp_set = new hash_set<string>();
            (*user_visit)[it->first] = tmp_set;
        }
        
        in = utils::ifstream_(trdata_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid1 = parts[1];
            pid2 = parts[4];
            if ((*user_visit)[uid]->find(pid1) == (*user_visit)[uid]->end()) {
                (*user_visit)[uid]->insert(pid1);
            }
            if ((*user_visit)[uid]->find(pid2) == (*user_visit)[uid]->end()) {
                (*user_visit)[uid]->insert(pid2);
            }
        }
        
        recommendation_result = new vector<vector<string> >();
        in = utils::ifstream_(tedata_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid1 = parts[1];
            day = parts[2];
            hour = parts[3];
            
            timeidx = utils::getTimeId(atoi(day.c_str()), atoi(hour.c_str()), 6);
            near_grids = utils::getNearGridsForPoi((*pois_latlng)[pid1],
                                ndimx, ndimy, size_lng, size_lat, true);
            candidate_pois = new vector<string>();
            for (vector<Coordinate*>::iterator it2=near_grids->begin();
                        it2!=near_grids->end(); it2++)
                candidate_pois->insert(candidate_pois->end(),
                                grids_pois[(*it2)->x][(*it2)->y].pois.begin(),
                                grids_pois[(*it2)->x][(*it2)->y].pois.end());
            inter_result = new vector<Rateval>();
            for (vector<string>::iterator it=candidate_pois->begin();
                    it!=candidate_pois->end(); it++) {
                c_pid = *it;
                if ((*user_visit)[uid]->find(uid) != (*user_visit)[uid]->end())
                    continue;
                score = predRating(user_ids[uid], poi_ids[pid1], poi_ids[c_pid], timeidx);
                //score = score*geo_a*pow(utils::directDistance((*pois_latlng)[pid1]->lat, (*pois_latlng)[pid1]->lng, (*pois_latlng)[c_pid]->lat, (*pois_latlng)[c_pid]->lng), geo_b);
                rateval = new Rateval();
                rateval->id = c_pid;
                rateval->score = score;
                inter_result->push_back(*rateval);
                delete rateval;
            }
            re_num = 0;
            sort(inter_result->begin(), inter_result->end(), utils::greaterCmp);
            onere_result = new vector<string>();
            for (vector<Rateval>::iterator it=inter_result->begin();
                    it!=inter_result->end(); it++) {
                onere_result->push_back(it->id);
                re_num += 1;
                if (re_num == re_topk)
                    break;
            }
            recommendation_result->push_back(*onere_result);
            printf("\rFinished Recommendation Pairs: %d", finished_num);
            finished_num++;

            // release memory
            delete inter_result;
            for (vector<Coordinate*>::iterator it1=near_grids->begin();
                    it1!=near_grids->end(); it1++)
                delete (*it1);
            near_grids->clear();
            vector<Coordinate*>(*near_grids).swap(*near_grids);
            delete near_grids;
            candidate_pois->clear();
            vector<string>(*candidate_pois).swap(*candidate_pois);
            delete candidate_pois;
            delete onere_result;
            parts.clear();
            vector<string>(parts).swap(parts);
        }
        printf("\nRecommendation Finished!\n");
        utils::write_submission(recommendation_result, submission_path);
        
        // release memory
        if (user_visit != NULL) {
            for (map<string, hash_set<string>*>::iterator it=user_visit->begin();
                    it!=user_visit->end(); it++)
                delete it->second;
            delete user_visit;
        }
        delete recommendation_result;
    }

    void loadModel() {
    
    }

    void saveModel() {
    
    }

};
