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
// Date: 2014/5/17                                               //
// Latent Collaborative Retrival with Pairwise Learning for      //
//   implicit feedback data.                                     //
// Factorization Format:                                         // 
//   Query*Item+User*Item                                        //
///////////////////////////////////////////////////////////////////

#pragma once

#include "../utils.hpp"

using namespace std;
using namespace __gnu_cxx;

struct TRAIN_PAIR{
    int uidx;
    int pidx;
    int qidx;
    int pidx1;
    int pidx2;
};
typedef struct TRAIN_PAIR TrainPair;

class LCR{
    ///////////////////////////////////////////////////////////////
    // Method control variable                                   //
    int niters;                                                  //
    int nsample;                                                 //
                                                                 //
    // Hyper-parameter setting                                   //
    int ndim;                                                    //
    double lr;                                                   //
    double u_reg;                                                //
    double p_reg;                                                //
    double q_reg;                                                //
    double bp_reg;                                               //
                                                                 //
    // Location related parameter                                //
    int* index_extent;                                           // 
    double grain_lat;                                            //
    double grain_lng;                                            //
    int ndimx;                                                   //
    int ndimy;                                                   //
    ///////////////////////////////////////////////////////////////

    Grid ** grids_pois;
    map<string, Poi*>* pois_latlng;

    map<string, int> user_ids;
    map<int, string> ruser_ids;
    map<string, int> poi_ids;
    map<int, string> rpoi_ids;
    map<string, int> query_ids;
    map<int, string> rquery_ids;
    map<string, string> poi_query;

    int n_users;
    int n_pois;
    int n_querys;
    double * total_factor;
    double ** user_factor;
    double ** poi_factor;
    double ** query_factor;
    double * poi_bias;
   
    int re_topk;
    bool bias_tag, restart_tag;
    char* trdata_path;
    char* vadata_path;
    char* tedata_path;
    char* query_path;
    char* user_factor_path;
    char* poi_factor_path;
    char* query_factor_path;
    char* poi_bias_path;

public:
    void paraInit(){
        // Method control variable 
        niters = 100;
        nsample = 5; 
    
        // Hyper-parameter setting 
        ndim = 20; 
        lr = 0.05; 
        u_reg = 0.1; 
        p_reg = 0.1; 
        q_reg = 0.1;
        bp_reg = 0.1; 
    
        // Location related parameter 
        index_extent = new int[4];
        index_extent[0] = -90, index_extent[1] = -180;
        index_extent[2] = 90, index_extent[3] = 180; 
        grain_lat = 0.05; 
        grain_lng = 0.05; 
        ndimx = int((index_extent[3]-index_extent[1])/grain_lng);
        ndimy = int((index_extent[2]-index_extent[0])/grain_lat);

        grids_pois = NULL;
        pois_latlng = NULL;

        n_users = 0;
        n_pois = 0;
        total_factor = NULL;
        user_factor = NULL;
        poi_factor = NULL;
        poi_bias = NULL;
   
        re_topk = -1;
        bias_tag = false;
        restart_tag = true;
        trdata_path = NULL;
        vadata_path = NULL;
        tedata_path = NULL;
        query_path = NULL;
        user_factor_path = NULL;
        poi_factor_path = NULL;
        poi_bias_path = NULL;
    }
    
    LCR(char* trdata_path, char* vadata_path, char* tedata_path,
            char* query_path, char* user_factor_path, char* poi_factor_path,
            char* query_factor_path, char* poi_bias_path, char* poi_path,
            char* grid_path, int data_num, bool bias_tag, bool restart_tag,
            int re_topk){
        paraInit();
        this->trdata_path = trdata_path;
        this->vadata_path = vadata_path;
        this->tedata_path = tedata_path;
        this->query_path  = query_path;
        this->user_factor_path = user_factor_path;
        this->poi_factor_path = poi_factor_path;
        this->query_factor_path = query_factor_path;
        this->poi_bias_path = poi_bias_path;
        this->bias_tag = bias_tag;
        this->restart_tag = restart_tag;
        this->re_topk = re_topk;

        printf("Model init.\n");
        utils::getIdMapRelation(trdata_path, &user_ids, &ruser_ids,
                &poi_ids, &rpoi_ids, &n_users, &n_pois);       
        printf("User: %d, Poi: %d\n", n_users, n_pois);
        utils::getPoiQueryInfo(query_path, &poi_query, &query_ids,
                &rquery_ids, n_querys);
        printf("Query: %d\n", n_querys);
        grids_pois = utils::loadGridInfo(grid_path, ndimx, ndimy);
        pois_latlng = utils::loadPoiInfo(&poi_ids, poi_path, data_num);

        if (restart_tag == true) {
            factorInit();
        } else {
            loadModel();
        }
        printf("Finishing initing.\n");
    }
    
    ~LCR() {
        if (!grids_pois) {
            for (int i=0; i<ndimx; i++)
                delete[] grids_pois[i];
            delete[] grids_pois;
        }
        
        if (!pois_latlng) {
            for (map<string, Poi*>::iterator it = pois_latlng->begin();
                    it != pois_latlng->end(); it++) {
                delete it->second;
            }
            delete pois_latlng;
        }
        
        if (total_factor)
            delete[] total_factor;
        if (user_factor)
            delete[] user_factor;
        if (poi_factor)
            delete[] poi_factor;

        user_ids.clear();
        map<string, int>(user_ids).swap(user_ids);
        ruser_ids.clear();
        map<int, string>(ruser_ids).swap(ruser_ids);
        poi_ids.clear();
        map<string, int>(poi_ids).swap(poi_ids);
        rpoi_ids.clear();
        map<int, string>(rpoi_ids).swap(rpoi_ids);
        query_ids.clear();
        map<string, int>(query_ids).swap(query_ids);
        rquery_ids.clear();
        map<int, string>(rquery_ids).swap(rquery_ids);
        poi_query.clear();
        map<string, string>(poi_query).swap(poi_query);
    }


    void factorInit() {
        int num_para=0;
        if (bias_tag)
            num_para = (n_users+n_pois+n_querys)*ndim+n_pois;
        else
            num_para = (n_users+n_pois+n_querys)*ndim;
        
        total_factor = new double[num_para];
        int ind = 0;
        user_factor = new double*[n_users];
        for (int u=0; u<n_users; u++) {
            user_factor[u] = total_factor+ind;
            utils::muldimGaussrand(user_factor[u], ndim);
            //utils::muldimUniform(user_factor[u], ndim);
            //uitls::muldimZero(user_factor[u], ndim);
            ind += ndim;
        }
        
        poi_factor = new double*[n_pois];
        for (int p=0; p<n_pois; p++) {
            poi_factor[p] = total_factor+ind;
            utils::muldimGaussrand(poi_factor[p], ndim);
            //utils::muldimUniform(poi_factor[p], ndim);
            //uitls::muldimZero(poi_factor[p], ndim);
            ind += ndim;
        }

        query_factor = new double*[n_querys];
        for (int q=0; q<n_querys; q++) {
            query_factor[q] = total_factor+ind;
            utils::muldimGaussrand(query_factor[q], ndim);
            //utils::muldimUniform(query_factor[q], ndim);
            //uitls::muldimZero(query_factor[q], ndim);
            ind += ndim;
        }
        
        if (bias_tag) {
            poi_bias = total_factor+ind;
            utils::muldimGaussrand(poi_bias, n_pois);
            //utils::muldimUniform(poi_bias, n_pois);
            //uitls::muldimZero(poi_bias, n_pois);
            ind += n_pois;
        }
    }

    vector<TrainPair*>* genTrainPairs() {
        string line;
        string uid, pid, pid1, pid2;
        vector<string> parts;
        vector<Coordinate*>* near_grids = NULL;
        vector<string>* candidate_pois = NULL;
        vector<string>* neg_samples = NULL;
        vector<TrainPair*>* tr_pairs = new vector<TrainPair*>();

        printf("Start generating training pairs.\n");
        //timeval t1, t2;
        int finished_num = 0;
        ifstream* in = utils::ifstream_(trdata_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid = parts[1];
            pid1 = parts[4];
            
            candidate_pois = new vector<string>();
            near_grids = utils::getNearGridsForPoi((*pois_latlng)[pid],
                                ndimx, ndimy, grain_lng, grain_lat, true);
            for (vector<Coordinate*>::iterator it2=near_grids->begin();
                    it2!=near_grids->end(); it2++)
                candidate_pois->insert(candidate_pois->end(),
                                grids_pois[(*it2)->x][(*it2)->y].pois.begin(),
                                grids_pois[(*it2)->x][(*it2)->y].pois.end());
            neg_samples = utils::genSamples(candidate_pois, pid1, nsample);
            
            for (vector<string>::iterator it2=neg_samples->begin();
                        it2!=neg_samples->end(); it2++) {
                pid2 = *it2;
                TrainPair * pair = new TrainPair();
                pair->uidx = user_ids[uid];
                pair->pidx = poi_ids[pid];
                pair->qidx = query_ids[poi_query[pid]];
                pair->pidx1 = poi_ids[pid1];
                pair->pidx2 = poi_ids[pid2];
                tr_pairs->push_back(pair);
                printf("\rGenerated training pairs: %d", finished_num);
                finished_num++;
            }
            delete candidate_pois;
            delete neg_samples;
            delete near_grids;
        }
        printf("\nGenerated Training Pair Num: %ld\n", tr_pairs->size());
        return tr_pairs;
    }


    void train() {
        double p_val=0.0, n_val=0.0, logit_loss=0.0;
        int uidx, pidx, qidx, pidx1, pidx2, finished_num;
        double * guser_factor = new double[ndim];
        double * gquery_factor = new double[ndim];
        double * gpoi_factor1 = new double[ndim];
        double * gpoi_factor2 = new double[ndim];
        double * t_diff = new double[ndim];
        double * t_sum = new double[ndim];
        
        vector<TrainPair*>* tr_pairs = genTrainPairs();
        for (int i=0; i<niters; i++) {
            finished_num = 0;
            random_shuffle(tr_pairs->begin(), tr_pairs->end());
            for (vector<TrainPair*>::iterator it=tr_pairs->begin();
                    it!=tr_pairs->end(); it++) {
                uidx = (*it)->uidx;
                pidx = (*it)->pidx;
                qidx = (*it)->qidx;
                pidx1 = (*it)->pidx1;
                pidx2 = (*it)->pidx2;
                if (bias_tag) {
                    p_val = utils::dot(user_factor[uidx], poi_factor[pidx1], ndim)
                          + utils::dot(query_factor[qidx], poi_factor[pidx1], ndim)
                          + poi_bias[(*it)->pidx1];
                    n_val = utils::dot(user_factor[uidx], poi_factor[pidx2], ndim)
                          + utils::dot(query_factor[qidx], poi_factor[pidx2], ndim)
                          + poi_bias[(*it)->pidx2];
                } else {
                    p_val = utils::dot(user_factor[uidx], poi_factor[pidx1], ndim)
                          + utils::dot(query_factor[qidx], poi_factor[pidx1], ndim);
                    n_val = utils::dot(user_factor[uidx], poi_factor[pidx2], ndim)
                          + utils::dot(query_factor[qidx], poi_factor[pidx2], ndim);
                }
                logit_loss = utils::logitLoss(p_val-n_val);
                
                // compute factor gradient
                #pragma omp parallel for
                for (int j=0; j<ndim; j++) {
                    t_diff[j] = poi_factor[pidx1][j]-poi_factor[pidx2][j];
                    guser_factor[j] = logit_loss*t_diff[j]
                                      -u_reg*user_factor[uidx][j];
                    gquery_factor[j] = logit_loss*t_diff[j]
                                       -p_reg*query_factor[qidx][j];
                    t_sum[j] = user_factor[uidx][j]+query_factor[qidx][j];
                    gpoi_factor1[j] = logit_loss*t_sum[j]
                                      -p_reg*poi_factor[pidx1][j];
                    gpoi_factor2[j] = -logit_loss*t_sum[j]
                                      -p_reg*poi_factor[pidx2][j];
                }
                // update factor
                #pragma omp parallel for
                for (int j=0; j<ndim; j++) {
                    user_factor[uidx][j] = user_factor[uidx][j] +lr*guser_factor[j];
                    query_factor[qidx][j] = query_factor[qidx][j] +lr*gquery_factor[j];
                    poi_factor[pidx1][j] = poi_factor[pidx1][j] +lr*gpoi_factor1[j];
                    poi_factor[pidx2][j] = poi_factor[pidx2][j] +lr*gpoi_factor2[j];
                }
                // compute and update poi bias if necessary
                if (bias_tag) {
                    poi_bias[pidx1] = poi_bias[pidx1]
                                    + lr*(logit_loss-bp_reg*poi_bias[pidx1]);
                    poi_bias[pidx2] = poi_bias[pidx2]
                                    + lr*(-logit_loss-bp_reg*poi_bias[pidx2]);
                }
                
                finished_num++;
                if (finished_num%100000 == 0) {
                    printf("\rCurrent Iteration: %d, Finished Training Pair Num: %d.", i+1, finished_num);
                    fflush(stdout);
                }
            }
            printf("Current Iteration: %d, AUC is %.6f...\n", i+1, evaluation(tr_pairs));
        }
        saveModel();
        
        //release memory
        for (vector<TrainPair*>::iterator it=tr_pairs->begin();
            it!=tr_pairs->end(); it++){
            delete *it;
            *it = NULL;
        }
        delete t_diff;
        delete t_sum;
        delete guser_factor;
        delete gquery_factor;
        delete gpoi_factor1;
        delete gpoi_factor2;
        delete tr_pairs; 
    }

    double evaluation(vector<TrainPair*>* tr_pairs) {
        int correct_num = 0;
        int uidx, pidx, qidx, pidx1, pidx2;
        double p_val, n_val;

        for (vector<TrainPair*>::iterator it=tr_pairs->begin();
                it!=tr_pairs->end(); it++) {
            uidx = (*it)->uidx;
            pidx = (*it)->pidx;
            qidx = (*it)->qidx;
            pidx1 = (*it)->pidx1;
            pidx2 = (*it)->pidx2;
            if (bias_tag) {
                p_val = utils::dot(user_factor[uidx], poi_factor[pidx1], ndim)
                      + utils::dot(query_factor[qidx], poi_factor[pidx1], ndim)
                      + poi_bias[(*it)->pidx1];
                n_val = utils::dot(user_factor[uidx], poi_factor[pidx2], ndim)
                      + utils::dot(query_factor[qidx], poi_factor[pidx2], ndim)
                      + poi_bias[(*it)->pidx2];
            } else {
                p_val = utils::dot(user_factor[uidx], poi_factor[pidx1], ndim)
                      + utils::dot(query_factor[qidx], poi_factor[pidx1], ndim);
                n_val = utils::dot(user_factor[uidx], poi_factor[pidx2], ndim)
                      + utils::dot(query_factor[qidx], poi_factor[pidx2], ndim);
            }
            if (p_val > n_val)
                correct_num += 1;
        }
        return 1.0*correct_num/(tr_pairs->size());
    }

    void recommendation(char* submission_path) {
        int re_num = 0, finished_num = 0;
        double score;
        string uid, pid1, c_pid, line;
        vector<string> parts;
        vector<Coordinate*>* near_grids = NULL;
        vector<string>* candidate_pois = NULL;
        vector<string>* neg_samples = NULL;
        vector<vector<string> >* recommendation_result = NULL;
        vector<string>* onere_result = NULL;
        vector<Rateval>* inter_result = NULL;
        Rateval* rateval = NULL;

        printf("\nStart recommendation!\n");
        recommendation_result = new vector<vector<string> >();
        ifstream* in = utils::ifstream_(tedata_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid1 = parts[1];
            near_grids = utils::getNearGridsForPoi((*pois_latlng)[pid1],
                                ndimx, ndimy, grain_lng, grain_lat, true);
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
                score = utils::dot(user_factor[user_ids[uid]], 
                            poi_factor[poi_ids[c_pid]], ndim)
                      + utils::dot(query_factor[query_ids[poi_query[pid1]]],
                            poi_factor[poi_ids[c_pid]], ndim);
                if (bias_tag)
                    score += poi_bias[poi_ids[c_pid]];
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
            finished_num++;
            recommendation_result->push_back(*onere_result);
            printf("\rFinished Recommendation Pairs: %d", finished_num);
            fflush(stdout);

            // release memory
            delete inter_result;
            delete candidate_pois;
            delete near_grids;
            delete neg_samples;
        }
        utils::write_submission(recommendation_result, submission_path);
        
        // release memory
        //for (vector<vector<string> >::iterator it=recommendation_result->begin();
        //        it!=recommendation_result->end(); it++) {
        //    it->clear();
        //    vector<string>(*it).swap(*it);
        //}
        delete recommendation_result;
    }

    void recommendationNewPOI(char* submission_path) {
        int re_num = 0, finished_num = 0;
        double score;
        string uid, pid1, pid2, c_pid, line;
        vector<string> parts;
        vector<Coordinate*>* near_grids = NULL;
        vector<string>* candidate_pois = NULL;
        vector<string>* neg_samples = NULL;
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
            near_grids = utils::getNearGridsForPoi((*pois_latlng)[pid1],
                                ndimx, ndimy, grain_lng, grain_lat, true);
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
                score = utils::dot(user_factor[user_ids[uid]], 
                            poi_factor[poi_ids[c_pid]], ndim)
                      + utils::dot(query_factor[query_ids[poi_query[pid1]]],
                            poi_factor[poi_ids[c_pid]], ndim);
                if (bias_tag)
                    score += poi_bias[poi_ids[c_pid]];
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
            delete candidate_pois;
            delete near_grids;
            delete neg_samples;
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

    void saveModel() {
        FILE* f = utils::fopen_(user_factor_path, "w");
        for (int u=0; u<n_users; u++)
            fwrite(user_factor[u], sizeof(double), ndim, f);
        fclose(f);

        f = utils::fopen_(poi_factor_path, "w");
        for (int p=0; p<n_pois; p++)
            fwrite(poi_factor[p], sizeof(double), ndim, f);
        fclose(f);
        
        f = utils::fopen_(query_factor_path, "w");
        for (int q=0; q<n_querys; q++)
            fwrite(query_factor[q], sizeof(double), ndim, f);
        fclose(f);

        if (bias_tag) {
            f = utils::fopen_(poi_bias_path, "w");
            fwrite(poi_bias, sizeof(double), n_pois, f);
            fclose(f);
        }
    }

    void loadModel() {
        int num_para=0;
        if (bias_tag)
            num_para = (n_users+n_pois+n_querys)*ndim+n_pois;
        else
            num_para = (n_users+n_pois+n_querys)*ndim;
        
        double * factor = new double[num_para];
        int ind = 0;
        user_factor = new double*[n_users];
        for (int u=0; u<n_users; u++) {
            user_factor[u] = factor+ind;
            ind += ndim;
        }
        poi_factor = new double*[n_pois];
        for (int p=0; p<n_pois; p++) {
            poi_factor[p] = factor+ind;
            ind += ndim;
        }
        query_factor = new double*[n_querys];
        for (int q=0; q<n_querys; q++) {
            query_factor[q] = total_factor+ind;
            ind += ndim;
        }
        poi_bias = new double[n_pois];
        if (bias_tag) {
            poi_bias = factor+ind;
            ind += n_pois;
        }
        
        FILE* f = utils::fopen_(user_factor_path, "r");
        for (int u=0; u<n_users; u++)
            utils::fread_(user_factor[u], sizeof(double), ndim, f);
        fclose(f);
        f = utils::fopen_(poi_factor_path, "r");
        for (int p=0; p<n_pois; p++)
            utils::fread_(poi_factor[p], sizeof(double), ndim, f);
        fclose(f);
        f = utils::fopen_(query_factor_path, "r");
        for (int q=0; q<n_querys; q++)
            utils::fread_(query_factor[q], sizeof(double), ndim, f);
        fclose(f);
        if (bias_tag) {
            f = utils::fopen_(poi_bias_path, "r");
            utils::fread_(poi_bias, sizeof(double), n_pois, f);
            fclose(f);
        }
    }
};
