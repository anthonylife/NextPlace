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
// Date: 2014/5/11                                               //
// Probabilistic Matric Factorization with Pairwise Learning for //
//   implicit feedback data.                                     //
///////////////////////////////////////////////////////////////////

#pragma once

#include "../utils.hpp"

using namespace std;

struct TRAIN_PAIR{
    int uidx;
    int pidx1;
    int pidx2;
};
typedef struct TRAIN_PAIR TrainPair;

class PMF{
    ///////////////////////////////////////////////////////////////
    // Method control variable                                   //
    int niters = 25;                                             //
    int nsample = 5;                                             //
                                                                 //
    // Hyper-parameter setting                                   //
    int ndim = 10;                                               //
    double lr = 0.1;                                             //
    double u_reg = 0.1;                                          //
    double p_reg = 0.1;                                          //
    double bp_reg = 0.1;                                         //
                                                                 //
    // Location related parameter                                //
    int[] index_extent = {-90, -180, 90, 180};                   // 
    double grain_lat = 0.05;                                     //
    double grain_lng = 0.05;                                     //
    int ndimx = int((index_extent[3]-index_extent[1])/grain_lng);//
    int ndimy = int((index_extent[2]-index_extent[0])/grain_lat);//
    ///////////////////////////////////////////////////////////////

    Grid ** grids_pois = NULL;
    map<string, Poi*>* pois_latlng = NULL;

    map<string, int> user_ids;
    map<int, string> ruser_ids;
    map<string, int> poi_ids;
    map<int, string> rpoi_ids;

    int n_users;
    int n_pois;
    double ** user_factor = NULL;
    double ** poi_factor = NULL;
    double * poi_bias = NULL;
    
    bool bias_tag;
    bool restart_tag;
    char* trdata_path;
    char* vadata_path;
    char* tedata_path;

public:
    PMF(char* trdata_path, char* vadata_path, char* tedata_path,
            char* poi_path, char* grid_path, int data_num, bool bias_tag, bool restart_tag){
        this->trdata_path = trdata_path;
        this->vadata_path = vadata_path;
        this->tedata_path = tedata_path;
        this->bias_tag = bias_tag;
        this->restart_tag = restart_tag;
        
        utils::getIdMapRelation(trdata_path, &user_ids, &ruser_ids, &poi_ids, &rpoi_ids, &n_users, &n_pois);       
        grids_pois = utils::loadGridInfo(grid_path, ndimx, ndimy);
        pois_latlng = utils::loadPoiInfo(poi_ids, poi_path, data_num);

        if (restart_tag == true) {
            factor_init();
        } else {
            loadModel();
        }
    
    }

    void factor_init() {
        int num_para;
        if (bias_tag)
            num_para = (n_users+n_pois)*ndim+n_pois;
        else
            num_para = (n_users+n_pois)*ndim;
        
        double * factor = new double[num_para];
        int idx = 0;
        for (int u=0; u<n_users; u++) {
            user_factor[u] = factor+ind;
            utils::muldimGaussrand(&user_factor[u], ndim);
            //utils::muldimUniform(&user_factor[u], ndim);
            //uitls::muldimZero(&user_factor[u], ndim);
            ind += ndim;
        }
        for (int p=0; p<n_pois; p++) {
            poi_factor[p] = factor+ind;
            utils::muldimGaussrand(&poi_factor[p], ndim);
            //utils::muldimUniform(&poi_factor[p], ndim);
            //uitls::muldimZero(&poi_factor[p], ndim);
            ind += ndim;
        }
        if (bias_tag) {
            poi_bias = factor+ind;
            utils::muldimGaussrand(&poi_bias, n_pois);
            //utils::muldimUniform(&poi_bias, n_pois);
            //uitls::muldimZero(&poi_bias, n_pois);
            ind += n_pois;
        }
    }

    vector<TrainPair*>* genTrainPairs() {
        string uid, pid1, pid2;
        vector<string> parts;
        vector<TrainPair*>* tr_pairs = new vector<TrainPair*>();
        map<string, set<string>*> *user_visit = new map<string, set<string>*>();
        vector<Coordinate*>* near_grids = NULL;
        vector<string>* candidate_pois = NULL;
        vector<string>* neg_samples = NULL;

        for (map<string, int>::iterator it=user_ids.begin(); it!=user_ids.end(); it++) {
            set<string>* tmp_set = new set<string>();
            (*user_visit)[it->first] = tmp_set;
        }

        ifstream* in = utils::ifstream_(trdata_path);
        while (getline(*in, line)) {
            parts = utils::split_str(line, ',');
            uid = parts[0];
            pid1 = parts[1];
            pid2 = parts[2];
            if ((*user_visit)[uid]->find(pid1) == (*user_visit)[uid]->end()) {
                (*user_visit)[uid]->insert(pid1);
            }
            if ((*user_visit)[uid]->find(pid2) == (*user_visit)[uid]->end()) {
                (*user_visit)[uid]->insert(pid2);
            }
        }
               
        for (map<string, set<string>*>::iterator it=user_visit->begin();
                it!=user_visit->end(); it++) {
            uid = it->first;
            for (set<string>::iterator it1=it->second->begin();
                    it1!=it->second->end(); it1++) {
                pid1 = *it1;
                candidate_pois = new vector<string>();
                near_grids = getNearGridsForPoi(pois_latlng[pid1],
                                ndimx, ndimy, grain_lng, grain_lat, true);
                for (vector<Coordinate*>::iterator it2=near_grids->begin();
                        it2!=near_grids->end(); it2++)
                    candidate_pois->insert(candidate_pois.end(),
                                    grids_pois[it2->x][it2->y].pois.begin(),
                                    grids_pois[it2->x][it2->y].pois.end());
                
                neg_samples = utils::genNegSamples(candidate_pois, it-second, nsample);
                for (vector<string>::iterator it2=neg_samples->begin();
                        it2!=neg_samples->end(); it2++) {
                    pid2 = *it2;
                    TrainPair * pair = new TrainPair();
                    pair->uidx = user_ids[uid];
                    pair->pidx1 = poi_ids[pid1];
                    pair->pidx2 = poi_ids[pid2];
                    tr_pairs->push_back(pair);
                }
                delete candidate_pois;
                delete neg_samples;
            }
        }

        if (user_visit != NULL)
            for (map<string, set*>::iterator it=user_visit.begin(); it!=user_ids.end(); it++)
                delete it->second;
            delete user_visit;
    
        return tr_pairs;
    }


    void train() {
        vector<TrainPair*>* tr_pairs = genTrainPairs();
        
        for (int i=0; i<niters; i++) {
            random_shuffle(tr_pairs->begin(), tr_pairs->end());
            for ()
        }

        for (vector<TrainPair*>::iterator it=tr_pairs.begin(); it!=tr_pairs.end(); it++)
            delete it;
        delete tr_pairs; 
    }


    ~PMF(){
        if (!grids ) {
            for (int i=0; i<ndimx; i++)
                delete[] grids[i];
            delete[] grids;
        }

        if (!pois_latlng) {
            for (map<string, Poi*>::iterator it = pois_latlng->begin(); it != pois_latlng->end; it++)
                delete it->second;
            delete pois_latlng;
        }

        user_ids.clear();
        map<string, int>(user_ids).swap(user_ids);
        ruser_ids.clear();
        map<int, string>(ruser_ids).swap(ruser_ids);
        pois_ids.clear();
        map<string, int>(poi_ids).swap(poi_ids);
        rpoi_ids.clear();
        map<int, string>(rpoi_ids).swap(rpoi_ids);
    }
};
