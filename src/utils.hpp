#pragma once

#include<iostream>
#include<algorithm>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

struct GRID {
    std::vector<string> pois;
};
typedef struct GRID Grid;

struct POI {
    float lat;
    float lng;
};
typedef struct POI Poi;

struct COORDINATE {
    int x;
    int y;
};
typedef struct COORDINATE Coordinate;


namespace utils{
    // data io
    FILE * fopen_(const char* p, const char* m);
    std::ifstream* ifstream_(const char* p);

    // loading POI related data
    Grid ** loadGridInfo(char* infile, int ndimx, int ndimy);

    std::map<std::string, Poi*>* loadPoiInfo(const std::map<std::string, int>* poi_ids, char* poi_path, int data_num);

    void getIdMapRelation(char* infile, std::map<std::string, int>* user_ids,
        std::map<int, std::string>* ruser_ids, std::map<std::string, int>* poi_ids,
        std::map<int, std::string>* rpoi_ids, int *n_users, int *n_pois); 
 
    Coordinate getGridIdxForPoi(Poi* latlng, int ndimx, int ndimy, double grain_lng, double grain_lat);
    
    Coordinate* checkBoundary(int x, int y, int ndimx, int ndimy);
    
    std::vector<Coordinate*>* getNearGridsForGrid(Coordinate grididx, int ndimx, int ndimy);
    
    std::vector<Coordinate*>* getNearGridsForPoi(Poi* latlng, int ndimx, int ndimy, double grain_lng, double grain_lat, bool tag);

    
    // mathematical functions
    inline double logitLoss(double x) {
        return (1-1.0/(1+exp(-x)));
    };

    // random number generator
    double gaussrand(double ep, double var);
    void muldimGaussrand(double ** factor, int ndim);
    void muldimUniform(double ** factor, int ndim);
    void muldimZero(double ** factor, int ndim);
    std::vector<std::string> genSamples(std::vector<std::string>* data, int nsample);

    // char array string split function
    std::vector<char *> split_str(char * in_str, char sep);
    // extract substr
    char * sub_str(int s_idx, int e_idx, char * raw_str);
    // string split function
    std::vector<std::string> split_str(std::string in_str, char sep);
    
    // count the line of file
    int cnt_file_line(std::string in_file);
    
    // allocate matrix memory
    int ** alloc_matrix(int xdim, int ydim);
    // allocate vector memory
    int * alloc_vector(int ndim);

    void pause();
}


FILE * utils::fopen_(const char* p, const char* m) {
    FILE *f = fopen(p, m);
    if (!f) {
        printf("Failed to open %s\n", p);
        exit(1);
    }
    return f;
}


std::ifstream* utils::ifstream_(const char* p){
    std::ifstream *in = new std::ifstream(p);
    if (!in) {
        printf("Failed to open %s\n", p);
        exit(1);
    }
    return in;
}

void utils::getIdMapRelation(char* infile, std::map<std::string, int>* user_ids,
        std::map<int, std::string>* ruser_ids, std::map<std::string, int>* poi_ids,
        std::map<int, std::string>* rpoi_ids, int *n_users, int *n_pois) {
    std::string uid, pid1, pid2;
    std::string line;
    std::vector<std::string> parts;
    
    *n_users=0, *n_pois=0;
    //FILE *f = fopen_(infile, "r");
    std::ifstream* in = utils::ifstream_(infile);
    while(std::getline(*in, line)){
        parts = utils::split_str(line, ',');
        uid = parts[0];
        pid1 = parts[1];
        pid2 = parts[4];
        if ((*user_ids).find(uid) == (*user_ids).end()) {
            (*user_ids)[uid] = *n_users;
            (*ruser_ids)[*n_users] = uid;
            *n_users++;
        }
        if ((*poi_ids).find(pid1) == (*poi_ids).end()) {
            (*poi_ids)[pid1] = *n_pois;
            (*rpoi_ids)[*n_pois] = pid1;
            *n_pois++;
        }
        if ((*poi_ids).find(pid2) == (*poi_ids).end()) {
            (*poi_ids)[pid2] = *n_pois;
            (*rpoi_ids)[*n_pois] = pid2;
            *n_pois++;
        }
    }
}


Grid ** utils::loadGridInfo(char* infile, int ndimx, int ndimy)
{
    int xidx, yidx;
    std::string line;
    std::vector<std::string> parts;
    
    Grid ** grids = new Grid*[ndimx];
    for (int i=0; i<ndimx; i++)
        grids[i] = new Grid[ndimy];
    
    //FILE *f = fopen_(infile, "r");
    std::ifstream* in = utils::ifstream_(infile);
    while(std::getline(*in, line)){
        parts = utils::split_str(line, ',');
        xidx = atoi(parts[0].c_str());
        yidx = atoi(parts[1].c_str());
        for (int i=2; i<parts.size(); i++)
            //grids[xidx][yidx].pois.push_back(atoi(parts[i].c_str()));
            grids[xidx][yidx].pois.push_back(parts[i]);
    }
    
    return grids;
}


std::map<std::string, Poi*>* utils::loadPoiInfo(const std::map<std::string, int>* poi_ids, char* poi_path, int data_num) {
    float lat, lng;
    std::string pid, line;
    std::vector<std::string> parts;
    std::map<std::string, Poi*>* pois_latlng = new std::map<std::string, Poi*>();
    
    //FILE *f = fopen_(poi_path, "r");
    std::ifstream* in = utils::ifstream_(poi_path);
    while(std::getline(*in, line)){
        if (data_num == 0 || data_num == 1) {
            parts = utils::split_str(line, '\t');
            pid = parts[4];
            lat = atof(parts[2].c_str());
            lng = atof(parts[3].c_str());
        }else if (data_num == 2){
            parts = utils::split_str(line, ',');
            pid = parts[0];
            lat = atof(parts[2].c_str());
            lng = atof(parts[3].c_str());
        }else {
            printf("Invalid choice of dataset!\n");
            exit(1);
        }
        Poi * poi = new Poi();
        poi->lat = lat;
        poi->lng = lng;
        (*pois_latlng)[pid] = poi;
    }
    return pois_latlng;
}


Coordinate utils::getGridIdxForPoi(Poi* latlng, int ndimx, int ndimy, double grain_lng, double grain_lat) {
    Coordinate grididx;
    grididx.x = int((latlng->lng+180)/grain_lng);
    grididx.y = int((latlng->lat+90)/grain_lat);
    if (latlng->lng == 180)
        grididx.x = ndimx-1;
    if (latlng->lat == 90)
        grididx.y = ndimy-1;
    return grididx;
}


Coordinate* utils::checkBoundary(int x, int y, int ndimx, int ndimy) {
    Coordinate* coordinate = new Coordinate();
    coordinate->x = x;
    coordinate->y = y;
    if (coordinate->x == -1)
        coordinate->x = ndimx-1;
    else if (coordinate->x == ndimx)
        coordinate->x = 0;
    if (coordinate->y == -1)
        coordinate->y = ndimy-1;
    else if (coordinate->y == ndimy)
        coordinate->y = 0;
    return coordinate;
}


std::vector<Coordinate*>* utils::getNearGridsForGrid(Coordinate grididx, int ndimx, int ndimy){
    std::vector<Coordinate*>* near_grids = new std::vector<Coordinate*>();
    Coordinate* coordinate = NULL;
    
    if (grididx.x < 0 || grididx.x >= ndimx || grididx.y < 0 || grididx.y >= ndimy){
        printf("Invalid grididx!\n");
        exit(1);
    }
    coordinate = utils::checkBoundary(grididx.x-1, griddix.y+1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x, griddix.y+1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x+1, griddix.y+1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x-1, griddix.y, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x+1, griddix.y, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x-1, griddix.y-1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x, griddix.y-1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    coordinate = utils::checkBoundary(grididx.x+1, griddix.y-1, ndimx, ndimy);
    near_grids.push_back(coordinate);
    return near_grids;
}


std::vector<Coordinate*>* utils::getNearGridsForPoi(Poi* latlng, int ndimx, int ndimy, double grain_lng, double grain_lat, bool tag) {
    std::vector<Coordinate*>* near_grids = NULL;
    
    Coordinate grididx = utils::getGridIdxForPoi(latlng, ndimx, ndimy, grain_lng, grain_lat);
    near_grids = utils::getNearGridsForGrid(grididx, ndimx, ndimy);
    if (tag)
        Coordinate* self_grid = &grididx;
        near_grids.push_back(self_grid);
    return near_grids;
}



double utils::gaussrand(double ep, double var) {
    double V1, V2, S;
    int phase = 0;
    double X;
                     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
                                         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
    X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    X = X*var + ep;
    return X;
}


void utils::muldimGaussrand(double ** factor, int ndim) {
    *factor = new double[ndim];
    for (int i=0; i<ndim; i++)
        (*factor)[i] = utils::gaussrand(0.0, 0.1);
}


void utils::muldimUniform(double ** factor, int ndim) {
    *factor = new double[ndim];
    for (int i=0; i<ndim; i++) {
        (*factor)[i] = 2*float(rand())/RAND_MAX-1;
    }
}


void utils::muldimZero(double ** factor, int ndim) {
    *factor = new double[ndim];
    for (int i=0; i<ndim; i++) {
        (*factor)[i] = 0.0;
    }
}


std::vector<std::string>* genNegSamples(std::vector<std::string>* data, std::set<std::string>* filter_samples, int nsample) {
    std::vector<std::string>* neg_samples = NULL;
    vector<std::string>::iterator it = NULL;
  
    it = data->begin();
    while(it != data->end()) {
        if (filter_samples->find(*it) != filter_samples->end())
            it = data->erase(it);
        else
            it++;
    }
    neg_samples = tmp_samples = utils::genSamples(data, tmp_nsample);
    return neg_samples;
}


std::vector<std::string>* utils::genSamples(std::vector<std::string>* data, int nsample) {
    int sampled_num;
    int sample_id;
    std::vector<std::string>* samples = new std::vector<std::string>();
    
    sampled_num = 0;
    while (sampled_num < nsample) {
        if (data.size() == 0)
            break;
        sample_id = rand()%data.size();
        samples.push_back(data[sample_id]);
        data.erase(data.begin()+sample_id);
        sampled_num++;
    }
    return samples;
}


std::vector<char *> utils::split_str(char * in_str, char sep){
    std::vector<char *> str_seg;
    
    int str_len = strlen(in_str);
    if(in_str[str_len - 1] == '\n')
        str_len = str_len - 1;
    
    int s_idx, e_idx;
    s_idx = 0;
    for (int i = 0; i < str_len; i++){
        if(in_str[i] == sep){
            e_idx = i-1;
            str_seg.push_back(utils::sub_str(s_idx, e_idx, in_str));
            s_idx = i+1;
        }
    }
    if (s_idx < str_len)
        str_seg.push_back(utils::sub_str(s_idx, str_len-1, in_str));

    return str_seg;
}

char * utils::sub_str(int s_idx, int e_idx, char * raw_str){
    //char new_str[e_idx+1-s_idx+1];  // first +1: right number, second +1: "\0"
    char * new_str = new char[e_idx+1-s_idx+1];
    memset(new_str, 0, e_idx+1-s_idx+1);

    for (int i=s_idx; i<=e_idx; i++)
        new_str[i-s_idx] = raw_str[i];

    return new_str;
}

std::vector<std::string> utils::split_str(std::string in_str, char sep){
    std::vector<std::string> str_seg;
    
    int str_len = in_str.length();
    if (in_str[str_len-1] == '\n')
        str_len = str_len - 1;
    
    int s_idx, e_idx;
    s_idx = 0;
    for (int i=0; i<str_len; i++){
        if(in_str[i] == sep){
            e_idx = i-1;
            str_seg.push_back(in_str.substr(s_idx, e_idx+1-s_idx));
            s_idx = i+1;
        }
    }
    if (s_idx < str_len)
        str_seg.push_back(in_str.substr(s_idx, str_len-s_idx));

    return str_seg;
}


int * utils::alloc_vector(int ndim){
    int * tmp_vector = NULL;
    
    tmp_vector = new int[ndim];
    for (int i=0; i<ndim; i++){
        tmp_vector[i] = 0;
    }

    return tmp_vector;
}

int ** utils::alloc_matrix(int xdim, int ydim){
    int ** tmp_matrix = NULL;

    tmp_matrix = new int*[xdim];
    for (int i=0; i<xdim; i++){
        tmp_matrix[i] = new int[ydim];
        for (int j=0; j<ydim; j++){
            tmp_matrix[i][j] = 0;
        }
    }

    return tmp_matrix;
}


void utils::pause(){
    printf("Pausing the program, type any character to restart...\n");
    getchar();
}

