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
// Running corresponding model (LCR)                             //
///////////////////////////////////////////////////////////////////

//#include "model.hpp"
#include "model1.hpp"

using namespace std;

string DATA_ROOT_PATH = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/data/";
string TRAIN_SUFFIX[3] = {"Brightkite_Train_Checkinpair.csv",
    "Gowalla1_Train_Checkinpair.csv", "Gowalla2_Train_Checkinpair.csv"};
string VALI_SUFFIX[3] = {"Brightkite_Vali_Checkinpair.csv", 
    "Gowalla1_Vali_Checkinpair.csv", "Gowalla2_Vali_Checkinpair.csv"};
//string TEST_SUFFIX[3] = {"Brightkite_Test_Checkinpair.csv", 
//    "Gowalla1_Test_Checkinpair.csv", "Gowalla2_Test_Checkinpair.csv"};
string TEST_SUFFIX[3] = {"Brightkite_Filter_Test_Checkinpair.csv",
    "Gowalla1_Filter_Test_Checkinpair.csv","Gowalla2_Filter_Test_Checkinpair.csv"};
string POI_SUFFIX[3] = {"Brightkite_totalCheckins.txt", 
    "Gowalla_totalCheckins.txt", "Gowalla_places.csv"};
string GRID_SUFFIX[3] = {"Brightkite_Grid_Place.csv",
    "Gowalla1_Grid_Place.csv", "Gowalla2_Grid_Place.csv"};


int ArgPos(char *str, int argc, char **argv) {
    int a;
    for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
        if (a == argc - 1) {
            printf("Argument missing for %s\n", str);
            exit(1);
        }
        return a;
    }
    return -1;
}


int main(int argc, char **argv) {
    int i;
    int a=0;
    char *b=NULL, *c=NULL;
    if (argc == 1) {
        printf("LCR v 0.1a\n");
        printf("\tExamples:\n");
        printf("./run -d 0(1,2) -r True(False) -b True(False)\n\n");
        return 0;
    }
    if ((i = ArgPos((char *)"-d", argc, argv)) > 0) a = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-r", argc, argv)) > 0) b = argv[i + 1];
    if ((i = ArgPos((char *)"-b", argc, argv)) > 0) c = argv[i + 1];
    if (a!=0 && a!=1 && a!=2) {
        printf("Invalid choice of dataset!\n");
        exit(1);
    }
    
    string submission_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/results/LCR_Result.dat";
    string user_factor_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/LCR_C/user_factor.model";
    string poi_factor_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/LCR_C/poi_factor.model";
    string query_factor_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/LCR_C/query_factor.model";
    string poi_bias_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/LCR_C/poi_bias.model";
    string trdata_path = DATA_ROOT_PATH + TRAIN_SUFFIX[a];
    string vadata_path = DATA_ROOT_PATH + VALI_SUFFIX[a];
    string tedata_path = DATA_ROOT_PATH + TEST_SUFFIX[a];
    string query_path = DATA_ROOT_PATH + "Gowalla_places.csv";
    string poi_path = DATA_ROOT_PATH + POI_SUFFIX[a];
    string grid_path = DATA_ROOT_PATH + GRID_SUFFIX[a];
    int re_topk = 10;
    int data_num = a;
    bool restart_tag, bias_tag;
    
    if (strcmp(b, (char *)"True") == 0)
        restart_tag = true;
    else if (strcmp(b, (char *)"False") == 0)
        restart_tag = false;
    else {
        printf("Invalid input of para -r\n");
        exit(1);
    }
    if (strcmp(c, (char *)"True") == 0)
        bias_tag = true;
    else if (strcmp(c, (char *)"False") == 0)
        bias_tag = false;
    else {
        printf("Invalid input of para -b\n");
        exit(1);
    }
   
    timeval start_t, end_t;
    utils::tic(start_t);
    
    LCR *lcr = new LCR((char *)trdata_path.c_str(),
                          (char *)vadata_path.c_str(),
                          (char *)tedata_path.c_str(),
                          (char *)query_path.c_str(),
                          (char *)user_factor_path.c_str(),
                          (char *)poi_factor_path.c_str(),
                          (char *)query_factor_path.c_str(),
                          (char *)poi_bias_path.c_str(),
                          (char *)poi_path.c_str(),
                          (char *)grid_path.c_str(),
                          data_num,
                          bias_tag,
                          restart_tag,
                          re_topk);
    if (restart_tag)
        lcr->train();
    lcr->recommendationNewPOI((char *)submission_path.c_str());
    
    utils::toc(start_t, end_t);

    return 0;
}
