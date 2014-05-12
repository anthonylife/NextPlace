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
// Date: 2014/5/12                                               //
// Running corresponding model (PMF)                             //
///////////////////////////////////////////////////////////////////

#include "model.hpp"

using namespace std;

const string DATA_ROOT_PATH = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/data/";
const string[3] TRAIN_SUFFIX = {"Brightkite_Train_Checkinpair.csv",
    "Gowalla1_Train_Checkinpair.csv", "Gowalla2_Train_Checkinpair.csv"};
const string[3] VALI_SUFFIX = {"Brightkite_Vali_Checkinpair.csv", 
    "Gowalla1_Vali_Checkinpair.csv", "Gowalla2_Vali_Checkinpair.csv"};
//const string[3] TEST_SUFFIX = {"Brightkite_Test_Checkinpair.csv", 
//    "Gowalla1_Test_Checkinpair.csv", "Gowalla2_Test_Checkinpair.csv"};
const string[3] TEST_SUFFIX = {"Brightkite_Filter_Test_Checkinpair.csv",
    "Gowalla1_Filter_Test_Checkinpair.csv","Gowalla2_Filter_Test_Checkinpair.csv"};
const string[3] POI_SUFFIX = {"Brightkite_totalCheckins.txt", 
    "Gowalla_totalCheckins.txt", "Gowalla_places.csv"};
const string[3] GRID_SUFFIX = {"Brightkite_Grid_Place.csv",
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
    int a;
    char* b, c;
    if (argc == 1) {
        printf("PMF v 0.1a\n");
        printf("\tExamples:\n");
        printf("./run -d 0(1,2) -r True(False) -b True(False)\n\n");
        return 0;
    }
    if ((i = ArgPos((char *)"-d", argc, argv)) > 0) a = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-r", argc, argv)) > 0) b = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-b", argc, argv)) > 0) c = atoi(argv[i + 1]);
    if (a!=0 || a!=1 || a!=2) {
        printf("Invalid choice of dataset!\n");
        exit(1);
    }
    
    const string submission_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/results/PMF_result.dat";
    const string user_factor_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/PMF_C/user_factor.model";
    const string poi_factor_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/PMF_C/poi_factor.model";
    const string poi_bias_path = "/home/anthonylife/Doctor/Code/MyPaperCode/NextLocation/src/PMF_C/poi_bias.model";
    string trdata_path = DATA_ROOT_PATH + TRAIN_SUFFIX[a];
    string vadata_path = DATA_ROOT_PATH + VALI_SUFFIX[a];
    string tedata_path = DATA_ROOT_PATH + TEST_SUFFIX[a];
    string poi_path = DATA_ROOT_PATH + POI_SUFFIX[a];
    string grid_path = DATA_ROOT_PATH + GRID_SUFFIX[a];
    const int re_topk = 10;
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
   
    PMF *pmf = new PMF(trdata_path.c_str(),
                       vadata_path.c_str(),
                       tedata_path.c_str(),
                       user_factor_path.c_str(),
                       poi_factor_path.c_str(),
                       poi_bias_path.c_str(),
                       poi_path.c_str(),
                       grid_path.c_str(),
                       data_num,
                       bias_tag,
                       restart_tag,
                       re_topk);
    if (restart_tag)
        pmf->train();
    pmf->recommendationNewPOI(submission_path.c_str());

    return 0;
}
