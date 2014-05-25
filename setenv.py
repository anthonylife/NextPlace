#!/usr/bin/env python
#encoding=utf8

import os, json

if __name__ == "__main__":
    settings = {}
    setting_file = "SETTINGS.json"
    wfp = open(setting_file, "w")
    settings["ROOT_PATH"] = os.getcwd() + "/"

    ### Configuration by user
    settings["SRC_DATA_FILE1_1"] = "data/Brightkite_totalCheckins.txt"
    settings["SRC_DATA_FILE1_2"] = "data/Brightkite_edges.txt"
    settings["SRC_DATA_FILE2_1"] = "data/Gowalla_totalCheckins.txt"
    settings["SRC_DATA_FILE2_2"] = "data/Gowalla_edges.txt"
    settings["SRC_DATA_FILE3_1"] = "data/Gowalla_checkins.csv"
    settings["SRC_DATA_FILE3_2"] = "data/Gowalla_friendship.csv"
    settings["SRC_DATA_FILE3_3"] = "data/Gowalla_places.csv"
    settings["CHECKIN_PAIR_FILE1"] = "data/Brightkite_checkinpair.csv"
    settings["CHECKIN_PAIR_FILE2"] = "data/Gowalla1_checkinpair.csv"
    settings["CHECKIN_PAIR_FILE3"] = "data/Gowalla2_checkinpair.csv"
    settings["FILTER_CHECKIN_PAIR_FILE1"] = "data/Brightkite_Filter_Checkinpair.csv"
    settings["FILTER_CHECKIN_PAIR_FILE2"] = "data/Gowalla1_Filter_Checkinpair.csv"
    settings["FILTER_CHECKIN_PAIR_FILE3"] = "data/Gowalla2_Filter_Checkinpair.csv"
    settings["FILTER_USER_VISIT_NUM"] = 10
    settings["FILTER_LOCATION_VISIT_NUM"] = 5
    settings["FILTER_USER_RECORD_NUM"] = 10
    settings["FILTER_LOCATION_RECORD_NUM"] = 5
    settings["TIME_SECONDS_THRESHOLD"] = 7200

    settings["TRAIN_PAIR_FILE1"] = "data/Brightkite_Train_Checkinpair.csv"
    settings["VALI_PAIR_FILE1"] = "data/Brightkite_Vali_Checkinpair.csv"
    settings["TEST_PAIR_FILE1"] = "data/Brightkite_Test_Checkinpair.csv"
    settings["FILTER_TEST_PAIR_FILE1"] = "data/Brightkite_Filter_Test_Checkinpair.csv"
    settings["CV_PAIR_FILE1"] = "data/Brightkite_CV_Checkinpair.csv"
    settings["TRAIN_PAIR_FILE2"] = "data/Gowalla1_Train_Checkinpair.csv"
    settings["VALI_PAIR_FILE2"] = "data/Gowalla1_Vali_Checkinpair.csv"
    settings["TEST_PAIR_FILE2"] = "data/Gowalla1_Test_Checkinpair.csv"
    settings["FILTER_TEST_PAIR_FILE2"] = "data/Gowalla1_Filter_Test_Checkinpair.csv"
    settings["CV_PAIR_FILE2"] = "data/Gowalla1_CV_Checkinpair.csv"
    settings["TRAIN_PAIR_FILE3"] = "data/Gowalla2_Train_Checkinpair.csv"
    settings["VALI_PAIR_FILE3"] = "data/Gowalla2_Vali_Checkinpair.csv"
    settings["TEST_PAIR_FILE3"] = "data/Gowalla2_Test_Checkinpair.csv"
    settings["FILTER_TEST_PAIR_FILE3"] = "data/Gowalla2_Filter_Test_Checkinpair.csv"
    settings["CV_PAIR_FILE3"] = "data/Gowalla2_CV_Checkinpair.csv"
    settings["TRAIN_RATIO"] = 0.7
    settings["VALI_RATIO"] = 0.1
    settings["TEST_RATIO"] = 0.2
    settings["CV_NUM"] = 10
    settings["NEAR_PLACE_FILE1"] = "data/Brightkite_Near_Place.csv"
    settings["NEAR_PLACE_FILE2"] = "data/Gowalla1_Near_Place.csv"
    settings["NEAR_PLACE_FILE3"] = "data/Gowalla2_Near_Place.csv"
    settings["DISTANCE_THRESHOLD"] = 10
    settings["GRID_PLACE_FILE1"] = "data/Brightkite_Grid_Place.csv"
    settings["GRID_PLACE_FILE2"] = "data/Gowalla1_Grid_Place.csv"
    settings["GRID_PLACE_FILE3"] = "data/Gowalla2_Grid_Place.csv"
    settings["GRID_LAT"] = 0.05
    settings["GRID_LNG"] = 0.05

    settings["TOPK"] = 10
    settings["FRIEND_SIMVAL_FILE1"] = "data/Brightkite_Friend_Simval.csv"
    settings["FRIEND_SIMVAL_FILE2"] = "data/Gowalla1_Friend_Simval.csv"
    settings["FRIEND_SIMVAL_FILE3"] = "data/Gowalla2_Friend_Simval.csv"

    settings["MAX_TOPK"] = 20
    settings["NUM_FRIEND"] = 10

    settings["NUM_LOGIT"] = 100000
    settings["MAX_LOGIT"] = 6

    # Popular-based Method
    settings["POPULAR_SUBMISSION_PATH"] = "results/Popular_Result.dat"

    # Personal Popular-based Method
    settings["PER_POPULAR_SUBMISSION_PATH"] = "results/Per_Popular_Result.dat"

    # PMF Method
    settings["PMF_USER_FACTOR_FILE"] = "user_factor.model"
    settings["PMF_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["PMF_POI_BIAS_FILE"] = "poi_bias.model"
    settings["PMF_INIT_GAUSSIAN"] = "gaussian"
    settings["PMF_INIT_ZERO"] = "zero"
    settings["PMF_SUBMISSION_PATH"] = "results/PMF_Result.dat"

    # FPMC Method
    settings["FPMC_USER_FACTOR_FILE"] = "user_factor.model"
    settings["FPMC_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["FPMC_POI_BIAS_FILE"] = "poi_bias.model"
    settings["FPMC_INIT_GAUSSIAN"] = "gaussian"
    settings["FPMC_INIT_ZERO"] = "zero"
    settings["FPMC_SUBMISSION_PATH"] = "results/FPMC_Result.dat"

    # SFPMC Method
    settings["SFPMC_USER_FACTOR_FILE"] = "user_factor.model"
    settings["SFPMC_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["SFPMC_POI_BIAS_FILE"] = "poi_bias.model"
    settings["SFPMC_INIT_GAUSSIAN"] = "gaussian"
    settings["SFPMC_INIT_ZERO"] = "zero"
    settings["SFPMC_SUBMISSION_PATH"] = "results/SFPMC_Result.dat"

    # LCR Method
    settings["LCR_USER_FACTOR_FILE"] = "user_factor.model"
    settings["LCR_USER_CONTEXT_FILE"] = "user_context.model"
    settings["LCR_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["LCR_QUERY_FACTOR_FILE"] = "category_factor.model"
    settings["LCR_POI_BIAS_FILE"] = "poi_bias.model"
    settings["LCR_CONTEXT_FACTOR_FILE"] = "context_factor.model"
    settings["LCR_INIT_GAUSSIAN"] = "gaussian"
    settings["LCR_INIT_ZERO"] = "zero"
    settings["LCR_SUBMISSION_PATH"] = "results/LCR_Result.dat"

    # SCR Method
    settings["SCR_USER_FACTOR_FILE"] = "user_factor.model"
    settings["SCR_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["SCR_QUERY_FACTOR_FILE"] = "category_factor.model"
    settings["SCR_POI_BIAS_FILE"] = "poi_bias.model"
    settings["SCR_INIT_GAUSSIAN"] = "gaussian"
    settings["SCR_INIT_ZERO"] = "zero"
    settings["SCR_SUBMISSION_PATH"] = "results/SCR_Result.dat"

    # LLCR Method
    settings["LLCR_USER_FACTOR_FILE"] = "user_factor.model"
    settings["LLCR_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["LLCR_QUERY_FACTOR_FILE"] = "category_factor.model"
    settings["LLCR_POI_BIAS_FILE"] = "poi_bias.model"
    settings["LLCR_INIT_GAUSSIAN"] = "gaussian"
    settings["LLCR_INIT_ZERO"] = "zero"
    settings["LLCR_SUBMISSION_PATH"] = "results/LLCR_Result.dat"

    # LSCR Method
    settings["LSCR_USER_FACTOR_FILE"] = "user_factor.model"
    settings["LSCR_POI_FACTOR_FILE"] = "poi_factor.model"
    settings["LSCR_QUERY_FACTOR_FILE"] = "category_factor.model"
    settings["LSCR_POI_BIAS_FILE"] = "poi_bias.model"
    settings["LSCR_INIT_GAUSSIAN"] = "gaussian"
    settings["LSCR_INIT_ZERO"] = "zero"
    settings["LSCR_SUBMISSION_PATH"] = "results/LSCR_Result.dat"

    # TFPMC Method
    settings["TFPMC_SUBMISSION_PATH"] = "results/TFPMC_Result.dat"

    # TLCR Method

    # TSCR Method
    settings["TLCR_SUBMISSION_PATH"] = "results/TLCR_Result.dat"


    # TLSCR Method
    settings["TLSCR_SUBMISSION_PATH"] = "results/TLSCR_Result.dat"

    # Evaluation method
    settings["F-score"] = "F-score"

    # Write result
    json.dump(settings, wfp, sort_keys=True, indent=4)
