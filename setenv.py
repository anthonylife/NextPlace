#!/usr/bin/env python
#encoding=utf8

import os, json

if __name__ == "__main__":
    settings = {}
    setting_file = "SETTINGS.json"
    wfp = open(setting_file, "w")
    settings["ROOT_PATH"] = os.getcwd() + "/"

    # Configuration
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
    settings["CV_PAIR_FILE1"] = "data/Brightkite_CV_Checkinpair.csv"
    settings["TRAIN_PAIR_FILE2"] = "data/Gowalla1_Train_Checkinpair.csv"
    settings["VALI_PAIR_FILE2"] = "data/Gowalla1_Vali_Checkinpair.csv"
    settings["TEST_PAIR_FILE2"] = "data/Gowalla1_Test_Checkinpair.csv"
    settings["CV_PAIR_FILE2"] = "data/Gowalla1_CV_Checkinpair.csv"
    settings["TRAIN_PAIR_FILE3"] = "data/Gowalla2_Train_Checkinpair.csv"
    settings["VALI_PAIR_FILE3"] = "data/Gowalla2_Vali_Checkinpair.csv"
    settings["TEST_PAIR_FILE3"] = "data/Gowalla2_Test_Checkinpair.csv"
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

    # Write result
    json.dump(settings, wfp, sort_keys=True, indent=4)
