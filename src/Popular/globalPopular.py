#!/usr/bin/env python
#encoding=utf8

#Copyright [2014] [Wei Zhang]

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

###################################################################
# Date: 2014/5/1                                                  #
# Recommendation Result based on global POI popularity            #
# Note: this popularity methods should consider location informat-#
#   ion. In this task, we restrict the candidate within nearest 9 #
#   grid.                                                         #
###################################################################

import sys, csv, json, argparse
sys.path.append("../")
from collections import defaultdict
from tool import loadGridInfo, loadPoiInfo, getGridIdxForPOI, getNearGrids, getNearGridsForPOI
from data_io import write_submission

with open("../../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


class PopularityRE:
    def __init__(self, trdata_path, vadata_path, tedata_path, poi_path, grid_path, data_num):
        self.trdata_path = trdata_path
        self.vadata_path = vadata_path
        self.tedata_path = tedata_path
        self.poi_path = poi_path
        self.grid_path = grid_path
        self.loadData(data_num)

    def loadData(self, data_num):
        self.grids_pois = loadGridInfo(self.grid_path)
        print len(self.grids_pois)
        self.pois_latlng = loadPoiInfo(self.poi_path, data_num)

    def train(self):
        self.pois_popularity = {}
        for entry in csv.reader(open(self.trdata_path)):
            pid1, pid2 = int(entry[1]), int(entry[4])
            if pid1 not in self.pois_popularity:
                self.pois_popularity[pid1] = 1
            else:
                self.pois_popularity[pid1] += 1
            if pid2 not in self.pois_popularity:
                self.pois_popularity[pid2] = 1
            else:
                self.pois_popularity[pid2] += 1

    def recommendation(self, submission_path):
        index_extent = (-90, -180, 90, 180)
        ndimx = int((index_extent[3]-index_extent[1])/settings["GRID_LNG"])
        ndimy = int((index_extent[2]-index_extent[0])/settings["GRID_LAT"])
        recommendation_result = {}
        for i, entry in enumerate(csv.reader(open(self.tedata_path))):
            uid, pid1 = int(entry[0]), int(entry[1])
            near_grids = getNearGridsForPOI(self.pois_latlng[pid1], ndimx, ndimy, True)
            candidate_pois = []
            for grididx in near_grids:
                candidate_pois += self.grids_pois[grididx[0]][grididx[1]]
            pois_score = [[poi, self.pois_popularity[poi]] for poi in candidate_pois]
            result = sorted(pois_score, key=lambda x:x[1], reverse=True)[:settings["MAX_TOPK"]]
            recommendation_result[i] = [pair[0] for pair in result]
            print i
        write_submission(recommendation_result, submission_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    if len(sys.argv) != 3:
        print 'Command e.g.: python globalPopular.py -d 0(1,2)'
        sys.exit(1)
    para = parser.parse_args()
    if para.data_num == 0:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE1"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE1"]
        tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE1"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE1"]
    elif para.data_num == 1:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE2"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE2"]
        tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE2"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE2"]
    elif para.data_num == 2:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE3"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE3"]
        tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE3"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_3"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE3"]
    submission_path = settings["ROOT_PATH"] + settings["POPULAR_SUBMISSION_PATH"]

    popularityRE = PopularityRE(trdata_path,
                                vadata_path,
                                tedata_path,
                                poi_path,
                                grid_path,
                                para.data_num)

    popularityRE.train()
    popularityRE.recommendation(submission_path)

if __name__ == "__main__":
    main()
