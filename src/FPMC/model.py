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
# Date: 2014/5/7                                                  #
# Factorized Personalized Markov Chain with Pairwise Learning for #
#   implicit feedback data.                                       #
###################################################################

import sys, csv, json, argparse, random
sys.path.append("../")
import numpy as np
from collections import defaultdict
from tool import loadGridInfo, loadPoiInfo, getNearGridsForPOI, getMulMapId, rZero, rGaussian
from tool import logitLoss
from data_io import write_submission

with open("../../SETTINGS.json") as fp:
    settings = json.loads(fp.read())

class FPMC():
    def __init__(self):
        # Method control variable
        self.niters = 50
        self.nsample = 5

        # Hyper-parameter setting
        self.ndim = 10
        self.lr = 0.05
        self.u_reg = 0.1
        self.p_reg = 0.1
        self.bp_reg = 0.1       #bias reg of poi

    def model_init(self, trdata_path, vadata_path, tedata_path,
            poi_path, grid_path, data_num, init_choice, bias_tag, restart_tag):
        self.trdata_path = trdata_path
        self.vadata_path = vadata_path
        self.tedata_path = tedata_path
        self.poi_path = poi_path
        self.grid_path = grid_path
        self.bias_tag = bias_tag

        self.grids_pois = loadGridInfo(self.grid_path)
        self.pois_latlng = loadPoiInfo(self.poi_path, data_num)
        init_para = None
        if init_choice == settings["FPMC_INIT_ZERO"]:
            init_para = rZero
        elif init_choice == settings["FPMC_INIT_GAUSSIAN"]:
            init_para = rGaussian
        else:
            print 'Choice of model initialization error.'
            sys.exit(1)

        if restart_tag == True:
            self.user_ids, self.ruser_ids, self.poi_ids, self.rpoi_ids = getMulMapId(self.trdata_path)
            self.user_factor = np.array([init_para(self.ndim) for j in xrange(len(self.user_ids))])
            self.poi_factor = np.array([init_para(self.ndim) for j in xrange(len(self.poi_ids))])
            if self.bias_tag == True:
                self.poi_bias = np.array([0.0 for i in xrange(self.poi_ids)])
        else:
            self.user_ids, self.ruser_ids, self.poi_ids, self.rpoi_ids, self.user_factor,\
                    self.poi_factor, self.poi_bias = self.load_model()

    def genTrainTriples(self):
        tr_triple = []
        index_extent = (-90, -180, 90, 180)
        ndimx = int((index_extent[3]-index_extent[1])/settings["GRID_LNG"])
        ndimy = int((index_extent[2]-index_extent[0])/settings["GRID_LAT"])
        idx = 0
        for entry in csv.reader(open(self.trdata_path)):
            sys.stdout.write("\rFINISHED TRIPLE NUM: %d. " % (idx+1))
            sys.stdout.flush()
            idx += 1
            uid, pid1, pid2 = int(entry[0]), int(entry[1]), int(entry[4])
            near_grids = getNearGridsForPOI(self.pois_latlng[pid2], ndimx, ndimy, True)
            candidate_pois = []
            for grididx in near_grids:
                candidate_pois += self.grids_pois[grididx[0]][grididx[1]]
            if self.nsample < len(candidate_pois):
                for pid in random.sample(set(candidate_pois)-set([pid2]), self.nsample):
                    tr_triple.append([self.user_ids[uid], self.poi_ids[pid1], self.poi_ids[pid2], self.poi_ids[pid]])
        print len(tr_triple)
        return tr_triple

    def train(self):
        self.tr_triples = self.genTrainTriples()
        for i in xrange(self.niters):
            random.shuffle(self.tr_triples)
            for j, triple in enumerate(self.tr_triples):
                if self.bias_tag == True:
                    pos_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[2]])\
                              + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[2]])\
                              + self.poi_bias[triple[2]]
                    neg_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[3]])\
                              + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[3]])\
                              + self.poi_bias[triple[3]]
                else:
                    pos_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[2]])\
                              + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[2]])
                    neg_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[3]])\
                              + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[3]])
                logit_loss = logitLoss(pos_score, neg_score)
                tmp_diff1 = self.poi_factor[triple[2]]-self.poi_factor[triple[3]]
                self.user_factor[triple[0]] = self.user_factor[triple[0]]+self.lr*(logit_loss*tmp_diff1-self.u_reg*self.user_factor[triple[0]])
                self.poi_factor[triple[1]] = self.poi_factor[triple[1]]+self.lr*(logit_loss*tmp_diff1-self.p_reg*self.poi_factor[triple[1]])
                tmp_diff2 = self.user_factor[triple[0]] + self.poi_factor[triple[1]]
                self.poi_factor[triple[2]] = self.poi_factor[triple[2]]+self.lr*(logit_loss*tmp_diff2-self.p_reg*self.poi_factor[triple[2]])
                self.poi_factor[triple[3]] = self.poi_factor[triple[3]]+self.lr*(logit_loss*(-tmp_diff2)-self.p_reg*self.poi_factor[triple[3]])
                if self.bias_tag == True:
                    self.poi_bias[triple[2]] = self.poi_bias[triple[2]]+self.lr*(logit_loss*-self.bp_reg*self.poi_bias[triple[2]])
                    self.poi_bias[triple[3]] = self.poi_bias[triple[3]]+self.lr*(logit_loss*-self.bp_reg*self.poi_bias[triple[3]])
                sys.stdout.write("\rFINISHED TRIPLE NUM: %d. " % (j+1))
                sys.stdout.flush()
            print "\nCurrent iteration %d, AUC is %f...\n" % (i+1, self.evaluation())
            #raw_input()
        self.tr_triples = None
        self.save_model()

    def evaluation(self):
        correct_triple = 0
        for triple in self.tr_triples:
            if self.bias_tag == True:
                pos_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[2]])\
                          + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[2]])\
                          + self.poi_bias[triple[2]]
                neg_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[3]])\
                          + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[3]])\
                          + self.poi_bias[triple[3]]
            else:
                pos_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[2]])\
                          + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[2]])
                neg_score = np.dot(self.user_factor[triple[0]], self.poi_factor[triple[3]])\
                          + np.dot(self.poi_factor[triple[1]], self.poi_factor[triple[3]])

            if pos_score > neg_score:
                correct_triple += 1
        return 1.0*correct_triple/len(self.tr_triples)

    def recommendation(self, submission_path):
        index_extent = (-90, -180, 90, 180)
        ndimx = int((index_extent[3]-index_extent[1])/settings["GRID_LNG"])
        ndimy = int((index_extent[2]-index_extent[0])/settings["GRID_LAT"])
        recommendation_result = {}
        cache_user_poi_score = defaultdict(dict)
        for i, entry in enumerate(csv.reader(open(self.tedata_path))):
            uid, pid1 = int(entry[0]), int(entry[1])
            near_grids = getNearGridsForPOI(self.pois_latlng[pid1], ndimx, ndimy, True)
            candidate_pois = []
            for grididx in near_grids:
                candidate_pois += self.grids_pois[grididx[0]][grididx[1]]
            result = []
            pois_score = []
            for c_pid in candidate_pois:
                if self.bias_tag == True:
                    score = np.dot(self.user_factor[self.user_ids[uid]],self.poi_factor[self.poi_ids[c_pid]])\
                          + np.dot(self.poi_factor[self.poi_ids[pid1]], self.poi_factor[self.poi_ids[c_pid]])\
                          + self.poi_bias[self.poi_ids[c_pid]]
                else:
                    score = np.dot(self.user_factor[self.user_ids[uid]],+ self.poi_factor[self.poi_ids[c_pid]])\
                          + np.dot(self.poi_factor[self.poi_ids[pid1]], self.poi_factor[self.poi_ids[c_pid]])
                    pois_score.append([c_pid, score])
                    cache_user_poi_score[uid][c_pid] = score
            result = sorted(pois_score, key=lambda x:x[1], reverse=True)[:settings["MAX_TOPK"]]
            recommendation_result[i] = [pair[0] for pair in result]
            sys.stdout.write("\rFINISHED PAIR NUM: %d. " % (i+1))
            sys.stdout.flush()
        write_submission(recommendation_result, submission_path)

    def recommendationNewPOI(self, submission_path):
        index_extent = (-90, -180, 90, 180)
        ndimx = int((index_extent[3]-index_extent[1])/settings["GRID_LNG"])
        ndimy = int((index_extent[2]-index_extent[0])/settings["GRID_LAT"])
        recommendation_result = {}
        cache_user_poi_score = defaultdict(dict)
        user_visited = defaultdict(set)
        for entry in csv.reader(open(self.trdata_path)):
            uid, pid1, pid2 = int(entry[0]), int(entry[1]), int(entry[4])
            user_visited[uid].add(pid1)
            user_visited[uid].add(pid2)
        for i, entry in enumerate(csv.reader(open(self.tedata_path))):
            uid, pid1 = int(entry[0]), int(entry[1])
            near_grids = getNearGridsForPOI(self.pois_latlng[pid1], ndimx, ndimy, True)
            candidate_pois = []
            for grididx in near_grids:
                candidate_pois += self.grids_pois[grididx[0]][grididx[1]]
            result = []
            pois_score = []
            for c_pid in set(candidate_pois)-user_visited[uid]:
                if self.bias_tag == True:
                    score = np.dot(self.user_factor[self.user_ids[uid]],self.poi_factor[self.poi_ids[c_pid]])\
                          + np.dot(self.poi_factor[self.poi_ids[pid1]], self.poi_factor[self.poi_ids[c_pid]])\
                          + self.poi_bias[self.poi_ids[c_pid]]
                else:
                    score = np.dot(self.user_factor[self.user_ids[uid]],+ self.poi_factor[self.poi_ids[c_pid]])\
                          + np.dot(self.poi_factor[self.poi_ids[pid1]], self.poi_factor[self.poi_ids[c_pid]])
                    pois_score.append([c_pid, score])
                    cache_user_poi_score[uid][c_pid] = score
            result = sorted(pois_score, key=lambda x:x[1], reverse=True)[:settings["MAX_TOPK"]]
            recommendation_result[i] = [pair[0] for pair in result]
            sys.stdout.write("\rFINISHED PAIR NUM: %d. " % (i+1))
            sys.stdout.flush()
        write_submission(recommendation_result, submission_path)

    def save_model(self):
        writer = csv.writer(open(settings["FPMC_USER_FACTOR_FILE"], "w"), lineterminator="\n")
        for i in xrange(len(self.user_factor)):
            writer.writerow([self.ruser_ids[i]]+list(self.user_factor[i]))
        writer = csv.writer(open(settings["FPMC_POI_FACTOR_FILE"], "w"), lineterminator="\n")
        for i in xrange(len(self.product_factor)):
            writer.writerow([self.rpoi_ids[i]]+list(self.poi_factor[i]))
        if self.bias_tag == True:
            writer = csv.writer(open(settings["FPMC_POI_BIAS_FILE"], "w"), lineterminator="\n")
            for i in xrange(len(self.poi_bias)):
                writer.writerow([self.rpoi_ids[i]]+list(self.poi_bias[i]))

    def load_model(self):
        self.user_ids = {}
        self.ruser_ids = {}
        u_cnt = 0
        for entry in csv.reader(open(settings["FPMC_USER_FACTOR_FILE"])):
            uid = int(entry[0])
            self.user_ids[uid] = u_cnt
            self.ruser_ids[u_cnt] = uid
            u_cnt += 1
        self.user_factor = np.array([rZero(xrange(self.ndim))
            for j in xrange(len(self.user_ids))])
        for entry in csv.reader(open(settings["FPMC_USER_FACTOR_FILE"])):
            uid = int(entry[0])
            self.user_factor[self.user_ids[uid]] = np.array(map(float, entry[1:]))

        self.poi_ids = {}
        self.rpoi_ids = {}
        p_cnt = 0
        for entry in csv.reader(open(settings["FPMC_POI_FACTOR_FILE"])):
            pid = int(entry[0])
            self.poi_ids[pid] = p_cnt
            self.rpoi_ids[p_cnt] = pid
            p_cnt += 1
        self.poi_factor = np.array([rZero(xrange(self.ndim))
            for j in xrange(len(self.poi_ids))])
        for entry in csv.reader(open(settings["FPMC_POI_FACTOR_FILE"])):
            pid = int(entry[0])
            self.poi_factor[self.poi_ids[pid]] = np.array(map(float, entry[1:]))

        if self.bias_tag == True:
            self.poi_bias = np.array([0.0 for i in xrange(self.poi_ids)])
            for entry in csv.reader(open(settings["FPMC_POI_BIAS_FILE"])):
                pid, bias = int(entry[0]), float(entry[1])
                self.poi_bias[self.poi_ids[pid]] = bias


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    parser.add_argument('-r', type=str, action='store',
            dest='restart', help='whether to retrain the model or not')
    parser.add_argument('-b', type=str, action='store',
            dest='bias', help='whether to retrain the model or not')
    #if len(sys.argv) != 7:
    #    print 'Command e.g.: python model.py -d 0(1,2) -r True(False) -b True(False)'
    #    sys.exit(1)
    para = parser.parse_args()
    if para.data_num == 0:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE1"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE1"]
        #tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE1"]
        tedata_path = settings["ROOT_PATH"] + settings["FILTER_TEST_PAIR_FILE1"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE1"]
    elif para.data_num == 1:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE2"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE2"]
        #tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE2"]
        tedata_path = settings["ROOT_PATH"] + settings["FILTER_TEST_PAIR_FILE2"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE2"]
    elif para.data_num == 2:
        trdata_path = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE3"]
        vadata_path = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE3"]
        #tedata_path = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE3"]
        tedata_path = settings["ROOT_PATH"] + settings["FILTER_TEST_PAIR_FILE3"]
        poi_path = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_3"]
        grid_path = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE3"]
    submission_path = settings["ROOT_PATH"] + settings["POPULAR_SUBMISSION_PATH"]


    restart_tag = True
    bias_tag = False
    if para.restart == 'False' or para.restart == 'F':
        restart_tag = False
    if para.bias == 'True' or para.bias == 'T':
        bias_tag = True

    fpmc = FPMC()
    fpmc.model_init(trdata_path,
                    vadata_path,
                    tedata_path,
                    poi_path,
                    grid_path,
                    para.data_num,
                    settings["FPMC_INIT_GAUSSIAN"],
                    bias_tag,
                    restart_tag)
    fpmc.train()
    #pmf.recommendation(submission_path)
    fpmc.recommendationNewPOI(submission_path)

if __name__ == "__main__":
    main()
