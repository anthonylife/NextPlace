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
# Date: 2014/5/3                                                  #
# Probabilistic Matric Factorization with Pairwise Learning for   #
#   implicit feedback data.                                       #
###################################################################

import sys, csv, json, argparse
sys.path.append("../")
import numpy as np
from tool import loadGridInfo, loadPoiInfo, getNearGridsForPOI, getMulMapId, rZero, rGaussian
from data_io import write_submission

with open("../../SETTINGS.json") as fp:
    settings = json.loads(fp.read())

class PMF():
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
            poi_path, grid_path, data_num, init_choice):
        self.trdata_path = trdata_path
        self.vadata_path = vadata_path
        self.tedata_path = tedata_path
        self.poi_path = poi_path
        self.grid_path = grid_path

        self.grids_pois = loadGridInfo(self.grid_path)
        self.pois_latlng = loadPoiInfo(self.poi_path, data_num)
        self.user_ids, self.ruser_ids, self.poi_ids, self.rpoi_ids = getMulMapId(self.trdata_path)

        self.user_factor = np.array()
        init_para = None
        if init_choice == settings["BPR_INIT_ZERO"]:
            init_para = rZero
        elif init_choice == settings["BPR_INIT_GAUSSIAN"]:
            init_para = rGaussian
        else:
            print 'Choice of model initialization error.'
            sys.exit(1)

