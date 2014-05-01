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
# Providing functions controling input and output for each algori-#
#   thms                                                          #
###################################################################

import sys, csv, json
import pickle

def write_submission(user_recommend_result, result_path):
    wfd = open(result_path, "w")
    id_sorted_result = sorted(user_recommend_result.items(), key=lambda x: x[0])
    for entry in id_sorted_result:
        wfd.write("%d\t" % entry[0])
        sorted_pids = sorted(entry[1])
        for i, pid in enumerate(sorted_pids):
            if i == len(sorted_pids)-1:
                wfd.write("%d\n" % pid)
            else:
                wfd.write("%d," % pid)
    wfd.close()

def save_model(model, out_path):
    pickle.dump(model, open(out_path, "w"))

def load_model(in_path):
    return pickle.load(open(in_path))

