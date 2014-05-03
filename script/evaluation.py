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
# Evaluation on file results                                      #
# Note:                                                           #
#   1. evaluation metrics including Precision, Recall, F-score.   #
###################################################################

import sys, csv, json, argparse
from collections import defaultdict

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


class Evaluation():
    def __init__(self):
        pass

    def evaluate(self, standard_rsult, recommendation_result, method, topk):
        if method == settings["F-score"]:
            result = self.evaluateFscore(standard_rsult, recommendation_result, topk)
            print 'Evaluation result: Precision=%.4f, Recall=%.4f, F-score=%.4f!\n' % result
        else:
            print 'Invalid choice of evaluation method!'
            sys.exit(1)

    def evaluateFscore(self, standard_result, recommendation_result, topk):
        if len(recommendation_result.keys()) != len(standard_result.keys()):
            print 'Number of recommendation items and number of standard results mismatch!'
            sys.exit(1)
        sum_pred_num = 0
        true_pred_num = 0
        true_standard_num = len(standard_result)
        for key in recommendation_result:
            #if len(recommendation_result[key]) != topk:
            #    print 'Not enough recommendation result for key %d!' % key
            #    print recommendation_result[key]
            #    sys.exit(1)
            for pid in recommendation_result[key]:
                if pid == standard_result[key]:
                    true_pred_num += 1
                sum_pred_num += 1
        precision = float(true_pred_num) / sum_pred_num
        recall = float(true_pred_num) / true_standard_num
        f_score = (2*precision*recall) / (precision+recall)
        return (precision, recall, f_score)


def loadResult(infile, choice, topk):
    if choice == 1:
        result = {}
        for i,line in enumerate(csv.reader(open(infile))):
            key = i
            pid = int(line[4])
            result[key] = pid
        return result
    elif choice == 2:
        result = defaultdict(list)
        for line in open(infile):
            parts = line.strip("\r\t\n").split("\t")
            key = int(parts[0])
            result[key] = map(int, parts[1].split(","))[:topk]
        return result
    else:
        print 'Invalid setting of result dataset format!'
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    parser.add_argument('-a', type=int, action='store',
            dest='algorithm_num', help='specify the algorithm which genenrate the \
                    recommendation results')
    parser.add_argument('-m', type=int, action='store',
            dest='eval_method', help='choose which method to evaluate the results')
    if len(sys.argv) != 7:
        print 'Command e.g.: python filterUserAndLocationByFreq.py -d 0(1,2)\
                -a 0(0,1,...) -m 0(1)'
        sys.exit(1)
    para = parser.parse_args()
    if para.data_num == 0:
        standard_result_file = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE1"]
    elif para.data_num == 1:
        standard_result_file = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE2"]
    elif para.data_num == 2:
        standard_result_file = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE3"]
    else:
        print 'Invalid choice of data set!'
        sys.exit(1)
    if para.algorithm_num == 0:
        recommendation_result_file = settings["ROOT_PATH"] + settings["POPULAR_SUBMISSION_PATH"]
    elif para.algorithm_num == 1:
        recommendation_result_file = settings["ROOT_PATH"] + settings["PER_POPULAR_SUBMISSION_PATH"]
    else:
        print 'Invalid choice of algorithm!'
        sys.exit(1)

    standard_result = loadResult(standard_result_file, 1, settings["TOPK"])
    recommendation_result = loadResult(recommendation_result_file, 2, settings["TOPK"])

    evaluation = Evaluation()
    if para.eval_method == 0:
        evaluation.evaluate(standard_result,
                            recommendation_result,
                            settings["F-score"],
                            settings["TOPK"])
    else:
        print 'Invalid choice of evalution method!'

if __name__ == "__main__":
    main()
