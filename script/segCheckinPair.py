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
# Date: 2014/4/26                                                 #
# Segment the whole dataset into training and test dataset.       #
# Currently, we plan to support 3 segmentation method:            #
#   1.segment the dataset into training and test set randomly     #
#   2.segment the dataset into training and test set based on time#
#   3.cross validation                                            #
# Note: same as previous works, we need to ensure that at least   #
#   one item and one user occurred in training data.              #
###################################################################

import sys, csv, json, argparse, random, math

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    parser.add_argument('-m', type=int, action='store',
            dest='seg_method', help='choose which method to segment data')
    if len(sys.argv) != 5:
        print 'Command e.g.: python segCheckinPair.py -d 0(1,2) -m 0(1,2)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE1"]
        train_outfile = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE1"]
        vali_outfile = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE1"]
        test_outfile = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE1"]
        cv_outfile = settings["ROOT_PATH"] + settings["CV_PAIR_FILE1"]
    elif para.data_num == 1:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE2"]
        train_outfile = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE2"]
        vali_outfile = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE2"]
        test_outfile = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE2"]
        cv_outfile = settings["ROOT_PATH"] + settings["CV_PAIR_FILE2"]
    elif para.data_num == 2:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE3"]
        train_outfile = settings["ROOT_PATH"] + settings["TRAIN_PAIR_FILE3"]
        vali_outfile = settings["ROOT_PATH"] + settings["VALI_PAIR_FILE3"]
        test_outfile = settings["ROOT_PATH"] + settings["TEST_PAIR_FILE3"]
        cv_outfile = settings["ROOT_PATH"] + settings["CV_PAIR_FILE3"]

    if para.seg_method == 0:
        tr_wfd = csv.writer(open(train_outfile, "w"), lineterminator="\n")
        va_wfd = csv.writer(open(vali_outfile, "w"), lineterminator="\n")
        te_wfd = csv.writer(open(test_outfile, "w"), lineterminator="\n")
        tr_uid = set([])
        tr_pid = set([])
        va_entry = []
        te_entry = []
        tag = False
        for entry in csv.reader(open(checkin_infile)):
            if not tag:
                tag=True
                continue
            s_ratio = random.random()
            uid, pid1, pid2 = entry[0], entry[1], entry[4]
            if s_ratio < settings["TRAIN_RATIO"]:
                tr_wfd.writerow(entry)
                tr_uid.add(uid)
                tr_pid.add(pid1)
                tr_pid.add(pid2)
            elif s_ratio < settings["TRAIN_RATIO"] + settings["VALI_RATIO"]:
                #va_wfd.writerow(entry)
                va_entry.append(entry)
            else:
                #te_wfd.writerow(entry)
                te_entry.append(entry)
        for entry in va_entry:
            uid, pid1, pid2 = entry[0], entry[1], entry[4]
            if uid not in tr_uid:
                tr_wfd.writerow(entry)
                tr_uid.add(uid)
            elif pid1 not in tr_pid:
                tr_wfd.writerow(entry)
                tr_pid.add(pid1)
            elif pid2 not in tr_pid:
                tr_wfd.writerow(entry)
                tr_pid.add(pid2)
            else:
                va_wfd.writerow(entry)
        del va_entry
        for entry in te_entry:
            uid, pid1, pid2 = entry[0], entry[1], entry[4]
            if uid not in tr_uid:
                tr_wfd.writerow(entry)
                tr_uid.add(uid)
            elif pid1 not in tr_pid:
                tr_wfd.writerow(entry)
                tr_pid.add(pid1)
            elif pid2 not in tr_pid:
                tr_wfd.writerow(entry)
                tr_pid.add(pid2)
            else:
                te_wfd.writerow(entry)
        del te_entry


    elif para.seg_method == 1:
        tr_wfd = csv.writer(open(train_outfile, "w"), lineterminator="\n")
        va_wfd = csv.writer(open(vali_outfile, "w"), lineterminator="\n")
        te_wfd = csv.writer(open(test_outfile, "w"), lineterminator="\n")
        dates = []
        tag = False
        for entry in csv.reader(open(checkin_infile)):
            if not tag:
                tag = True
                continue
            dates.append(entry[7])
        dates = sorted(dates)
        tr_date = dates[int(math.floor(len(dates)*settings["TRAIN_RATIO"]))]
        va_date = dates[int(math.floor(len(dates)*(settings["VALI_RATIO"])\
                +settings["TRAIN_RATIO"]))]
        for entry in csv.reader(open(checkin_infile)):
            if entry[7] < tr_date:
                tr_wfd.writerow(entry)
            elif entry[7] < va_date:
                va_wfd.writerow(entry)
            else:
                te_wfd.writerow(entry)

    elif para.seg_method == 2:
        cv_wfd = csv.writer(open(cv_outfile, "w"), lineterminator="\n")
        cv_data = [[] for i in xrange(settings["CV_NUM"])]
        tag = False
        for entry in csv.reader(open(checkin_infile)):
            if not tag:
                tag = True
                continue
            r_int = random.randint(0, settings["CV_NUM"]-1)
            cv_data[r_int].append(entry)
        for one_cv in cv_data:
            cv_wfd.writerow(len(one_cv))
            cv_wfd.writerows(one_cv)

    else:
        print 'Invalid choice of segmentation method.'
        sys.exit(1)

if __name__ == "__main__":
    main()

