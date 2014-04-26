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
# Date: 2014/4/25                                                 #
# Count the number of pairs each user have                        #
###################################################################

import sys, csv, json, argparse, pylab
import numpy as np

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    parser.add_argument('-t', type=int, action='store',
            dest='data_type', help='choose which format of data')

    if len(sys.argv) != 5:
        print 'Command e.g.: python genConsequtiveCheckinPair.py -d 0(1,2) -t 0(1)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        if para.data_type == 0:
            checkinpair_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        elif para.data_type == 1:
            checkinpair_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE1"]
    elif para.data_num == 1:
        if para.data_type == 0:
            checkinpair_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        elif para.data_type == 1:
            checkinpair_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE2"]
    elif para.data_num == 2:
        if para.data_type == 0:
            checkinpair_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_1"]
        elif para.data_type == 1:
            checkinpair_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE3"]

    uid_cnt = {}
    for line in csv.reader(open(checkinpair_infile)):
        if para.data_type == 1:
            uid = line[0]
        elif para.data_type == 0:
            uid = line[1]
        if uid not in uid_cnt:
            uid_cnt[uid] = 1
        else:
            uid_cnt[uid] += 1

    cnt_num = {}
    for key in uid_cnt:
        num = uid_cnt[key]
        if num not in cnt_num:
            cnt_num[num] = 1
        else:
            cnt_num[num] += 1

    cnt_num = sorted(cnt_num.items(), key=lambda x:x[0])
    keys = [entry[0] for entry in cnt_num]
    vals = [entry[1] for entry in cnt_num]
    width = 0.2
    pylab.xticks(np.array(keys[:50])+width/2.0, keys, rotation=45)
    pylab.bar(keys[:50], vals[:50], width, color='r')
    pylab.show()
    print keys[50], vals[50]


if __name__ == "__main__":
    main()

