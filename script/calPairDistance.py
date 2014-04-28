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
# Date: 2014/4/28                                                 #
# Calculate teh distance between consequtive checkin location     #
###################################################################

import sys, csv, json, argparse, math, pylab
sys.path.append("../geopy-0.95.1/")
from geopy import distance
import numpy as np

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())
distance.distance = distance.GreatCircleDistance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    if len(sys.argv) != 3:
        print 'Command e.g.: python calPairDistance.py -d 0(1,2)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE1"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
    elif para.data_num == 1:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE2"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
    elif para.data_num == 2:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE3"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_3"]

    loc_latlng = {}
    if para.data_num == 0 or para.data_num == 1:
        pass
    elif para.data_num == 2:
        tag = False
        for line in open(location_infile):
            if not tag:
                tag = True
                continue
            entry = line.strip("\r\t\n").replace("\,", " ").split(",")
            locid, lat, lng = int(entry[0]), float(entry[2]), float(entry[3])
            loc_latlng[locid] = [lat, lng]

    dis_cnt = {}
    tag = False
    for line in csv.reader(open(checkin_infile)):
        if not tag:
            tag = True
            continue
        locid1, locid2 = int(line[1]), int(line[4])
        dis = int(math.ceil(distance.distance(loc_latlng[locid1], loc_latlng[locid2]).miles))
        if dis not in dis_cnt:
            dis_cnt[dis] = 1
        else:
            dis_cnt[dis] += 1

    dis_cnt = sorted(dis_cnt.items(), key=lambda x:x[0])
    keys = [entry[0] for entry in dis_cnt]
    vals = [entry[1] for entry in dis_cnt]
    width = 0.2
    pylab.xticks(np.array(keys[:50])+width/2.0, keys, rotation=45)
    pylab.bar(keys[:50], vals[:50], width, color='r')
    pylab.show()
    print keys[50], vals[50]
    dis_cnt = sorted(dis_cnt, key=lambda x:x[1])
    print dis_cnt[-1][0], dis_cnt[-1][1]

if __name__ == "__main__":
    main()
