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
# Date: 2014/4/29                                                 #
# Calculate teh distance between consequtive checkin location.    #
#   Note: directly calculating distance by scanning all nodes     #
#         is too slow. Currently, we utilize QuadTree to speed    #
#         up the search process                                   #
###################################################################

import sys, csv, json, argparse, math, time
sys.path.append("../geopy-0.95.1/")
from geopy import distance
from quadtree import Quadtree

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())
distance.distance = distance.GreatCircleDistance


def findNearPlaceByQuadtree(loc_latlng, query, quadtree, dis_threshold):
    answer = []
    stack = []
    subtree = quadtree
    while True:
        if subtree['nodes'] == 4:
            for i in xrange(4):
                sontree = subtree['nodes'][i]
                if query[0] >= sontree["bounds"][0] and\
                        query[0] <= sontree["bounds"][2] and\
                        query[1] >= sontree["bounds"][1] and\
                        query[1] <= sontree["bounds"][3]:
                    stack.append([subtree, i])
                    subtree = sontree
                    break
        else:
            for pid in subtree['ids']:
                dis = distance.distance(loc_latlng[pid], query).miles
                if dis < dis_threshold:
                    answer.append(pid)
            break
    for ii in xrange(len(stack)-1, -1, -1):
        entry = stack[ii]
        sub_answer = []
        for i in xrange(4):
            if i != entry[1]:
                sub_answer += findNearPlaceByQuadtree(loc_latlng,
                                                     query,
                                                     entry[0]['nodes'][i],
                                                     dis_threshold)
        if len(sub_answer) == 0:
            break
        else:
            answer += sub_answer
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    if len(sys.argv) != 3:
        print 'Command e.g.: python findNearPlace.py -d 0(1,2)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        nearplace_outfile = settings["ROOT_PATH"] + settings["NEAR_PLACE_FILE1"]
    elif para.data_num == 1:
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        nearplace_outfile = settings["ROOT_PATH"] + settings["NEAR_PLACE_FILE2"]
    elif para.data_num == 2:
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_3"]
        nearplace_outfile = settings["ROOT_PATH"] + settings["NEAR_PLACE_FILE3"]
    else:
        print 'Invalid choice of data set'
        sys.exit(1)

    loc_latlng = {}
    try:
        for entry in csv.reader(open(location_infile, 'rU')):
            pid, lat, lng = int(entry[0]), float(entry[2]), float(entry[3])
            loc_latlng[pid] = (lat, lng)
    except:
        print entry
        sys.exit(1)

    # directly scanning all POIs to get answer, which is too slow
    '''writer = csv.writer(open(nearplace_outfile, "w"), lineterminator="\r\n")
    pids = loc_latlng.keys()
    for i in xrange(len(pids)):
        pid1 = pids[i]
        near_place = []
        for j in xrange(len(pids)):
            pid2 = pids[j]
            dis = distance.distance(loc_latlng[pid1], loc_latlng[pid2]).miles
            if dis < settings["DISTANCE_THRESHOLD"]:
                near_place.append(pid2)
        writer.writerow([pid1] + near_place)
        print i'''

    # quad tree
    index_extent = (-180, -90, 180, 90)
    index = Quadtree(index_extent)
    for pid in loc_latlng:
        index.add(pid, loc_latlng[pid])

    for pid in loc_latlng:
        start_time = time.clock()
        pid_set = findNearPlaceByQuadtree(loc_latlng,
                                          loc_latlng[pid],
                                          index.struct(),
                                          settings["DISTANCE_THRESHOLD"])
        end_time = time.clock()
        print "Time Cost: %f(s)" % (end_time-start_time)
        raw_input()
        print len(pid_set)
        raw_input()


if __name__ == "__main__":
    main()
