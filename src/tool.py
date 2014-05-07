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
# Providing all tool functions for each algorithms                #
###################################################################

import sys, csv, json, math, random
from collections import defaultdict

with open("../../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


def loadGridInfo(infile):
    grids = defaultdict(dict)
    for entry in csv.reader(open(infile)):
        xidx,yidx,pids = int(entry[0]), int(entry[1]), map(int, entry[2:])
        grids[xidx][yidx] = pids
    return grids


def loadPoiInfo(infile, data_num):
    pois = {}
    if data_num == 0 or data_num == 1:
        for entry in open(infile):
            parts = entry.strip("\r\t\n").split("\t")
            pid, lat, lng = int(parts[4]), float(parts[2]), float(parts[3])
            pois[pid] = (lat, lng)
    elif data_num == 2:
        for entry in csv.reader(open(infile)):
            pid, lat, lng = int(entry[0]), float(entry[2]), float(entry[3])
            pois[pid] = (lat, lng)
    return pois


def getGridIdxForPOI(latlng, ndimx, ndimy):
    x_idx = int(math.floor((latlng[1]+180)/settings["GRID_LNG"]))
    y_idx = int(math.floor((latlng[0]+90)/settings["GRID_LAT"]))
    if latlng[0] == 90:
        y_idx = ndimy-1
    if latlng[1] == 180:
        x_idx = ndimx-1
    return (x_idx, y_idx)


def checkBoundary(grididx, ndimx, ndimy):
    revised_grididx = [0, 0]
    revised_grididx[0] = grididx[0]
    revised_grididx[1] = grididx[1]
    if grididx[0] == -1:
        revised_grididx[0] = ndimx-1
    elif grididx[0] == ndimx:
        revised_grididx[0] = 0
    if grididx[1] == -1:
        revised_grididx[1] = ndimy-1
    elif grididx[1] == ndimy:
        revised_grididx[1] = 0
    return revised_grididx


def getNearGrids(grididx, ndimx, ndimy):
    if grididx[0] < 0 or grididx[0] >= ndimx or grididx[1] < 0 or grididx[1] >= ndimy:
        print 'Invalid of grididx!'
        sys.exit(1)
    top_left = checkBoundary((grididx[0]-1, grididx[1]+1), ndimx, ndimy)
    above = checkBoundary((grididx[0], grididx[1]+1), ndimx, ndimy)
    top_right = checkBoundary((grididx[0]+1, grididx[1]+1), ndimx, ndimy)
    left = checkBoundary((grididx[0]-1, grididx[1]), ndimx, ndimy)
    right = checkBoundary((grididx[0]+1, grididx[1]), ndimx, ndimy)
    bottom_left = checkBoundary((grididx[0]-1, grididx[1]-1), ndimx, ndimy)
    below = checkBoundary((grididx[0], grididx[1]-1), ndimx, ndimy)
    bottom_right = checkBoundary((grididx[0]+1, grididx[1]-1), ndimx, ndimy)
    near_grids = [top_left, above, top_right, left, right, bottom_left, below, bottom_right]
    return near_grids


def getNearGridsForPOI(latlng, ndimx, ndimy, tag):
    grididx = getGridIdxForPOI(latlng, ndimx, ndimy)
    near_grids = getNearGrids(grididx, ndimx, ndimy)
    if tag:
        near_grids.append(grididx)
        return near_grids
    else:
        return near_grids


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


def getMulMapId(infile):
    user_ids = {}
    ruser_ids = {}
    poi_ids = {}
    rpoi_ids = {}
    u_cnt = 0
    p_cnt = 0
    for entry in csv.reader(open(infile)):
        uid, pid1, pid2 = int(entry[0]), int(entry[1]), int(entry[4])
        if uid not in user_ids:
            user_ids[uid] = u_cnt
            ruser_ids[u_cnt] = uid
            u_cnt += 1
        if pid1 not in poi_ids:
            poi_ids[pid1] = p_cnt
            rpoi_ids[p_cnt] = pid1
            p_cnt += 1
        if pid2 not in poi_ids:
            poi_ids[pid2] = p_cnt
            rpoi_ids[p_cnt] = pid2
            p_cnt += 1
    return user_ids, ruser_ids, poi_ids, rpoi_ids


def rZero(k):
    return [0.0 for i in range(k)]


def rGaussian(k):
    factor = [random.normalvariate(0, 0.01) for i in xrange(k)]
    for i in xrange(len(factor)):
        if factor[i] > 1:
            factor[i] = 1
        elif factor[i] < -1:
            factor[i] = -1
    return factor


def logitLoss(pos_score, neg_score):
    return (1-1.0/(1+math.exp(-(pos_score-neg_score))))


def test():
    # Testing function checkBoundary()
    print "[0, 3599]"
    print getNearGrids([0, 3599], 7200, 3600)


if __name__ == "__main__":
    test()
