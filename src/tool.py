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

import sys, csv, json, math
from collections import defaultdict

with open("../../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


def loadGridInfo(infile):
    grids = defaultdict(dict)
    for entry in csv.reader(infile):
        xidx,yidx,pids = int(entry[0]), int(entry[1]), map(int, entry[2:])
        grids[xidx][yidx] = pids


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
    revised_grididx = [0.0, 0.0]
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
    top_left = checkBoundary((grididx[0]-1, grididx[1]+1), ndimx, ndimy)
    above = checkBoundary((grididx[0], grididx[1]+1), ndimx, ndimy)
    top_right = checkBoundary((grididx[0]+1, grididx[1]+1), ndimx, ndimy)
    left = checkBoundary((grididx[0]-1, grididx[1]), ndimx, ndimy)
    right = checkBoundary((grididx[0]+1, grididx[1]), ndimx, ndimy)
    bottom_left = checkBoundary((grididx[0]-1, grididx[1]-1), ndimx, ndimy)
    below = checkBoundary((grididx[0], grididx[1]-1), ndimx, ndimy)
    bottom_right = checkBoundary((grididx[0]+1, grididx[1]+1), ndimx, ndimy)

    near_grids = [top_left, above, top_right, left, right, bottom_left, below, bottom_right]
    return near_grids


def getNearGridsForPOI(latlng, ndimx, ndimy, tag):
    grididx = getGridIdxForPOI(latlng, ndimx, ndimy)
    near_grids = getNearGrids(grididx, ndimx, ndimy)
    if tag:
        return near_grids.append(grididx)
    else:
        return near_grids
