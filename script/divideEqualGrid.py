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
# Divide the whole earth into equally sized grids.                #
#   Note: x-->longitude, y-->latitude.                            #
###################################################################

import sys, csv, json, argparse, math
from geopy import distance

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())
distance.distance = distance.GreatCircleDistance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    if len(sys.argv) != 3:
        print 'Command e.g.: python findNearPlace.py -d 0(1,2)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE1"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        grid_place_outfile = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE1"]
    elif para.data_num == 1:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE2"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        grid_place_outfile = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE2"]
    elif para.data_num == 2:
        checkin_infile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE3"]
        location_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_3"]
        grid_place_outfile = settings["ROOT_PATH"] + settings["GRID_PLACE_FILE3"]
    else:
        print 'Invalid choice of data set'
        sys.exit(1)

    locid_set = set([])
    for i, entry in enumerate(csv.reader(open(checkin_infile))):
        if i==0:
            continue
        pid1, pid2 = int(entry[1]), int(entry[4])
        locid_set.add(pid1)
        locid_set.add(pid2)

    loc_latlng = {}
    try:
        for entry in csv.reader(open(location_infile, 'rU')):
            pid, lat, lng = int(entry[0]), float(entry[2]), float(entry[3])
            if pid in locid_set:
                loc_latlng[pid] = (lat, lng)
    except:
        print entry
        sys.exit(1)
    locid_set.clear()

    # Grids construction
    index_extent = (-90, -180, 90, 180)
    ndimx = int((index_extent[3]-index_extent[1])/settings["GRID_LNG"])
    ndimy = int((index_extent[2]-index_extent[0])/settings["GRID_LAT"])
    grids = {}
    for i in xrange(ndimx):
        grids[i] = {}
        for j in xrange(ndimy):
            grids[i][j] = []
    for locid in loc_latlng:
        x_idx = int(math.floor((loc_latlng[locid][1]+180)/settings["GRID_LNG"]))
        y_idx = int(math.floor((loc_latlng[locid][0]+90)/settings["GRID_LAT"]))
        if loc_latlng[locid][0] == 90:
            y_idx = ndimy-1
        if loc_latlng[locid][1] == 180:
            x_idx = ndimx-1
        if x_idx in grids and y_idx in grids[x_idx]:
            grids[x_idx][y_idx].append(locid)
        else:
            print locid, loc_latlng[locid][1], loc_latlng[locid][0], x_idx, y_idx
            print 'Incorrect construction of grids.'
            sys.exit(1)

    writer = csv.writer(open(grid_place_outfile, "w"), lineterminator="\n")
    for x_idx in grids:
        for y_idx in grids[x_idx]:
            writer.writerow([x_idx, y_idx] + list(grids[x_idx][y_idx]))


if __name__ == "__main__":
    main()
