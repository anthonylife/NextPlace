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
# Generate consequtive checkin pair for training and testing      #
###################################################################

import sys, csv, json, argparse, datetime

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())
dt = datetime.datetime.now()


def genPair(entries):
    pairs = []
    for i in xrange(len(entries)-1):
        time1 = dt.strptime(entries[i][3], '%Y-%m-%d %H:%M:%S')
        time2 = dt.strptime(entries[i+1][3],'%Y-%m-%d %H:%M:%S')
        time_diff = time2-time1
        if time_diff.seconds < settings["TIME_SECONDS_THRESHOLD"]:
            uid = entries[i][0]
            pid1 = entries[i][1]
            pid2 = entries[i+1][1]
            day1 = time1.weekday()
            day2 = time2.weekday()
            hour1 = time1.hour
            hour2 = time2.hour
            if pid1 != pid2:
                pairs.append([uid, pid1, day1, hour1, pid2, day2, hour2, entries[i][3]])
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')

    if len(sys.argv) != 3:
        print 'Command e.g.: python genConsequtiveCheckinPair.py -d 0(1,2)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        checkin_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE1_1"]
        checkinpair_outfile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE1"]
    elif para.data_num == 1:
        checkin_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE2_1"]
        checkinpair_outfile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE2"]
    elif para.data_num == 2:
        checkin_infile = settings["ROOT_PATH"] + settings["SRC_DATA_FILE3_1"]
        checkinpair_outfile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE3"]

        last_userid = "-1"
        cache_data = []
        writer = csv.writer(open(checkinpair_outfile, "w"), lineterminator="\n")
        writer.writerow(["uid", "pid1", "day1", "hour1", "pid2", "day2", "hour2"])
        for entry in csv.reader(open(checkin_infile)):
            now_userid = entry[1]
            if last_userid == "-1" or now_userid == last_userid:
                cache_data.append(entry[1:])
            else:
                cache_data = sorted(cache_data, key=lambda x:x[3])
                output_result = genPair(cache_data)
                writer.writerows(output_result)
                cache_data = []
                cache_data.append(entry[1:])
            last_userid = now_userid
    else:
        print 'Choice of file invalid!'
        sys.exit(1)


if __name__ == "__main__":
    main()
