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
# Filter users and locations by the specified frequency           #
###################################################################

import sys, csv, json, argparse
from collections import defaultdict

with open("../SETTINGS.json") as fp:
    settings = json.loads(fp.read())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=int, action='store',
            dest='data_num', help='choose which data set to use')
    parser.add_argument('-m', type=int, action='store',
            dest='filter_method', help='choose which method to filter data')
    if len(sys.argv) != 5:
        print 'Command e.g.: python genConsequtiveCheckinPair.py -d 0(1,2) -m 0(1)'
        sys.exit(1)

    para = parser.parse_args()
    if para.data_num == 0:
        checkin_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE1"]
        checkin_outfile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE1"]
    elif para.data_num == 1:
        checkin_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE2"]
        checkin_outfile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PATR_FILE2"]
    elif para.data_num == 2:
        checkin_infile = settings["ROOT_PATH"] + settings["CHECKIN_PAIR_FILE3"]
        checkin_outfile = settings["ROOT_PATH"] + settings["FILTER_CHECKIN_PAIR_FILE3"]

    uid_set = set([])
    pid_set = set([])
    if para.filter_method == 0:
        pid_uid = defaultdict(set)
        uid_pid = defaultdict(set)
        tag = False
        for line in csv.reader(open(checkin_infile)):
            if not tag:
                tag = True
                continue
            entry = map(int, line[:-1])
            uid, pid1, pid2 = entry[0], entry[1], entry[4]
            pid_uid[pid1].add(uid)
            pid_uid[pid2].add(uid)
            uid_pid[uid].add(pid1)
            uid_pid[uid].add(pid2)
        removed_pid = set([])
        removed_uid = set([])
        while True:
            removed_pid.clear()
            for pid in pid_uid:
                pid_uid[pid] = pid_uid[pid] - removed_uid
                if len(pid_uid[pid]) < settings["FILTER_LOCATION_VISIT_NUM"]:
                    removed_pid.add(pid)
            for pid in removed_pid:
                pid_uid.pop(pid)
            removed_uid.clear()
            for uid in uid_pid:
                uid_pid[uid] = uid_pid[uid]-removed_pid
                if len(uid_pid[uid]) < settings["FILTER_USER_VISIT_NUM"]:
                    removed_uid.add(uid)
            for uid in removed_uid:
                uid_pid.pop(uid)
            if len(removed_uid) == 0:
                uid_set = set(uid_pid.keys())
                pid_set = set(pid_uid.keys())
                uid_pid = None
                pid_uid = None
                break
    elif para.filter_method == 1:
        data = [entry for entry in csv.reader(open(checkin_infile))]
        data = [map(int, entry[:-1]) for entry in data[1:]]
        pid_record = defaultdict(set)
        uid_record = defaultdict(set)
        for i, entry in enumerate(data):
            uid, pid1, pid2 = entry[0], entry[1], entry[4]
            uid_record[uid].add(i)
            pid_record[pid1].add(i)
            pid_record[pid2].add(i)
        removed_record = set([])
        removed_pid = set([])
        removed_uid = set([])
        while True:
            for idx in removed_record:
                uid, pid1, pid2 = data[idx][0], data[idx][1], data[idx][4]
                pid_record[pid1] = pid_record[pid1] - set([idx])
                pid_record[pid2] = pid_record[pid2] - set([idx])
            removed_record.clear()
            for pid in pid_record:
                if len(pid_record[pid]) < settings["FILTER_LOCATION_RECORD_NUM"]:
                    removed_pid.add(pid)
                    for idx in pid_record[pid]:
                        removed_record.add(idx)
            for pid in removed_pid:
                pid_record.pop(pid)
            removed_pid.clear()
            for uid in uid_record:
                if len(uid_record[uid]) < settings["FILTER_USER_RECORD_NUM"]:
                    removed_uid.add(uid)
                    for idx in uid_record[uid]:
                        removed_record.add(idx)
            for uid in removed_uid:
                uid_record.pop(uid)
            removed_uid.clear()
            print "Removed Record Number: %d" % len(removed_record)
            if len(removed_record) == 0:
                uid_set = set(uid_record.keys())
                pid_set = set(pid_record.keys())
                data = None
                uid_record = None
                pid_record = None
                break

    tag = False
    with open(checkin_outfile, "w") as wfp:
        writer = csv.writer(wfp, lineterminator="\n")
        for entry in csv.reader(open(checkin_infile)):
            if not tag:
                tag = True
                writer.writerow(entry)
            else:
                uid, pid1, pid2 = map(int, [entry[0], entry[1], entry[4]])
                if uid in uid_set and pid1 in pid_set and pid2 in pid_set:
                    writer.writerow(entry)

if __name__ == "__main__":
    main()

