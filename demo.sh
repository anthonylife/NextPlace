#!/bin/sh

python setenv.py

echo 'Data preprocessing...'
cd script/
# generate consequtive checkin pair
python genConsequtiveCheckinPair.py -d 2
# filter checkin pair, user and location by specified frequency
python filterUserAndLocationByFreq.py -d 2 -m 1
# segment data for training and evaluation
python segCheckinPair.py -d 2 -m 0
# divide map space into equal-sized grids
python divideEqualGrid.py -d 2


echo 'Runing diverse next location prediction algorithm...'
# running global popularity based method
cd ../src/Popular/
python globalPopular.py -d 2
cd ../../script/
python evalution.py -d 2 -a 0 -m 0

# running personal popularity based method
cd ../src/PerPopular/
python personalPopular.py -d 2
cd ../../script/
python evalution.py -d 2 -a 1 -m 0

