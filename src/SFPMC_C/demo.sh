#/bin/sh

make
./run -d 2 -r True -b False

#valgrind --tool=memcheck --leak-check=full --show-reachable=yes --error-limit=no --log-file=valgrind.log ./demo.sh
