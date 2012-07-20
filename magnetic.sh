#!/bin/bash

scp magnetic.cpp lyubimov@cluster.dvfu.ru:~/magnetic/
ssh lyubimov@cluster.dvfu.ru "cd ~/magnetic/; mpic++ -o magnetic.o  magnetic.cpp; mpiexec -np 8 ./magnetic.o"
sftp -r lyubimov@cluster.dvfu.ru:magnetic/initial.txt  ./data/
sftp -r lyubimov@cluster.dvfu.ru:magnetic/final.txt  ./data/

./generate.sh