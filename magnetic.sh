#!/bin/bash

### даунлоад через SSH (c 10.1.2.3 на локальную машину)
### sftp -r backuper@10.1.2.3:/BACKUP/*  /backup/current/

### аплоад через SSH (c локальной машины на 10.1.2.3)
scp magnetic.cpp lyubimov@cluster.dvfu.ru:~/magnetic/
ssh lyubimov@cluster.dvfu.ru "cd ~/magnetic/; mpic++ -o magnetic.o  magnetic.cpp; mpiexec -np 4 ./magnetic.o"
sftp -r lyubimov@cluster.dvfu.ru:~/magnetic/initial.txt  ./data/
sftp -r lyubimov@cluster.dvfu.ru:~/magnetic/final.txt  ./data/