#!/bin/bash

### даунлоад через SSH (c 10.1.2.3 на локальную машину)
### sftp -r backuper@10.1.2.3:/BACKUP/*  /backup/current/

### аплоад через SSH (c локальной машины на 10.1.2.3)
scp /home/zlo/Desktop/* backuper@10.1.2.3:/backup/zlo/