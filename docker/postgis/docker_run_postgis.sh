#! /bin/bash
docker run --name tin_db -v $(pwd):/tmp -w /tmp -p 5501:5432 -d mdillon/postgis:latest
sleep 5
docker exec tin_db bash /tmp/init_db.sh
