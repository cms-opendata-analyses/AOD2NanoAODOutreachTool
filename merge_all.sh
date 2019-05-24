#!/bin/bash

BASE_PATH=/path/to/opendata/files

for FOLDER in $(ls $BASE_PATH | grep -v .root)
do
    ./merge_jobs.py ${BASE_PATH}/${FOLDER} &
done

wait
