#!/usr/bin/env bash

# This is a script transfers files from s3 to local directory

while [[ $# -gt 0 ]]
do
key="$1"

case ${key} in
    -f|--file_bucket)
    BUCKETPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -o|-output)
    SAVEDIR="$2"
    shift # past argument
    shift # past value
    ;;

esac
done


# echo BUCKETPATH  = "${BUCKETPATH}"
# echo SAVEDIR = "${SAVEDIR}"

echo `aws s3 cp s3://${BUCKETPATH} ${SAVEDIR}`

