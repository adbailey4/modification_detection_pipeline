#!/usr/bin/env bash

# This is a script transfers files to s3 from local directory

while [[ $# -gt 0 ]]
do
key="$1"

case ${key} in
    -f|--file)
    FILEPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -b|--bucket)
    BUCKET="$2"
    shift # past argument
    shift # past value
    ;;

esac
done


echo FILEPATH  = "${FILEPATH}"
echo BUCKET = "${BUCKET}"

echo `aws s3 cp ${FILEPATH} s3://${BUCKET}`

