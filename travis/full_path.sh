#!/bin/sh

TARGET_FILE=$1

while [ "$TARGET_FILE" != "" ]; do
    cd `dirname $TARGET_FILE`
    FILENAME=`basename $TARGET_FILE`
    TARGET_FILE=`readlink $FILENAME`
done

echo `pwd -P`/$FILENAME
