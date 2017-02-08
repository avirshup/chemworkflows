#!/bin/bash

version=$1

if [ -z $version ]; then
   echo "No version passed ..."
   exit 127
fi

docker-make --all --repo docker.io/avirshup/mst: --tag $version

if [ $? -ne 0 ]; then
   echo "Build failed."
   exit 1
fi

./test.sh $version

