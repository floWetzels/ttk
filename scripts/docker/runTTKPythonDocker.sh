#!/bin/bash

PARAVIEW_VERSION=5.9.1
TTK_VERSION=0.9.9

#if [ -z "$1" ]; then
#    PARAVIEW_VERSION="$1"
#fi
#if [ -z "$2" ]; then
#    PARAVIEW_VERSION="$2"
#fi

docker run -it --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-${TTK_VERSION} pvpython ${@:1}
