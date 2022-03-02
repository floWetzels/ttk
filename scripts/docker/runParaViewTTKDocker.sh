#!/bin/bash

PARAVIEW_PATH=$1
#PARAVIEW_VERSION=5.10.0

if [ -z "${PARAVIEW_PATH}" ];

then
    if ! [ -x "$(command -v paraview)" ]; then
        echo "No ParaView installation found. Please provide path to ParaView as first argument."
        echo "Usage:"
        echo "  $0 <Path to ParaView binary>"
        exit 0
    fi
    PV_VERSION_OUTPUT = $(paraview --version)
    PARAVIEW_VERSION = "${PV_VERSION_OUTPUT:17}"
    echo "$PARAVIEW_VERSION"
    
else
    PV_VERSION_OUTPUT=$($PARAVIEW_PATH --version)
    PARAVIEW_VERSION="${PV_VERSION_OUTPUT:17}"
    echo "$PARAVIEW_VERSION"

fi

DOCKER_ID=`docker run -d --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-master`
if ! [ -n "${DOCKER_ID}" ]; then
    DOCKER_ID=`docker run -d --rm -p 11111:11111 -v "${HOME}:/home/${USER}/" --user ${UID} topologytoolkit/ttk:${PARAVIEW_VERSION}-dev`
    if ! [ -n "${DOCKER_ID}" ]; then
        echo "No docker image for ParaView version ${PARAVIEW_VERSION} and TTK-master or TTK-dev found."
        echo "Please use a different ParaView version."
        exit 0
    fi
fi

${PARAVIEW_PATH} --server-url=cs://127.0.0.1:11111 ${@:2}
docker kill ${DOCKER_ID} &> /dev/null
