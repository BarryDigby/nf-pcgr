#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers
# Mimicks Maxime Garcia's implementation of VEP containers

build_push() {
    PCGR_VERSION=$1
    GENOME=$2
    BUNDLE_VERSION=$3

    docker build \
        . \
        -t barryd237/pcgr:${PCGR_VERSION}.${GENOME} \
        --build-arg PCGR_VERSION=${PCGR_VERSION} \
        --build-arg GENOME=${GENOME} \
        --build-arg BUNDLE_VERSION=${BUNDLE_VERSION}

    docker push barryd237/pcgr:${PCGR_VERSION}.${GENOME}
}

build_push "1.0.3"    "grch37"           "20220203"
build_push "1.0.3"    "grch38"           "20220203"
