#! /usr/bin/bash -e

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd ~/alice/data/000vtx-alien/
root -l -q run_primary_vertexer_ITS.C++\("${1:-false},${2:-false},${3:--1},${4:-1}"\)
echo test returned $?
popd
