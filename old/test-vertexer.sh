#! /usr/bin/bash -e

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd /mnt/PythiaPbPb/MB/001/
root -l -q /home/alidock/alice/O2/macro/run_primary_vertexer_ITS.C++\("${1:-false},${2:-false},${3:--1},${4:-1}"\)
echo test returned $?
popd
