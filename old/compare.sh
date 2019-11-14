#! /usr/bin/bash -ex

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd ~/alice/data/000vtx-alien/
root -l -q compare.C++\(\)
echo test returned $?
popd
