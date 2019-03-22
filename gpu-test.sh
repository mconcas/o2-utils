#! /usr/bin/bash -ex

if [[ -z "$1" ]];
 then
 evt=
else
 evt=$1
fi

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd ~/alice/data/000vtx-alien/
cuda-memcheck root -l -q run_vert_gpu.C++\($evt\)
echo test returned $?
popd
