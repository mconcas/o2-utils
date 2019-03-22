#! /usr/bin/bash -ex

if [[ -z "$1" ]];
 then
 evt=
else
 evt=$1
fi

pushd ~/alice/data/000vtx-alien/
cuda-memcheck root -l -q run_vert_gpu.C++\($evt\)
echo test returned $?
popd
