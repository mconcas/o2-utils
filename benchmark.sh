#! /usr/bin/bash -ex

if [[ -z "$1" ]];
 then
 evt=5
else
 evt=$1
fi

pushd ~/alice/data/000vtx-alien/
nvprof --profile-child-processes root -l -q run_vert_gpu.C++\($evt\)
echo test returned $?
popd
