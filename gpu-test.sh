#! /usr/bin/bash -ex

if [[ -z "$1" ]];
  then
    evt=
  else
    startfrom=$1
fi

if [[ -z "$2" ]];
  then
    endat=
  else
    endat=",$2"
fi

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd ~/alice/data/000vtx-alien/
cuda-memcheck root -l -q run_vert_gpu.C++\($startfrom$endat\)
echo test returned $?
popd
