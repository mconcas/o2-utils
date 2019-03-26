#! /usr/bin/bash -ex

if [[ -z "$1" ]];
then
  mcvalid=
else
  mcvalid=true
  if [[ -z "$2" ]];
  then
    evt=
  else
    startfrom=",$2"
    if [[ -z "$3" ]];
    then
      endat=
    else
      endat=",$3"
    fi
  fi
fi

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd ~/alice/data/000vtx-alien/
cuda-memcheck root -l -q run_vert_gpu.C++\($mcvalid$startfrom$endat\)
echo test returned $?
popd
