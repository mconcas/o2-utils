#! /usr/bin/bash -ex

if [[ -z "$1" ]];
 then
 evt=5
else
 evt=$1
fi

pushd ~/alice/data/000vtx-alien/
root -l -q run_vert_ca_its.C+\($evt\)
popd
