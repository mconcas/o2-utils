#! /usr/bin/bash -ex

if [[ -z "$1" ]];
 then
 evt=5
else
 evt=$1
fi

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

# pushd ~/alice/data/000vtx-alien/
root -l /home/alidock/util/VertexerDataBrowser.C++
echo test returned $?
# popd
