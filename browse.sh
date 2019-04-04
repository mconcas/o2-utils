#! /usr/bin/bash -ex

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

# pushd ~/alice/data/000vtx-alien/
root -l /home/alidock/util/VertexerDataBrowser.C++
echo test returned $?
# popd
