#! /usr/bin/bash -e

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

# pushd ~/alice/data/000vtx-alien/
root -l /home/alidock/util/VertexerDataBrowser.C++\("\"/home/alidock/alice/data/000vtx-alien/vertexer_${1:?choose either serial or gpu}_data.root"\"\)
echo test returned $?
# popd
