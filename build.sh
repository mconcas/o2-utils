#! /usr/bin/bash

export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)

pushd /home/alidock/alice/sw/BUILD/O2-latest-$(git --git-dir=/home/alidock/alice/O2/.git rev-parse --abbrev-ref HEAD)/O2 &&
ninja -j12 install
popd
