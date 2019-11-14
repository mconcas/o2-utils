#!/bin/bash
export ALIBUILD_WORK_DIR="/home/alidock/alice/sw"
eval $(alienv load O2/latest 2> /dev/null)
env
cd /home/alidock/alice/data/000vtx-alien
exec root.exe -l -q run_vert_gpu.C++\(14\)
#root -l -q /home/alidock/alice/data/000vtx-alien/run_vert_gpu.C++\(5\)
