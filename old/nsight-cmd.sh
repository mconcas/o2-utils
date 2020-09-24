#!/bin/bash
export ALIBUILD_WORK_DIR="/home/mconcas/alice/sw"
eval $(alienv load O2/latest-pr-NO-debug-cuda10.2-o2 2> /dev/null)
env
cd /data1/ruben/pbpbVtx
# exec root.exe -l -q run_primary_vertexer_ITS.C++\(GPUDataTypes::DeviceType::CUDA\)
exec root.exe -l -q run_primary_vertexer_ITS.C++\(GPUDataTypes::DeviceType::HIP\)
#root -l -q /home/alidock/alice/data/000vtx-alien/run_vert_gpu.C++\(5\)
