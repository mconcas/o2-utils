#!/bin/bash
# cd /data1/ruben/pbpbVtx

# root -l -b -q run_primary_vertexer_ITS.C++\(0.005,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_0005.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.010,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_001.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.020,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_002.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.030,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_003.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.040,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_004.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.050,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_005.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.060,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_006.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.070,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_007.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.080,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_008.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.090,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_009.root
# root -l -b -q run_primary_vertexer_ITS.C++\(0.100,GPUDataTypes::DeviceType::CPU,false,20,1\) && mv dbg_ITSVertexerCPU.root phiCutVariationData/dbg_ITSVertexerCPU_PhiCut_010.root

# cd -

cd /data1/ruben/pbpbVtx

# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.001,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_001.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.002,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_002.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.003,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_003.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.004,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_004.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.005,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_005.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.006,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_006.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.007,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_007.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.008,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_008.root
# # root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.009,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_009.root
root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.010,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_010.root
root -l -b -q run_primary_vertexer_ITS.C++\(0.100,0.011,GPUDataTypes::DeviceType::CPU,false\) && mv dbg_ITSVertexerCPU.root tanLambdaCutVariationData/dbg_ITSVertexerCPU_tanLambdaCut_011.root

cd -