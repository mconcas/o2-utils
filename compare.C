#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsParameters/GRPObject.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TFile.h"

#include <iostream>
#include <vector>

int compare(const int inspEvt = -1, const int numEvents = 1,
            const std::string inputClustersITS = "o2clus_its.root",
            const std::string paramfilename = "o2sim_par.root")
{
    const std::string path = "./";
    TChain itsClusters("o2sim");
    itsClusters.AddFile((path + inputClustersITS).data());

    // Setup Runtime DB
    TFile paramFile(paramfilename.data());
    paramFile.Get("FairGeoParSet");
    auto gman = o2::ITS::GeometryTGeo::Instance();
    gman->fillMatrixCache(
        o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                            o2::TransformType::L2G)); // request cached transforms

    // get clusters
    std::vector<o2::ITSMFT::Cluster> *clusters = nullptr;
    itsClusters.SetBranchAddress("ITSCluster", &clusters);

    // get labels
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels = nullptr;
    itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);

    const int stopAt =
        (inspEvt == -1) ? itsClusters.GetEntries() : inspEvt + numEvents;
    std::uint32_t roFrame{0};
    o2::ITS::ROframe frame(-123);

    for (int iEvent = (inspEvt == -1) ? 0 : inspEvt; iEvent < stopAt; ++iEvent)
    {
        std::cout << "Event " << iEvent << ": \n";
        itsClusters.GetEntry(iEvent);
        int nclLeft = clusters->size();
        int nclUsed =
            o2::ITS::IOUtils::loadROFrameData(roFrame, frame, clusters, labels);
        if (nclUsed)
        {
            // use the frame
            std::cout << "Hello, world!\n";
        }
    }
    return 0;
}
#endif