#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <memory>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtuple.h>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DetectorsBase/GeometryManager.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"

#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#endif

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;
using namespace o2::gpu;

float signnum_c(float x)
{
    if (x > 0.0)
        return 1.0;
    if (x < 0.0)
        return -1.0;
    return x;
}

int plop(const std::string inputClustersITS = "o2clus_its.root",
         const std::string inputGRP = "o2sim_grp.root",
         // const std::string simfilename = "o2sim_Kine.root",
         const std::string path = "./")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCMYK);
    const auto grp = o2::parameters::GRPObject::loadFrom(path + inputGRP);
    const bool isITS = grp->isDetReadOut(o2::detectors::DetID::ITS);
    const bool isContITS = grp->isDetContinuousReadOut(o2::detectors::DetID::ITS);
    std::cout << "ITS is in " << (isContITS ? "CONTINUOS" : "TRIGGERED") << " readout mode" << std::endl;
    TChain itsClusters("o2sim");
    itsClusters.AddFile((path + inputClustersITS).data());

    o2::base::GeometryManager::loadGeometry(path);
    o2::its::GeometryTGeo *geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                                              o2::TransformType::L2G)); // request cached transforms

    // // Get event header
    // TChain mcHeaderTree("o2sim");
    // mcHeaderTree.AddFile((path + simfilename).data());
    // o2::dataformats::MCEventHeader *mcHeader = nullptr;
    // if (!mcHeaderTree.GetBranch("MCEventHeader."))
    // {
    //     LOG(FATAL) << "Did not find MC event header in the input header file.";
    // }
    // mcHeaderTree.SetBranchAddress("MCEventHeader.", &mcHeader);

    if (!itsClusters.GetBranch("ITSCluster"))
    {
        LOG(FATAL) << "Did not find ITS clusters branch ITSClusters in the input tree";
    }
    std::vector<o2::itsmft::Cluster> *clusters = nullptr;
    itsClusters.SetBranchAddress("ITSCluster", &clusters);

    if (!itsClusters.GetBranch("ITSClustersROF"))
    {
        LOG(FATAL) << "Did not find ITS clusters branch ITSClustersROF in the input tree";
    }
    std::vector<o2::itsmft::ROFRecord> *rofs = nullptr;
    itsClusters.SetBranchAddress("ITSClustersROF", &rofs);
    itsClusters.GetEntry(0);

    TH2F *zr = new TH2F("zr", "Clusters on the longitudinal IB section;z (cm);r (cm)", 1000, -15.f, 15.f, 1000, -5.f, 5.f);
    TH2F *tr = new TH2F("tr", "Clusters on the transverse IB section;x (cm);y (cm)", 1000, -5.f, 5.f, 1000, -5.f, 5.f);
    for (size_t iROfCount{0}; iROfCount < rofs->size(); ++iROfCount)
    {
        auto &rof = (*rofs)[iROfCount];
        o2::its::ROframe frame(iROfCount); // to get meaningful roframeId
        std::cout << "ROframe: " << iROfCount << std::endl;
        // int nclUsed = o2::its::ioutils::loadROFrameData(rof, frame, gsl::span(clusters->data(), clusters->size()), labels);
        frame.clear();
        o2::its::GeometryTGeo *geom = o2::its::GeometryTGeo::Instance();
        geom->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));
        int clusterId{0};

        auto first = rof.getFirstEntry();
        gsl::span<const o2::itsmft::Cluster> span = gsl::span(clusters->data(), clusters->size());
        auto clusters_in_frame = rof.getROFData(span);
        for (auto &c : clusters_in_frame)
        {
            int layer = geom->getLayer(c.getSensorID());

            /// Clusters are stored in the tracking frame
            auto xyz = c.getXYZGloRot(*geom);
            // std::cout << "\tcluster " << clusterId << " " << xyz.x() << " " << xyz.y() << " " << xyz.z() << " " << std::endl;
            zr->Fill(xyz.z(), signnum_c(xyz.y()) * TMath::Sqrt(xyz.x() * xyz.x() + xyz.y() * xyz.y()));
            tr->Fill(xyz.x(), xyz.y());
            ++clusterId;
        }
    }
    TCanvas *czr = new TCanvas("zrc", "zrc", 1600, 1600);
    czr->cd();
    zr->SetDirectory(0);
    zr->Draw("colz");
    czr->SaveAs("/home/mconcas/cernbox/thesis_pictures/ZRdistClusters.png", "r");
    TCanvas *ctr = new TCanvas("ctr", "ctr", 1600, 1600);
    ctr->cd();
    tr->SetDirectory(0);
    tr->Draw("colz");
    ctr->SaveAs("/home/mconcas/cernbox/thesis_pictures/TransverseClusters.png", "r");
    return 0;
}