#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <unordered_map>
#include <list>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"

void postProcess(const string fileName = "/data1/ruben/pbpbVtx/dbg_ITSVertexerCPU.root",
                 const string fileNameLabels = "/data1/ruben/pbpbVtx/label2Track0.root")
{
    // ROOT::EnableImplicitMT();
    auto vertInfo = ROOT::RDataFrame("verticesInfo", fileName);
    auto histEvtIds = vertInfo.Histo1D({"recEventID", "Reconstructed Events", 151, -1.5, 149.5}, "eventId");
    // histEvtIds->DrawClone();

    // residuals
    std::unordered_map<int, std::array<float, 3>> umap;
    auto labels2Tracks = ROOT::RDataFrame("Labels2Tracks", fileNameLabels);
    auto funzione = [&umap](std::vector<o2::MCCompLabel> labels, std::vector<o2::MCTrack> t) {
        for (int i{0}; i < (int)labels.size(); ++i)
        {
            if (t[i].getMotherTrackId() == -1)
            {
                umap.emplace(labels[i].getEventID(), std::array<float, 3>{(float)t[i].GetStartVertexCoordinatesX(), (float)t[i].GetStartVertexCoordinatesY(), (float)t[i].GetStartVertexCoordinatesZ()});
            }
        }
    };
    auto hResX = new TH1F("resX", "Residual X", 300u, -0.02f, 0.02f);
    auto hResY = new TH1F("resY", "Residual Y", 300u, -0.02f, 0.02f);
    auto hResZ = new TH1F("resZ", "Residual Z", 300u, -0.02f, 0.02f);

    labels2Tracks.Foreach(funzione, {"MCLabels", "Tracks"});
    auto funzione2 = [&](int evtId, float x, float y, float z) {
        auto it = umap.find(evtId);
        if (it != umap.end())
        {
            hResX->Fill(x - it->second[0]);
            hResY->Fill(y - it->second[1]);
            hResZ->Fill(z - it->second[2]);
            umap.erase(it);
        }
    };
    vertInfo.Foreach(funzione2, {"eventId", "xCoord", "yCoord", "zCoord"});
    auto canvasX = new TCanvas("resX", "resX", 800, 600);
    hResX->Draw();
    auto canvasY = new TCanvas("resY", "resY", 800, 600);
    hResY->Draw();
    auto canvasZ = new TCanvas("resZ", "resZ", 800, 600);
    hResZ->Draw();
    canvasX->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResX.png", "r");
    canvasY->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResY.png", "r");
    canvasZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResZ.png", "r");
}

#endif