#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <memory>
#include <TChain.h>
#include <TGraph.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TBranch.h>
#include <TTreeReader.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <TMultiGraph.h>

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

#include "GPUO2Interface.h"
#include "GPUReconstruction.h"
#include "GPUChainITS.h"

#endif

namespace o2
{
namespace its
{

void arrangeClusters(ROframe *event, std::array<std::vector<Cluster>, constants::its::LayersNumberVertexer> &mClusters)
{

    for (int iLayer{0}; iLayer < constants::its::LayersNumberVertexer; ++iLayer)
    {
        const auto &currentLayer{event->getClustersOnLayer(iLayer)};
        const size_t clustersNum{currentLayer.size()};
        if (clustersNum > 0)
        {
            if (clustersNum > mClusters[iLayer].capacity())
            {
                mClusters[iLayer].reserve(clustersNum);
            }
            for (unsigned int iCluster{0}; iCluster < clustersNum; ++iCluster)
            {
                mClusters[iLayer].emplace_back(iLayer, currentLayer.at(iCluster));
            }
            std::sort(mClusters[iLayer].begin(), mClusters[iLayer].end(), [](Cluster &cluster1, Cluster &cluster2) {
                return cluster1.indexTableBinIndex < cluster2.indexTableBinIndex;
            });
        }
    }
}
} // namespace its
} // namespace o2

void plotClusters(const int startAt,
                  const int stopAt,
                  std::vector<o2::itsmft::ROFRecord> *rofs,
                  std::vector<o2::itsmft::Cluster> *clusters,
                  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels)
{
    std::array<std::vector<o2::its::Cluster>, o2::its::constants::its::LayersNumberVertexer> itsclusters;
    gStyle->SetOptStat(0);
    TH1F *histClus0Phi =
        new TH1F("Layer 0", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus0R =
        new TH1F("Layer 0", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus0Z =
        new TH1F("Layer 0", "Z coordinate, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -15.5f, 15.5f);

    // Setup histograms
    histClus0Phi->SetMinimum(0);
    histClus0Phi->SetFillColor(TColor::GetColor("#29335C"));
    histClus0Phi->SetLineColor(TColor::GetColor("#29335C"));
    histClus0Phi->SetFillStyle(1001);
    histClus0Phi->SetMinimum(0);
    histClus0R->SetFillColor(TColor::GetColor("#29335C"));
    histClus0R->SetLineColor(TColor::GetColor("#29335C"));
    histClus0R->SetFillStyle(1001);
    histClus0Z->SetFillColor(TColor::GetColor("#29335C"));
    histClus0Z->SetLineColor(TColor::GetColor("#29335C"));
    histClus0Z->SetFillStyle(1001);

    TH1F *histClus1Phi =
        new TH1F("Layer 1", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus1R =
        new TH1F("Layer 1", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus1Z =
        new TH1F("Layer 1", "Z coordinate Z, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -16.5f, 16.5f);

    histClus1Phi->SetFillColor(TColor::GetColor("#E4572E"));
    histClus1Phi->SetLineColor(TColor::GetColor("#E4572E"));
    histClus1Phi->SetFillStyle(1001);
    histClus1R->SetFillColor(TColor::GetColor("#E4572E"));
    histClus1R->SetLineColor(TColor::GetColor("#E4572E"));
    histClus1R->SetFillStyle(1001);
    histClus1Z->SetFillColor(TColor::GetColor("#E4572E"));
    histClus1Z->SetLineColor(TColor::GetColor("#E4572E"));
    histClus1Z->SetFillStyle(1001);

    TH1F *histClus2Phi =
        new TH1F("Layer 2", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus2R =
        new TH1F("Layer 2", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus2Z =
        new TH1F("Layer 2", "Z coordinate Z, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -16.5f, 16.5f);
    histClus1Phi->SetMinimum(0);
    histClus2Phi->SetFillColor(TColor::GetColor("#F3A712"));
    histClus2Phi->SetLineColor(TColor::GetColor("#F3A712"));
    histClus2Phi->SetFillStyle(1001);
    histClus2R->SetFillColor(TColor::GetColor("#F3A712"));
    histClus2R->SetLineColor(TColor::GetColor("#F3A712"));
    histClus2R->SetFillStyle(1001);
    histClus2Z->SetMinimum(0);
    histClus2Z->SetFillColor(TColor::GetColor("#F3A712"));
    histClus2Z->SetLineColor(TColor::GetColor("#F3A712"));
    histClus2Z->SetFillStyle(1001);

    TH2D *rphi = new TH2D("R-#phi", "R vs #phi, 150 events PbPb minBias; #phi (rad); R (cm)", 400, 0.f, TMath::TwoPi(), 400, 2.f, 4.2f);

    for (int iROfCount{startAt}; iROfCount < stopAt; ++iROfCount)
    {
        for (auto &layer : itsclusters)
            layer.clear();
        auto &rof = (*rofs)[iROfCount];
        o2::its::ROframe frame{iROfCount}; // to get meaningful roframeId
        // std::cout << "ROframe: " << iROfCount << std::endl;
        int nclUsed = o2::its::ioutils::loadROFrameData(rof, frame, clusters, labels);
        o2::its::arrangeClusters(&frame, itsclusters);

        for (auto &cluster : itsclusters[0])
        {
            histClus0Phi->Fill(cluster.phiCoordinate);
            histClus0R->Fill(cluster.rCoordinate);
            histClus0Z->Fill(cluster.zCoordinate);
            rphi->Fill(cluster.phiCoordinate, cluster.rCoordinate);
        }
        for (auto &cluster : itsclusters[1])
        {
            histClus1Phi->Fill(cluster.phiCoordinate);
            histClus1R->Fill(cluster.rCoordinate);
            histClus1Z->Fill(cluster.zCoordinate);
            rphi->Fill(cluster.phiCoordinate, cluster.rCoordinate);
        }
        for (auto &cluster : itsclusters[2])
        {
            histClus2Phi->Fill(cluster.phiCoordinate);
            histClus2R->Fill(cluster.rCoordinate);
            histClus2Z->Fill(cluster.zCoordinate);
            rphi->Fill(cluster.phiCoordinate, cluster.rCoordinate);
        }
    }

    auto canvasClustersPhi = new TCanvas("ClustersPhi", "Clusters data phi", 1300, 1000);
    histClus0Phi->SetDirectory(0);
    canvasClustersPhi->cd();
    histClus0Phi->GetYaxis()->SetMaxDigits(2);
    histClus0Phi->Draw();
    histClus1Phi->SetDirectory(0);
    histClus1Phi->GetYaxis()->SetMaxDigits(2);
    histClus1Phi->Draw("same");
    histClus2Phi->SetDirectory(0);
    histClus2Phi->GetYaxis()->SetMaxDigits(2);
    histClus2Phi->Draw("same");

    // Legend
    gStyle->SetLegendTextSize(0.);
    auto legendPhi = new TLegend(0.3, 0.12, 0.7, 0.32);
    legendPhi->AddEntry(histClus0Phi, Form("%s: entries %d ", histClus0Phi->GetName(), (int)histClus0Phi->GetEntries()), "f");
    legendPhi->AddEntry(histClus1Phi, Form("%s: entries %d ", histClus1Phi->GetName(), (int)histClus1Phi->GetEntries()), "f");
    legendPhi->AddEntry(histClus2Phi, Form("%s: entries %d ", histClus2Phi->GetName(), (int)histClus2Phi->GetEntries()), "f");
    legendPhi->Draw();

    /////////////////////////////////
    auto canvasClustersR = new TCanvas("ClustersR", "Clusters data R", 1300, 1000);
    canvasClustersR->cd();
    histClus0R->SetDirectory(0);
    histClus0R->GetYaxis()->SetMaxDigits(0);
    histClus0R->Draw();
    histClus1R->SetDirectory(0);
    histClus1R->GetYaxis()->SetMaxDigits(0);
    histClus1R->Draw("same");
    histClus2R->SetDirectory(0);
    histClus2R->GetYaxis()->SetMaxDigits(0);
    histClus2R->Draw("same");

    // Legend
    gStyle->SetLegendTextSize(0.);
    auto legendR = new TLegend(0.3, 0.68, 0.7, 0.88);
    legendR->AddEntry(histClus0R, Form("%s: entries %d ", histClus0R->GetName(), (int)histClus0R->GetEntries()), "f");
    legendR->AddEntry(histClus1R, Form("%s: entries %d ", histClus1R->GetName(), (int)histClus1R->GetEntries()), "f");
    legendR->AddEntry(histClus2R, Form("%s: entries %d ", histClus2R->GetName(), (int)histClus2R->GetEntries()), "f");
    legendR->Draw();

    ///////////////////////////////
    auto canvasClustersZ = new TCanvas("ClustersZ", "Clusters data Z", 1600, 1000);
    canvasClustersZ->cd();
    histClus0Z->SetDirectory(0);
    histClus0Z->GetYaxis()->SetMaxDigits(2);
    histClus0Z->Draw();
    histClus1Z->SetDirectory(0);
    histClus1Z->GetYaxis()->SetMaxDigits(2);
    histClus1Z->Draw("same");
    histClus2Z->SetDirectory(0);
    histClus2Z->GetYaxis()->SetMaxDigits(2);
    histClus2Z->Draw("same");

    // Legend
    gStyle->SetLegendTextSize(0.);
    auto legendZ = new TLegend(0.3, 0.12, 0.7, 0.32);
    legendZ->AddEntry(histClus0Z, Form("%s: entries %d ", histClus0Z->GetName(), (int)histClus0Z->GetEntries()), "f");
    legendZ->AddEntry(histClus1Z, Form("%s: entries %d ", histClus1Z->GetName(), (int)histClus1Z->GetEntries()), "f");
    legendZ->AddEntry(histClus2Z, Form("%s: entries %d ", histClus2Z->GetName(), (int)histClus2Z->GetEntries()), "f");
    legendZ->Draw();

    ////////////////////
    auto canvasRphi = new TCanvas("R vs phi", "R vs phi", 1600, 1000);
    canvasRphi->SetGridy();
    canvasRphi->cd();
    rphi->SetDirectory(0);
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    rphi->Draw("colz");

    canvasClustersPhi->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersPhi.png", "r");
    canvasClustersR->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersR.png", "r");
    canvasClustersZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersZ.png", "r");
    canvasRphi->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersRPhi.png", "r");
    canvasRphi->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersRPhi.png", "r");
}

std::unordered_map<o2::MCCompLabel, o2::MCTrack> getLabelToTrackMap(TFile *file)
{
    std::unordered_map<o2::MCCompLabel, o2::MCTrack> umap;
    TTreeReader readerl2T("Labels2Tracks", file);
    TTreeReaderValue<std::vector<o2::MCCompLabel>> labels(readerl2T, "MCLabels");
    TTreeReaderValue<std::vector<o2::MCTrack>> tracks(readerl2T, "Tracks");
    while (readerl2T.Next())
    {
        for (size_t i{0}; i < (*labels).size(); ++i)
        {
            umap.emplace((*labels)[i], (*tracks)[i]);
        }
    }

    return umap;
}

void plotDBGCPU(TFile *dbgCPUFile, TFile *l2tiFile)
{
    gStyle->SetOptStat(0);
    TTreeReader readerST("selectedTracklets", dbgCPUFile);
    TTreeReaderValue<float> deltaTanLambdaST(readerST, "deltaTanlambda");
    TTreeReaderValue<float> deltaPhiST(readerST, "deltaPhi");
    TTreeReaderValue<float> c0z(readerST, "cluster0z");
    TTreeReaderValue<float> c1z(readerST, "cluster1z");
    TTreeReaderValue<float> c2z(readerST, "cluster2z");
    TTreeReaderValue<float> c0r(readerST, "cluster0r");
    TTreeReaderValue<float> c1r(readerST, "cluster1r");
    TTreeReaderValue<float> c2r(readerST, "cluster2r");
    TTreeReaderValue<float> c0phi(readerST, "cluster0phi");
    TTreeReaderValue<float> c1phi(readerST, "cluster1phi");
    // TTreeReaderValue<float> c2phi(readerST, "cluster2phi");
    TTreeReaderValue<o2::MCCompLabel> labels0(readerST, "lblClus0");
    TTreeReaderValue<o2::MCCompLabel> labels2(readerST, "lblClus2");
    auto umap = getLabelToTrackMap(l2tiFile);

    // Histos
    TH1F *deltaPhi = new TH1F("deltaPhi", "Clusters: #Delta#phi, 150 events PbPb minBias; #Delta#phi (rad); N_{entries}", 200, 0.f, 0.05f);
    TH1F *deltaZ = new TH1F("deltaZ", "Clusters: #DeltaZ, 150 events PbPb minBias; #DeltaZ (cm); N_{entries}", 200, -35.f, 35.f);
    TH1F *pTdist = new TH1F("pTdistribution", "#pi^{#pm} #it{p}_{T} distribution, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); N_{entries}", 400, 0.f, 5.f);
    TH2F *cluDeltaPhiPt = new TH2F("cluDeltaPhiPt", "Clusters: #it{p}_{T} vs #Delta#phi, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); #Delta#phi (rad)", 400, 0.f, 5.f, 400, 0.f, 0.01f);
    TH2F *cluDeltaZPt = new TH2F("cluDeltaZPt", "Clusters: #DeltaZ vs #it{p}_{T}, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); #DeltaZ (cm)", 400, 0.f, 5.f, 400, -10.f, 10.f);
    TH2F *trkDeltaPhiPt = new TH2F("trkDeltaPhiPt", "Tracklets: |#Delta#phi| vs #it{p}_{T}, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); |#Delta#phi| (rad)", 400, 0.f, 5.f, 400, 0.f, 0.01 * TMath::Pi());
    TH2F *trkDeltaLambdaPt = new TH2F("trkDeltaLambdaPt", "Tracklets: |#Delta#lambda| vs #it{p}_{T}, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); |#Delta#lambda| (rad)", 400, 0.f, 5.f, 400, 0.f, 0.01 * TMath::Pi());

    // TH2F *test = new TH2F("test", "#pi_{tracklets}^{#pm} #it{p}_{T} vs #Delta#lambda, 150 events PbPb minBias; #Delta#phi; #Delta#lambda", 4000, -TMath::Pi(), TMath::Pi(), 4000, -TMath::Pi(), TMath::Pi());

    while (readerST.Next())
    {
        auto it = umap.find(*labels0);

        if (it != umap.end() && it->second.getMotherTrackId() == -1)
        {

            // Only primaries

            deltaZ->Fill((*c0z) - (*c1z));
            cluDeltaPhiPt->Fill(it->second.GetPt(), TMath::Abs((*c0phi) - (*c1phi)));
            cluDeltaZPt->Fill(it->second.GetPt(), (*c0z) - (*c1z));
            trkDeltaPhiPt->Fill(it->second.GetPt(), TMath::Abs(*deltaPhiST));
            trkDeltaLambdaPt->Fill(it->second.GetPt(), TMath::Abs(TMath::ATan2((*c1z - *c0z), (*c1r - *c0r)) - TMath::ATan2((*c2z - *c1z), (*c2r - *c1r))));
            // Only pions
            if (TMath::Abs(it->second.GetPdgCode()) == 211)
            {
                pTdist->Fill(it->second.GetPt());
                deltaPhi->Fill(TMath::Abs((*c0phi) - (*c1phi)));
            }
            // Only specific Pt range
            // if (it->second.GetPt() < 0.2 && it->second.GetPt() > 0.1)
            // {
            //     test->Fill(*deltaPhiST, TMath::ATan2((*c1z - *c0z), (*c1r - *c0r)) - TMath::ATan2((*c2z - *c1z), (*c2r - *c1r)), it->second.GetPt());
            // }
            // }
            umap.erase(it);
        }
    }

    // Draw section
    gStyle->SetPalette(kBird);

    auto canvasDeltaPhi = new TCanvas("deltaPhi", "deltaPhi", 800, 600);
    canvasDeltaPhi->SetGrid();
    canvasDeltaPhi->cd();
    canvasDeltaPhi->SetLogy();
    deltaPhi->SetDirectory(0);
    deltaPhi->SetLineColor(kBlack);
    deltaPhi->SetFillColor(kBlue - 8);
    deltaPhi->SetFillStyle(3015);
    deltaPhi->GetYaxis()->SetNdivisions(10);
    deltaPhi->GetXaxis()->SetNdivisions(20);
    deltaPhi->GetYaxis()->SetMaxDigits(2);
    deltaPhi->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendDeltaPhi = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendDeltaPhi->AddEntry(pTdist, Form("Entries: %d ", (int)deltaPhi->GetEntries()), "LF");
    legendDeltaPhi->Draw();

    auto canvasDZ = new TCanvas("deltaZ", "deltaZ", 800, 600);
    canvasDZ->SetLogy();
    canvasDZ->SetGrid();
    canvasDZ->cd();
    deltaZ->SetDirectory(0);
    deltaZ->SetLineColor(kBlack);
    deltaZ->SetFillColor(kBlue - 8);
    deltaZ->SetFillStyle(3015);
    deltaZ->GetYaxis()->SetMaxDigits(2);
    deltaZ->Draw();

    auto canvasPt = new TCanvas("PiPt", "PiPt", 800, 600);
    // canvasPt->SetLogy();
    canvasPt->SetGrid();
    canvasPt->cd();
    canvasPt->SetLogy();
    pTdist->SetDirectory(0);
    pTdist->SetLineColor(kBlack);
    pTdist->SetFillColor(kBlue - 8);
    pTdist->SetFillStyle(3015);
    pTdist->GetYaxis()->SetMaxDigits(2);
    pTdist->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendPtdist = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendPtdist->AddEntry(pTdist, Form("Entries: %d ", (int)pTdist->GetEntries()), "LF");
    legendPtdist->Draw();

    auto canvasPtPhi = new TCanvas("deltaPtPhi", "deltaPtPhi", 800, 600);
    // canvasPtPhi->SetLogy();
    canvasPtPhi->SetGrid();
    canvasPtPhi->cd();
    cluDeltaPhiPt->SetDirectory(0);
    cluDeltaPhiPt->Draw("colz");

    auto canvasClusPtZ = new TCanvas("cluDeltaPtZ", "cluDeltaPtZ", 800, 600);
    canvasClusPtZ->SetGrid();
    canvasClusPtZ->cd();
    cluDeltaZPt->SetDirectory(0);
    cluDeltaZPt->Draw("colz");

    auto canvasTrackPtDeltaLambda = new TCanvas("trkDeltaLambdaPt", "trkDeltaLambdaPt", 800, 600);
    canvasTrackPtDeltaLambda->SetGrid();
    canvasTrackPtDeltaLambda->cd();
    // canvasTrackPtDeltaLambda->SetLogy();
    trkDeltaLambdaPt->SetDirectory(0);
    trkDeltaLambdaPt->Draw("colz");

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendTrkDeltaLambdaPt = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendTrkDeltaLambdaPt->AddEntry(trkDeltaLambdaPt, Form("Entries: %d ", (int)trkDeltaLambdaPt->GetEntries()), "l");
    legendTrkDeltaLambdaPt->Draw();

    auto canvasTrackPtPhi = new TCanvas("trkDeltaZPt", "trkDeltaZPt", 800, 600);
    canvasTrackPtPhi->SetGrid();
    canvasTrackPtPhi->cd();
    trkDeltaPhiPt->SetDirectory(0);
    trkDeltaPhiPt->Draw("colz");

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendTrkDeltaPhiPt = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendTrkDeltaPhiPt->AddEntry(trkDeltaPhiPt, Form("Entries: %d ", (int)trkDeltaPhiPt->GetEntries()), "l");
    legendTrkDeltaPhiPt->Draw();

    // auto canvasTest = new TCanvas("test", "test", 800, 600);
    // canvasTest->SetGrid();
    // canvasTest->cd();
    // test->SetDirectory(0);
    // test->Draw("colz");

    canvasDeltaPhi->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersDeltaPhi.png", "r");
    canvasDZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersDeltaZeta.png", "r");
    canvasPt->SaveAs("/home/mconcas/cernbox/thesis_pictures/clusters0Pt.png", "r");
    canvasPtPhi->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersPtvsDeltaPhi.png", "r");
    canvasClusPtZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/clustersPtvsDeltaZ.png", "r");
    canvasTrackPtPhi->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackletsPtvsDeltaPhi.png", "r");
    canvasTrackPtDeltaLambda->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackletsPtvsDeltaZ.png", "r");
}

void plotPhiCutVariation(TFile *l2tiFile, std::vector<TFile *> fileVector)
{
    float x[10] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
    float good_v[10];
    float fake_v[10];
    float good_norm_v[10];
    int counter{0};
    auto umap = getLabelToTrackMap(l2tiFile);
    int primaries{0};
    for(auto& l: umap) {
        if(l.second.getMotherTrackId() == -1) {
            primaries++;
        }
    }
    const int norm = umap.size();
    for (auto fileptr : fileVector)
    {
        int fake = 0;
        int good = 0;

        TTreeReader readerPhiCV("selectedTracklets", fileptr);
        TTreeReaderValue<unsigned char> validated(readerPhiCV, "isValidated");
        while (readerPhiCV.Next())
        {
            if (*validated)
                ++good;
            else
                ++fake;
        }

        good_v[counter] = good / (float)(good + fake);
        fake_v[counter] = fake / (float)(good + fake);
        good_norm_v[counter] = good / (float)(primaries);
        ++counter;
    }

    auto goodFake = new TCanvas("goodFake", "goodFake", 800, 600);
    goodFake->SetGrid();
    goodFake->cd();
    TGraph *grGood = new TGraph(10, x, good_v);
    grGood->SetLineColor(TColor::GetColor("#F9627D"));
    grGood->SetMarkerStyle(22);
    grGood->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood->SetLineWidth(2);
    TGraph *grGoodNorm = new TGraph(10, x, good_norm_v);
    grGoodNorm->SetFillColor(TColor::GetColor("#5B3758"));
    TGraph *grFake = new TGraph(10, x, fake_v);
    grFake->SetLineColor(TColor::GetColor("#83B692"));
    grFake->SetMarkerStyle(23);
    grFake->SetMarkerColor(TColor::GetColor("#83B692"));
    grFake->SetLineWidth(2);
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Fake/Good tracklets vs #Delta#phi selection, 150 PbPb events minBias; #Delta#phi (rad); efficiency");
    gPad->Modified();
    mg->GetXaxis()->SetNdivisions(20);
    mg->GetXaxis()->SetLimits(0.f, 0.11f);
    mg->SetMinimum(0.f);
    mg->SetMaximum(1.f);
    mg->Add(grGoodNorm, "APB");
    mg->Add(grGood, "APL");
    mg->Add(grFake, "APL");
    mg->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendGoodVsFake = new TLegend(0.5, 0.5, 0.85, 0.68);
    legendGoodVsFake->AddEntry(grGood, "real/found ", "lp");
    legendGoodVsFake->AddEntry(grFake, "fake/found ", "lp");
    legendGoodVsFake->AddEntry(grGoodNorm, "real/MC primaries ", "f");
    legendGoodVsFake->Draw();
    goodFake->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackletingEfficiencies.png", "r");
}

int plots(const int inspEvt = -1,
          const int numEvents = 1,
          const std::string inputClustersITS = "o2clus_its.root",
          const std::string inputGRP = "o2sim_grp.root",
          const std::string simfilename = "o2sim.root",
          const std::string paramfilename = "O2geometry.root",
          const std::string dbgcpufilename = "dbg_ITSVertexerCPU_PureMC.root",
          const std::string labl2trackinfofilename = "label2Track0.root",
          const std::string phiCutVariationDir = "phiCutVariationData",
          const std::string path = "./")
{
    // file load and stuff
    const auto grp = o2::parameters::GRPObject::loadFrom(path + inputGRP);
    const bool isITS = grp->isDetReadOut(o2::detectors::DetID::ITS);
    const bool isContITS = grp->isDetContinuousReadOut(o2::detectors::DetID::ITS);
    std::cout << "ITS is in " << (isContITS ? "CONTINUOS" : "TRIGGERED") << " readout mode" << std::endl;
    TChain itsClusters("o2sim");
    itsClusters.AddFile((path + inputClustersITS).data());

    // Setup Runtime DB
    TFile paramFile((path + paramfilename).data());
    paramFile.Get("FAIRGeom");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                                              o2::TransformType::L2G)); // request cached transforms

    // o2sim
    TChain mcHeaderTree("o2sim");
    mcHeaderTree.AddFile((path + simfilename).data());
    o2::dataformats::MCEventHeader *mcHeader = nullptr;
    mcHeaderTree.SetBranchAddress("MCEventHeader.", &mcHeader);

    std::vector<o2::itsmft::Cluster> *clusters = nullptr;
    itsClusters.SetBranchAddress("ITSCluster", &clusters);

    TChain itsClustersROF("ITSClustersROF");
    itsClustersROF.AddFile((path + inputClustersITS).data());

    std::vector<o2::itsmft::ROFRecord> *rofs = nullptr;
    itsClustersROF.SetBranchAddress("ITSClustersROF", &rofs);
    itsClustersROF.GetEntry(0);

    o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels = nullptr;
    itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);

    // dbg CPU
    TFile dbgCPUFile((path + dbgcpufilename).data());

    // labels to tracks
    TFile l2tiFile((path + labl2trackinfofilename).data());

    std::vector<TFile *> phiDBGFiles;
    for (int i{1}; i < 11; ++i)
    {
        phiDBGFiles.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_PhiCut_0%02d.root", path.data(), phiCutVariationDir.data(), i)));
    }
    
    // config
    const int stopAt = (inspEvt == -1) ? rofs->size() : inspEvt + numEvents;
    const int startAt = (inspEvt == -1) ? 0 : inspEvt;

    itsClusters.GetEntry(0);
    mcHeaderTree.GetEntry(0);

    // bC01->GetEntry(0);

    // Hereafter: direct calls to plotting functions
    // plotClusters(startAt, stopAt, rofs, clusters, labels);
    // plotDBGCPU(&dbgCPUFile, &l2tiFile);
    plotPhiCutVariation(&l2tiFile, phiDBGFiles);

    return 0;
}