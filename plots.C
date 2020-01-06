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
#include <THStack.h>

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

// #include "GPUO2Interface.h"
// #include "GPUReconstruction.h"
// #include "GPUChainITS.h"

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

const int kBlueC = TColor::GetColor("#1f78b4");
const int kBlueCT = TColor::GetColorTransparent(kBlueC, 0.9);
const int kRedC = TColor::GetColor("#e31a1c");
const int kRedCT = TColor::GetColorTransparent(kRedC, 0.2);
const int kPurpleC = TColor::GetColor("#911eb4");
const int kPurpleCT = TColor::GetColorTransparent(kPurpleC, 0.2);
const int kOrangeC = TColor::GetColor("#ff7f00");
const int kOrangeCT = TColor::GetColorTransparent(kOrangeC, 0.5);
const int kGreenC = TColor::GetColor("#33a02c");
const int kGreenCT = TColor::GetColorTransparent(kGreenC, 0.2);
const int kMagentaC = TColor::GetColor("#f032e6");
const int kMagentaCT = TColor::GetColorTransparent(kMagentaC, 0.2);
const int kYellowC = TColor::GetColor("#ffe119");
const int kYellowCT = TColor::GetColorTransparent(kYellowC, 0.2);
const int kBrownC = TColor::GetColor("#b15928");
const int kBrownCT = TColor::GetColorTransparent(kBrownC, 0.2);

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

void plotClusters(const int startAt,
                  const int stopAt,
                  std::vector<o2::itsmft::ROFRecord> *rofs,
                  std::vector<o2::itsmft::Cluster> *clusters,
                  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels)
{
    std::array<std::vector<o2::its::Cluster>, o2::its::constants::its::LayersNumberVertexer> itsclusters;
    gStyle->SetOptStat(0);
    TH1F *histClus0Phi =
        new TH1F("Layer0Phi", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus0R =
        new TH1F("Layer0R", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus0Z =
        new TH1F("Layer0Z", "Z coordinate, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -15.5f, 15.5f);

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
        new TH1F("Layer1Phi", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus1R =
        new TH1F("Layer1R", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus1Z =
        new TH1F("Layer1Z", "Z coordinate Z, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -16.5f, 16.5f);

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
        new TH1F("Layer2Phi", "Azimuthal angle #phi, 150 events PbPb minBias;#phi (rad); N_{clusters}", 400, 0.f, TMath::TwoPi());
    TH1F *histClus2R =
        new TH1F("Layer2R", "Radial coordinate R, 150 events PbPb minBias;R (cm); N_{clusters}", 400, 2.f, 4.5f);
    TH1F *histClus2Z =
        new TH1F("Layer2Z", "Z coordinate Z, 150 events PbPb minBias;z (cm); N_{clusters}", 500, -16.5f, 16.5f);
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
    TH1F *deltaPhi = new TH1F("deltaPhi", "#Delta#phi for correlated clusters L_{0}-L_{1}; #Delta#phi (rad); N_{entries}", 200, 0.f, 0.05f);
    TH1F *deltaZ = new TH1F("deltaZ", "#DeltaZ for correlated clusters L_{0}-L_{1}; #DeltaZ (cm); N_{entries}", 200, -35.f, 35.f);
    TH1F *pTdist = new TH1F("pTdistribution", "#pi^{#pm}: #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); N_{entries}", 400, 0.f, 5.f);
    TH2F *cluDeltaPhiPt = new TH2F("cluDeltaPhiPt", "Clusters: #it{p}_{T} vs #Delta#phi, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); #Delta#phi (rad)", 400, 0.f, 5.f, 400, 0.f, 0.01f);
    TH2F *cluDeltaZPt = new TH2F("cluDeltaZPt", "Clusters: #DeltaZ vs #it{p}_{T}, 150 events PbPb minBias; #it{p}_{T} (GeV/#it{c}); #DeltaZ (cm)", 400, 0.f, 5.f, 400, -10.f, 10.f);
    TH2F *trkDeltaPhiPt = new TH2F("trkDeltaPhiPt", "|#Delta#phi| vs #it{p}_{T} for correlated primary tracklets; #it{p}_{T} (GeV/#it{c}); |#Delta#phi| (rad)", 1000, 0.f, 5.f, 600, 0.f, 0.1 * TMath::Pi());
    TH2F *trkDeltaTanLambdaPt = new TH2F("trkDeltaTanLambdaPt", "|#Deltatan#lambda| vs #it{p}_{T} for correlated primary tracklets; #it{p}_{T} (GeV/#it{c}); |#Deltatan#lambda|", 1000, 0.f, 5.f, 1000, 0.f, 0.1 * TMath::Pi());

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
            trkDeltaTanLambdaPt->Fill(it->second.GetPt(), TMath::Abs(TMath::Tan(TMath::ATan2((*c1z - *c0z), (*c1r - *c0r))) - TMath::Tan(TMath::ATan2((*c2z - *c1z), (*c2r - *c1r)))));
            deltaPhi->Fill(TMath::Abs((*c0phi) - (*c1phi)));
            // Only pions
            if (TMath::Abs(it->second.GetPdgCode()) == 211)
            {
                pTdist->Fill(it->second.GetPt());
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

    auto canvasDeltaPhi = new TCanvas("deltaPhi", "deltaPhi", 800, 800);
    canvasDeltaPhi->SetGrid();
    canvasDeltaPhi->cd();
    canvasDeltaPhi->SetLogy();
    deltaPhi->SetDirectory(0);
    deltaPhi->SetLineColor(kBlack);
    deltaPhi->SetFillColor(kOrangeCT);
    // deltaPhi->SetFillStyle(3015);
    deltaPhi->GetYaxis()->SetNdivisions(10);
    deltaPhi->GetXaxis()->SetNdivisions(10);
    deltaPhi->GetYaxis()->SetMaxDigits(2);
    deltaPhi->GetYaxis()->SetMaxDigits(3);
    deltaPhi->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendDeltaPhi = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendDeltaPhi->SetHeader("150 events PbPb MB");
    legendDeltaPhi->AddEntry(deltaPhi, Form("Entries: %d ", (int)deltaPhi->GetEntries()), "LF");
    legendDeltaPhi->AddEntry(deltaPhi, Form("#LT#Delta#phi#GT=%2.4f rad", deltaPhi->GetMean()), "");
    legendDeltaPhi->Draw();

    auto canvasDZ = new TCanvas("deltaZ", "deltaZ", 800, 800);
    canvasDZ->SetLogy();
    canvasDZ->SetGrid();
    canvasDZ->cd();
    deltaZ->SetDirectory(0);
    deltaZ->SetLineColor(kBlack);
    deltaZ->SetFillColor(kOrangeCT);
    // deltaZ->SetFillStyle(3015);
    deltaZ->GetYaxis()->SetMaxDigits(2);
    deltaZ->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legenddeltaZ = new TLegend(0.6, 0.78, 0.87, 0.88);
    legenddeltaZ->SetHeader("150 events PbPb MB");
    legenddeltaZ->AddEntry(deltaZ, Form("Entries: %d ", (int)deltaZ->GetEntries()), "LF");
    legenddeltaZ->AddEntry(deltaZ, Form("#LT#it{p}_{T}#GT= %2.2f GeV/#it{c}", deltaZ->GetMean()), "");
    legenddeltaZ->Draw();

    auto canvasPt = new TCanvas("PiPt", "PiPt", 800, 600);
    // canvasPt->SetLogy();
    canvasPt->SetGrid();
    canvasPt->cd();
    canvasPt->SetLogy();
    pTdist->SetDirectory(0);
    pTdist->SetLineColor(kBlack);
    pTdist->SetFillColor(kBlueCT);
    // pTdist->SetFillStyle(3015);
    pTdist->GetYaxis()->SetMaxDigits(2);
    pTdist->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendPtdist = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendPtdist->SetHeader("150 events PbPb MB");
    legendPtdist->AddEntry(pTdist, Form("Entries: %d ", (int)pTdist->GetEntries()), "LF");
    legendPtdist->AddEntry(pTdist, Form("#LT#it{p}_{T}#GT= %2.2f GeV/#it{c}", pTdist->GetMean()), "");
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

    auto canvasTrackPtDeltaLambda = new TCanvas("trkDeltaTanLambdaPt", "trkDeltaTanLambdaPt", 1200, 1200);
    canvasTrackPtDeltaLambda->SetGrid();
    canvasTrackPtDeltaLambda->cd();
    canvasTrackPtDeltaLambda->SetLogy();
    trkDeltaTanLambdaPt->SetDirectory(0);
    trkDeltaTanLambdaPt->Draw("colz");

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendTrkDeltaTanLambdaPt = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendTrkDeltaTanLambdaPt->SetHeader("150 events PbPb MB");
    legendTrkDeltaTanLambdaPt->AddEntry(trkDeltaTanLambdaPt, Form("Entries: %d ", (int)trkDeltaTanLambdaPt->GetEntries()), "l");
    legendTrkDeltaTanLambdaPt->Draw();

    auto canvasTrackPtPhi = new TCanvas("trkDeltaZPt", "trkDeltaZPt", 1200, 1200);
    canvasTrackPtPhi->SetGrid();
    canvasTrackPtPhi->cd();
    canvasTrackPtPhi->SetLogy();
    trkDeltaPhiPt->SetDirectory(0);
    trkDeltaPhiPt->Draw("colz");

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendTrkDeltaPhiPt = new TLegend(0.6, 0.78, 0.87, 0.88);
    legendTrkDeltaPhiPt->SetHeader("150 events PbPb MB");
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
    canvasTrackPtDeltaLambda->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackletsPtvsDeltaTanLambda.png", "r");
}

void plotPhiCutVariation(TFile *l2tiFile, TFile *dbgCPUFile, TFile *dbgCPUFileSingle, std::vector<TFile *> fileVectorTrackleting,
                         std::vector<TFile *> fileVectorSelection)
{
    gStyle->SetBarWidth(0.5);
    std::array<float, 11> x{0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};

    std::array<float, 11> good_v, fake_v, good01_v, fake01_v, good12_v, fake12_v, good_norm_v, good01_norm_v, good12_norm_v;
    int counter{0}, counter01{0}, counter12{0};
    auto umap = getLabelToTrackMap(l2tiFile);
    int primaries{0};
    TTreeReader readerPMC01("combinatorics01", dbgCPUFileSingle);
    TTreeReader readerPMC12("combinatorics12", dbgCPUFileSingle);

    TH1F *histReference = new TH1F("", "", 50, 0.f, 5.f);

    auto norm01{readerPMC01.GetEntries()};
    auto norm12{readerPMC12.GetEntries()};
    for (auto &l : umap)
    {
        if (l.second.getMotherTrackId() == -1 /* && TMath::Abs(l.second.GetPdgCode()) == 211*/)
        {
            primaries++;
            histReference->Fill(l.second.GetPt());
        }
    }
    std::vector<TH1F *> histos(11);
    std::vector<TH1F *> histosDivided(11);
    THStack *hs = new THStack("deltaphiStack", "");

    for (auto fileptr : fileVectorSelection)
    {
        auto umaptmp = getLabelToTrackMap(l2tiFile);
        histos[counter] = new TH1F(Form("histodeltaphi%d", counter), Form("histo%d", counter), 50, 0.f, 5.f);
        int fake{0}, good{0};
        TTreeReader readerPhiCV("selectedTracklets", fileptr);
        TTreeReaderValue<o2::MCCompLabel> labels(readerPhiCV, "lblClus0");
        TTreeReaderValue<unsigned char> validated(readerPhiCV, "isValidated");
        while (readerPhiCV.Next())
        {
            if (*validated)
            {
                auto it = umaptmp.find(*labels);
                if (it != umaptmp.end())
                {

                    histos[counter]->Fill(it->second.GetPt());
                    umaptmp.erase(it);
                }
                ++good;
            }
            else
                ++fake;
        }
        histosDivided[counter] = (TH1F *)histos[counter]->Clone(Form("histoPhi%d", counter));
        histosDivided[counter]->Divide(histos[counter], histReference);
        hs->Add(histosDivided[counter]);
        good_v[counter] = good / (float)(good + fake);
        fake_v[counter] = fake / (float)(good + fake);
        good_norm_v[counter] = good / (float)(primaries);
        ++counter;
    }

    auto deltaPhiPtEff = new TCanvas("deltaphiPtEfficiencies", "deltaphiPtEfficiencies", 800, 600);
    hs->Draw("nostack");

    auto selEfficienciesCanvas = new TCanvas("selEfficiencies", "selEfficiencies", 800, 600);
    selEfficienciesCanvas->SetGrid();
    selEfficienciesCanvas->cd();
    TGraph *grGood = new TGraph(11, x.data(), good_v.data());
    grGood->SetLineColor(TColor::GetColor("#F9627D"));
    grGood->SetMarkerStyle(22);
    grGood->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood->SetLineWidth(2);
    TGraph *grGoodNorm = new TGraph(11, x.data(), good_norm_v.data());
    grGoodNorm->SetFillColor(kPurpleCT);
    grGoodNorm->SetFillStyle(3001);
    TGraph *grFake = new TGraph(11, x.data(), fake_v.data());
    grFake->SetLineColor(TColor::GetColor("#83B692"));
    grFake->SetMarkerStyle(23);
    grFake->SetMarkerColor(TColor::GetColor("#83B692"));
    grFake->SetLineWidth(2);
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Tracklet selection: #Deltatan#lambda=0.01; #Delta#phi (rad); efficiency");
    gPad->Modified();
    mg->GetXaxis()->SetNdivisions(21);
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
    legendGoodVsFake->SetHeader("150 events PbPb MB", "C");
    legendGoodVsFake->AddEntry(grGood, "validated/found ", "lp");
    legendGoodVsFake->AddEntry(grFake, "fake/found ", "lp");
    legendGoodVsFake->AddEntry(grGoodNorm, "validated/MC_{primaries} ", "f");
    legendGoodVsFake->Draw();
    selEfficienciesCanvas->SaveAs("/home/mconcas/cernbox/thesis_pictures/selectionEfficiencies.png", "r");

    TMultiGraph *mg2 = new TMultiGraph();
    std::vector<TGraph *> graphs;
    graphs.resize(11);
    mg2->SetTitle("Tracklet finding:  #Deltatan#lambda=0.01; #it{p}_{T} (GeV/#it{c}); efficiency");
    mg2->SetMinimum(0.f);
    for (auto iHisto{0}; iHisto < 11; ++iHisto)
    {
        graphs[iHisto] = new TGraph(histosDivided[iHisto]);
        graphs[iHisto]->SetMarkerStyle(23);
    }

    mg2->SetMaximum(1.);
    graphs[0]->SetMarkerColor(kRedCT);
    graphs[0]->SetLineColor(kRedCT);
    graphs[5]->SetMarkerColor(kBlueCT);
    graphs[5]->SetLineColor(kBlueCT);
    graphs[10]->SetMarkerColor(kGreenCT);
    graphs[10]->SetLineColor(kGreenCT);

    mg2->Add(graphs[0], "APL");
    mg2->Add(graphs[5], "APL");
    mg2->Add(graphs[10], "APL");

    auto deltaPhiPtEffGraph = new TCanvas("deltaphiPtEfficienciesAsGraphs", "deltaphiPtEfficienciesAsGraphs", 800, 800);
    deltaPhiPtEffGraph->SetGrid();
    mg2->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendEffPhiCut = new TLegend(0.6, 0.14, 0.85, 0.38);
    legendEffPhiCut->SetHeader("150 events PbPb MB", "C");
    legendEffPhiCut->AddEntry(graphs[0], "#Delta#phi=0.005 (rad)", "lp");
    legendEffPhiCut->AddEntry(graphs[5], "#Delta#phi=0.05 (rad)", "lp");
    legendEffPhiCut->AddEntry(graphs[10], "#Delta#phi=0.1 (rad)", "lp");
    legendEffPhiCut->Draw();
    deltaPhiPtEffGraph->SaveAs("/home/mconcas/cernbox/thesis_pictures/phiCutEfficiencyTanLambdaFixed.png", "r");

    //
    for (auto fileptr : fileVectorTrackleting)
    {
        long int fake01{0}, fake12{0}, good01{0}, good12{0};
        TTreeReader readerComb01("combinatorics01", fileptr);
        TTreeReader readerComb12("combinatorics12", fileptr);
        TTreeReaderValue<unsigned char> c01validated(readerComb01, "isValidated");
        TTreeReaderValue<unsigned char> c12validated(readerComb12, "isValidated");
        TTreeReaderValue<o2::MCCompLabel> labels0(readerComb01, "lblClus0");
        TTreeReaderValue<o2::MCCompLabel> labels2(readerComb12, "lblClus2");

        while (readerComb01.Next())
        {
            if (*c01validated)
            {
                ++good01;
            }
            else
                ++fake01;
        }
        while (readerComb12.Next())
        {
            if (*c12validated)
            {
                ++good12;
            }
            else
                ++fake12;
        }

        good01_v[counter01] = good01 / (double)(good01 + fake01);
        fake01_v[counter01] = fake01 / (double)(good01 + fake01);
        good01_norm_v[counter01] = good01 / (double)(norm01);

        good12_v[counter12] = good12 / (double)(good12 + fake12);
        fake12_v[counter12] = fake12 / (double)(good12 + fake12);
        good12_norm_v[counter12] = good12 / (double)(norm12);
        ++counter01;
        ++counter12;
    }
    // -----------------------------------
    auto goodFake01 = new TCanvas("goodFake01", "goodFake01", 800, 800);
    goodFake01->SetGrid();
    goodFake01->cd();
    TGraph *grGood01 = new TGraph(11, x.data(), good01_v.data());
    grGood01->SetLineColor(TColor::GetColor("#F9627D"));
    grGood01->SetMarkerStyle(22);
    grGood01->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood01->SetLineWidth(2);
    TGraph *grGoodNorm01 = new TGraph(11, x.data(), good01_norm_v.data());
    grGoodNorm01->SetFillColor(kBlueCT);
    grGoodNorm01->SetFillStyle(3001);
    TGraph *grFake01 = new TGraph(11, x.data(), fake01_v.data());
    grFake01->SetLineColor(TColor::GetColor("#83B692"));
    grFake01->SetMarkerStyle(23);
    grFake01->SetMarkerColor(TColor::GetColor("#83B692"));
    grFake01->SetLineWidth(2);
    TMultiGraph *mg01 = new TMultiGraph();
    mg01->SetTitle("Tracklet finder L_{0}-L_{1}; #Delta#phi (rad); efficiency");
    gPad->Modified();
    mg01->GetXaxis()->SetNdivisions(21);
    mg01->GetXaxis()->SetLimits(0.f, 0.11f);
    mg01->SetMinimum(0.f);
    // mg01->SetMaximum(1.f);
    mg01->Add(grGoodNorm01, "APB");
    mg01->Add(grGood01, "APL");
    mg01->Add(grFake01, "APL");
    mg01->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendGoodVsFake01 = new TLegend(0.5, 0.5, 0.85, 0.68);
    legendGoodVsFake01->SetHeader("150 events PbPb MB", "C");
    legendGoodVsFake01->AddEntry(grGood01, "validated/found ", "lp");
    legendGoodVsFake01->AddEntry(grFake01, "fake/found ", "lp");
    legendGoodVsFake01->AddEntry(grGoodNorm01, "validated/MC_{matches}", "f");
    legendGoodVsFake01->Draw();
    goodFake01->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackleting01Efficiencies.png", "r");

    // -----------------------------------
    auto goodFake12 = new TCanvas("goodFake12", "goodFake12", 800, 800);
    goodFake12->SetGrid();
    goodFake12->cd();
    TGraph *grGood12 = new TGraph(11, x.data(), good01_v.data());
    grGood12->SetLineColor(TColor::GetColor("#F9627D"));
    grGood12->SetMarkerStyle(22);
    grGood12->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood12->SetLineWidth(2);
    TGraph *grGoodNorm12 = new TGraph(11, x.data(), good01_norm_v.data());
    grGoodNorm12->SetFillColor(kBrownCT);
    grGoodNorm12->SetFillStyle(3001);
    TGraph *grFake12 = new TGraph(11, x.data(), fake01_v.data());
    grFake12->SetLineColor(TColor::GetColor("#83B692"));
    grFake12->SetMarkerStyle(23);
    grFake12->SetMarkerColor(TColor::GetColor("#83B692"));
    grFake12->SetLineWidth(2);
    TMultiGraph *mg12 = new TMultiGraph();
    mg12->SetTitle("Tracklet finder L_{1}-L_{2}; #Delta#phi (rad); efficiency");
    gPad->Modified();
    mg12->GetXaxis()->SetNdivisions(21);
    mg12->GetXaxis()->SetLimits(0.f, 0.11f);
    mg12->SetMinimum(0.f);
    // mg12->SetMaximum(1.f);
    mg12->Add(grGoodNorm12, "APB");
    mg12->Add(grGood12, "APL");
    mg12->Add(grFake12, "APL");
    mg12->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendGoodVsFake12 = new TLegend(0.5, 0.5, 0.85, 0.68);
    legendGoodVsFake12->SetHeader("150 events PbPb MB", "C");
    legendGoodVsFake12->AddEntry(grGood12, "validated/found ", "lp");
    legendGoodVsFake12->AddEntry(grFake12, "fake/found ", "lp");
    legendGoodVsFake12->AddEntry(grGoodNorm12, "validated/MC_{matches} ", "f");
    legendGoodVsFake12->Draw();
    goodFake12->SaveAs("/home/mconcas/cernbox/thesis_pictures/trackleting12Efficiencies.png", "r");
}

void plotTanLambdaVariation(TFile *l2tiFile, std::vector<TFile *> fileVectorSelection, std::vector<TFile *> fileVectorSelectionNarrow)
{
    gStyle->SetBarWidth(0.5);
    float tanLambdaCuts[11] = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011};
    float good_v[11], fake_v[11], good_norm_v[11];
    int counter{0};
    auto umap = getLabelToTrackMap(l2tiFile);
    int primaries{0};
    TH1F *histReference = new TH1F("", "", 50, 0.f, 5.f);
    for (auto &l : umap)
    {
        if (l.second.getMotherTrackId() == -1)
        {
            primaries++;
            histReference->Fill(l.second.GetPt());
        }
    }
    std::vector<TH1F *> histos(11);
    std::vector<TH1F *> histosDivided(11);
    THStack *hs = new THStack("tanLambdaStack", "");
    for (auto fileptr : fileVectorSelection)
    {
        auto umaptmp = getLabelToTrackMap(l2tiFile);
        histos[counter] = new TH1F(Form("histo%d", counter), Form("histo%d", counter), 50, 0.f, 5.f);
        int fake{0}, good{0};
        TTreeReader readerPhiCV("selectedTracklets", fileptr);
        TTreeReaderValue<unsigned char> validated(readerPhiCV, "isValidated");
        TTreeReaderValue<o2::MCCompLabel> labels(readerPhiCV, "lblClus0");
        while (readerPhiCV.Next())
        {
            if (*validated)
            {
                ++good;
                auto it = umaptmp.find(*labels);
                if (it != umaptmp.end())
                {
                    histos[counter]->Fill(it->second.GetPt());
                    umaptmp.erase(it);
                }
            }
            else
                ++fake;
        }
        histosDivided[counter] = (TH1F *)histos[counter]->Clone(Form("histo%d", counter));
        histosDivided[counter]->Divide(histos[counter], histReference);
        hs->Add(histosDivided[counter]);
        good_v[counter] = good / (float)(good + fake);
        fake_v[counter] = fake / (float)(good + fake);
        good_norm_v[counter] = good / (float)(primaries);
        ++counter;
    }

    auto tanLambdaPtEff = new TCanvas("tanLambdaPtEfficiencies", "tanLambdaPtEfficiencies", 800, 600);
    hs->Draw("nostack");

    auto goodFake = new TCanvas("goodFake", "goodFake", 800, 600);
    goodFake->SetGrid();
    goodFake->cd();
    TGraph *grGood = new TGraph(11, tanLambdaCuts, good_v);
    grGood->SetLineColor(TColor::GetColor("#F9627D"));
    grGood->SetMarkerStyle(22);
    grGood->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood->SetLineWidth(2);
    TGraph *grGoodNorm = new TGraph(11, tanLambdaCuts, good_norm_v);
    grGoodNorm->SetFillColor(kPurpleCT);
    grGoodNorm->SetFillStyle(3001);
    TGraph *grFake = new TGraph(11, tanLambdaCuts, fake_v);
    grFake->SetLineColor(TColor::GetColor("#83B692"));
    grFake->SetMarkerStyle(23);
    grFake->SetMarkerColor(TColor::GetColor("#83B692"));
    grFake->SetLineWidth(2);
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Tracklet selection: #Delta#phi=0.1; #Deltatan#lambda; efficiency");
    gPad->Modified();
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
    legendGoodVsFake->SetHeader("150 events PbPb MB", "C");
    legendGoodVsFake->AddEntry(grGood, "validated/found ", "lp");
    legendGoodVsFake->AddEntry(grFake, "fake/found ", "lp");
    legendGoodVsFake->AddEntry(grGoodNorm, "validated/MC_{primaries} ", "f");
    legendGoodVsFake->Draw();
    goodFake->SaveAs("/home/mconcas/cernbox/thesis_pictures/selectionEfficienciesTanLambda.png", "r");

    TMultiGraph *mg2 = new TMultiGraph();
    std::vector<TGraph *> graphs;
    graphs.resize(11);
    mg2->SetTitle("Tracklet finding:  #Delta#phi=0.1; #it{p}_{T} (GeV/#it{c}); efficiency");
    mg2->SetMinimum(0.f);
    for (auto iHisto{0}; iHisto < 11; ++iHisto)
    {
        graphs[iHisto] = new TGraph(histosDivided[iHisto]);
        graphs[iHisto]->SetMarkerStyle(23);
    }

    mg2->SetMaximum(1.);
    graphs[0]->SetMarkerColor(kRedCT);
    graphs[0]->SetLineColor(kRedCT);
    graphs[1]->SetMarkerColor(kBlueCT);
    graphs[1]->SetLineColor(kBlueCT);
    graphs[2]->SetMarkerColor(kPurpleCT);
    graphs[2]->SetLineColor(kPurpleCT);
    graphs[3]->SetMarkerColor(kOrangeCT);
    graphs[3]->SetLineColor(kOrangeCT);
    graphs[4]->SetMarkerColor(kMagentaCT);
    graphs[4]->SetLineColor(kMagentaCT);
    graphs[5]->SetMarkerColor(kBrownCT);
    graphs[5]->SetLineColor(kBrownCT);
    graphs[9]->SetMarkerColor(kGreenCT);
    graphs[9]->SetLineColor(kGreenCT);

    mg2->Add(graphs[0], "APL");
    mg2->Add(graphs[1], "APL");
    mg2->Add(graphs[2], "APL");
    mg2->Add(graphs[3], "APL");
    mg2->Add(graphs[4], "APL");
    mg2->Add(graphs[5], "APL");
    mg2->Add(graphs[9], "APL");

    auto deltaTanLambdaPtCanvas = new TCanvas("deltatanlambdaPtEfficienciesAsGraphs", "deltatanlambdaPtEfficienciesAsGraphs", 800, 800);
    deltaTanLambdaPtCanvas->SetGrid();
    mg2->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendEffTanlambdaCutGraph = new TLegend(0.6, 0.14, 0.85, 0.38);
    legendEffTanlambdaCutGraph->SetHeader("150 events PbPb MB", "C");
    legendEffTanlambdaCutGraph->AddEntry(graphs[0], "#Deltatan#lambda=0.001", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[1], "#Deltatan#lambda=0.002", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[2], "#Deltatan#lambda=0.003", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[3], "#Deltatan#lambda=0.004", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[4], "#Deltatan#lambda=0.005", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[4], "#Deltatan#lambda=0.006", "lp");
    legendEffTanlambdaCutGraph->AddEntry(graphs[9], "#Deltatan#lambda=0.010", "lp");
    legendEffTanlambdaCutGraph->Draw();
    deltaTanLambdaPtCanvas->SaveAs("/home/mconcas/cernbox/thesis_pictures/tanlLambdaEfficiencyPhiCutFixed.png", "r");

    // ------------------------------------------ narrow phi cut
    float good_narrow_v[11], fake_narrow_v[11], good_narrow_norm_v[11];
    int counter_narrow{0};
    for (auto fileptr : fileVectorSelectionNarrow)
    {
        auto umaptmp = getLabelToTrackMap(l2tiFile);
        int fake{0}, good{0};
        TTreeReader readerPhiCV("selectedTracklets", fileptr);
        TTreeReaderValue<unsigned char> validated(readerPhiCV, "isValidated");
        TTreeReaderValue<o2::MCCompLabel> labels(readerPhiCV, "lblClus0");
        while (readerPhiCV.Next())
        {
            if (*validated)
            {
                ++good;
                auto it = umaptmp.find(*labels);
                if (it != umaptmp.end())
                {
                    histos[counter_narrow]->Fill(it->second.GetPt());
                    umaptmp.erase(it);
                }
            }
            else
                ++fake;
        }
        good_narrow_v[counter_narrow] = good / (float)(good + fake);
        fake_narrow_v[counter_narrow] = fake / (float)(good + fake);
        good_narrow_norm_v[counter_narrow] = good / (float)(primaries);
        ++counter_narrow;
    }
    auto goodFake_narrow = new TCanvas("goodFake_narrow", "goodFake_narrow", 800, 600);
    goodFake_narrow->SetGrid();
    goodFake_narrow->cd();
    TGraph *grGood_narrow = new TGraph(11, tanLambdaCuts, good_narrow_v);
    grGood_narrow->SetLineColor(TColor::GetColor("#F9627D"));
    grGood_narrow->SetMarkerStyle(22);
    grGood_narrow->SetMarkerColor(TColor::GetColor("#F9627D"));
    grGood_narrow->SetLineWidth(2);
    TGraph *grGood_narrow_Norm = new TGraph(11, tanLambdaCuts, good_narrow_norm_v);
    grGood_narrow_Norm->SetFillColor(kPurpleCT);
    grGood_narrow_Norm->SetFillStyle(3001);
    TGraph *gr_narrow_Fake = new TGraph(11, tanLambdaCuts, fake_narrow_v);
    gr_narrow_Fake->SetLineColor(TColor::GetColor("#83B692"));
    gr_narrow_Fake->SetMarkerStyle(23);
    gr_narrow_Fake->SetMarkerColor(TColor::GetColor("#83B692"));
    gr_narrow_Fake->SetLineWidth(2);
    TMultiGraph *mg_narrow = new TMultiGraph();
    mg_narrow->SetTitle("Tracklet selection: #Delta#phi=0.005; #Deltatan#lambda; efficiency");
    gPad->Modified();
    mg_narrow->SetMinimum(0.f);
    mg_narrow->SetMaximum(1.f);
    mg_narrow->Add(grGood_narrow_Norm, "APB");
    mg_narrow->Add(grGood_narrow, "APL");
    mg_narrow->Add(gr_narrow_Fake, "APL");
    mg_narrow->Draw("a");

    // Legend
    gStyle->SetLegendTextSize(0.);
    gStyle->SetLegendBorderSize(1);
    auto legendGoodVsFake_narrow = new TLegend(0.5, 0.3, 0.85, 0.48);
    legendGoodVsFake_narrow->SetHeader("150 events PbPb MB", "C");
    legendGoodVsFake_narrow->AddEntry(grGood_narrow, "validated/found ", "lp");
    legendGoodVsFake_narrow->AddEntry(gr_narrow_Fake, "fake/found ", "lp");
    legendGoodVsFake_narrow->AddEntry(grGood_narrow_Norm, "validated/MC_{primaries} ", "f");
    legendGoodVsFake_narrow->Draw();
    goodFake_narrow->SaveAs("/home/mconcas/cernbox/thesis_pictures/selectionEfficienciesTanLambdaNarrow.png", "r");
}

void plotZCorrelations(TFile *fileptr)
{
    gStyle->SetOptStat(0);
    TTreeReader readerPMC01("combinatorics01", fileptr);
    TTreeReader readerPMC12("combinatorics12", fileptr);
    TTreeReaderValue<float> z0(readerPMC01, "c0z");
    TTreeReaderValue<float> z1(readerPMC01, "c1z");
    TTreeReaderValue<float> z11(readerPMC12, "c1z");
    TTreeReaderValue<float> z21(readerPMC12, "c2z");

    TH2F *hist01 = new TH2F("z0z1Correlations", "z_{0}-z_{1} correlations, #Delta#phi=0.005 rad; z_{0} (cm); z_{1} (cm)", 1000, -14.5f, 14.5f, 400, -14.5f, 14.5f);
    while (readerPMC01.Next())
    {
        hist01->Fill(*z0, *z1);
    }
    TH2F *hist12 = new TH2F("z1z2Correlations", "z_{1}-z_{2} correlations, #Delta#phi=0.005 rad; z_{1} (cm); z_{2} (cm)", 1000, -14.5f, 14.5f, 400, -14.5f, 14.5f);
    while (readerPMC12.Next())
    {
        hist12->Fill(*z11, *z21);
    }
    auto canvasz0z1 = new TCanvas("z0z1Correlations", "z0z1Correlations", 800, 800);
    canvasz0z1->SetGrid();
    canvasz0z1->cd();
    hist01->SetDirectory(0);
    hist01->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendCorrZ01 = new TLegend(0.6, 0.28, 0.87, 0.38);
    legendCorrZ01->SetHeader("PbPb single event", "C");
    legendCorrZ01->AddEntry(hist01, Form("Tracklets: %d ", (int)hist01->GetEntries()), "l");
    legendCorrZ01->Draw();

    auto canvasz1z2 = new TCanvas("z1z2Correlations", "z1z2Correlations", 800, 800);
    canvasz1z2->SetGrid();
    canvasz1z2->cd();
    hist12->SetDirectory(0);
    hist12->Draw();

    // Legend
    gStyle->SetLegendBorderSize(1);
    auto legendCorrZ12 = new TLegend(0.6, 0.28, 0.87, 0.38);
    legendCorrZ12->SetHeader("PbPb single event", "C");
    legendCorrZ12->AddEntry(hist12, Form("Tracklets: %d ", (int)hist12->GetEntries()), "l");
    legendCorrZ12->Draw();

    canvasz0z1->SaveAs("/home/mconcas/cernbox/thesis_pictures/z0z1correlations.png", "r");
    canvasz1z2->SaveAs("/home/mconcas/cernbox/thesis_pictures/z1z2correlations.png", "r");
}

int plots(const int inspEvt = -1,
          const int numEvents = 1,
          const std::string inputClustersITS = "o2clus_its.root",
          const std::string inputGRP = "o2sim_grp.root",
          const std::string simfilename = "o2sim.root",
          const std::string paramfilename = "O2geometry.root",
          const std::string dbgcpufilename = "dbg_ITSVertexerCPU_150_PureMC.root",
          const std::string dbgcpufilenameSingle505 = "dbg_ITSVertexerCPU_PureMC_single_505.root",
          const std::string labl2trackinfofilename = "label2Track0.root",
          const std::string phiCutVariationSelectionDir = "phiCutVariationLargeTanLambda/150evts",
          const std::string phiCutVariationTrackletingDir = "phiCutVariationData/single505",
          const std::string tanLambdaCutVariationSelectionDir = "tanLambdaCutVariationData/150evts",
          const std::string tanlambdawithnarrowphipath = "tanlambdaVariation/150evts",
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

    // dbg CPU single event
    TFile dbgcpufileSingle505((path + dbgcpufilenameSingle505).data());

    // labels to tracks
    TFile l2tiFile((path + labl2trackinfofilename).data());

    std::vector<TFile *> phiDBGFilesSelection;
    phiDBGFilesSelection.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_PhiCut_0005.root", path.data(), phiCutVariationSelectionDir.data())));
    for (int i{1}; i < 11; ++i)
    {
        phiDBGFilesSelection.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_PhiCut_0%02d.root", path.data(), phiCutVariationSelectionDir.data(), i)));
    }

    std::vector<TFile *> phiDBGFilesTrackleting;
    phiDBGFilesTrackleting.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_PhiCut_0005.root", path.data(), phiCutVariationTrackletingDir.data())));
    for (int i{1}; i < 11; ++i)
    {
        phiDBGFilesTrackleting.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_PhiCut_0%02d.root", path.data(), phiCutVariationTrackletingDir.data(), i)));
    }

    std::vector<TFile *> tanLambdaDBGFilesTrackleting;
    for (int i{1}; i < 12; ++i)
    {
        tanLambdaDBGFilesTrackleting.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_tanLambdaCut_%03d.root", path.data(), tanLambdaCutVariationSelectionDir.data(), i)));
    }

    std::vector<TFile *> tanlambdaVariationfiles;
    tanlambdaVariationfiles.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_tanLambdaCut_0005.root", path.data(), tanlambdawithnarrowphipath.data())));
    for (int i{1}; i < 11; ++i)
    {
        tanlambdaVariationfiles.push_back(TFile::Open(Form("%s%s/dbg_ITSVertexerCPU_tanLambdaCut_%03d.root", path.data(), tanlambdawithnarrowphipath.data(), i)));
    }
    // config
    const int stopAt = (inspEvt == -1) ? rofs->size() : inspEvt + numEvents;
    const int startAt = (inspEvt == -1) ? 0 : inspEvt;

    itsClusters.GetEntry(0);
    mcHeaderTree.GetEntry(0);

    // Hereafter: direct calls to plotting functions
    // plotClusters(startAt, stopAt, rofs, clusters, labels);
    // plotDBGCPU(&dbgCPUFile, &l2tiFile);
    plotPhiCutVariation(&l2tiFile, &dbgCPUFile, &dbgcpufileSingle505, phiDBGFilesTrackleting, phiDBGFilesSelection);
    // plotZCorrelations(phiDBGFilesTrackleting[0]);
    plotTanLambdaVariation(&l2tiFile, tanLambdaDBGFilesTrackleting, tanlambdaVariationfiles);

    return 0;
}