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
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"

#include <vector>

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

void postProcess(const std::string fileName = "/data1/ruben/pbpbVtx/dbg_ITSVertexerCPU.root",
                 const std::string fileNameLabels = "/data1/ruben/pbpbVtx/label2Track0.root")
{
    gStyle->SetOptStat(0);
    // ROOT::EnableImplicitMT();
    auto vertInfo = ROOT::RDataFrame("verticesInfo", fileName);
    auto histEvtIds = vertInfo.Histo1D({"EventIDs", "Reconstructed vertices in ROFrames; event ID; reconstructed", 151, -1.5, 149.5}, "eventId");
    auto canvasEvtID = new TCanvas("EvtID", "EventID", 800, 600);
    histEvtIds->SetFillColor(kBrownC);
    histEvtIds->SetLineColor(kBrownC);
    histEvtIds->DrawClone();
    auto lEvtIds = new TLegend(0.70, 0.2, 0.98, 0.33);
    lEvtIds->SetHeader("150 events PbPb MB, 545 ROframes", "C");
    lEvtIds->AddEntry("EventIDs", "Reconstructed EvtID", "f");
    lEvtIds->SetBorderSize(1);
    lEvtIds->SetTextSize(0.022);
    lEvtIds->Draw();
    canvasEvtID->SaveAs("/home/mconcas/cernbox/thesis_pictures/eventsid.png", "r");

    auto histPurity = vertInfo.Histo1D({"Purity", "Purity of reconstructed vertices; #it{purity}; counts", 300, 0.9f, 1.f}, "purity");
    histPurity->SetFillColor(kRedC);
    histPurity->SetLineColor(kRedC);
    auto canvasPurity = new TCanvas("Purity", "Purity", 800, 600);
    histPurity->DrawClone();
    auto lPurity = new TLegend(0.2, 0.7, 0.55, 0.83);
    lPurity->SetHeader("150 events PbPb MB, 545 ROframes", "C");
    lPurity->AddEntry("Purity", "#it{purity} = #frac{validated lines}{used lines}", "f");
    lPurity->SetBorderSize(1);
    lPurity->SetTextSize(0.022);
    lPurity->Draw();
    canvasPurity->SaveAs("/home/mconcas/cernbox/thesis_pictures/verticespurity.png", "r");

    // residuals
    gStyle->SetOptStat(1);
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
    auto hResX = new TH1F("resX", "Residuals X; #DeltaX (cm); counts", 300u, -0.014f, 0.014f);
    auto hResY = new TH1F("resY", "Residuals Y; #DeltaY (cm); counts", 300u, -0.014f, 0.014f);
    auto hResZ = new TH1F("resZ", "Residuals Z; #DeltaZ (cm); counts", 300u, -0.014f, 0.014f);
    auto vResX = std::vector<double>{};
    auto vResY = std::vector<double>{};
    auto vResZ = std::vector<double>{};
    auto vContribs = std::vector<double>{};
    auto vPurity = std::vector<double>{};

    labels2Tracks.Foreach(funzione, {"MCLabels", "Tracks"});
    auto funzione2 = [&](int evtId, float x, float y, float z, int contributors, float purity) {
        auto it = umap.find(evtId);
        if (it != umap.end())
        {
            hResX->Fill(x - it->second[0]);
            vResX.push_back((double)TMath::Abs(x - it->second[0]));
            hResY->Fill(y - it->second[1]);
            vResY.push_back((double)TMath::Abs(y - it->second[1]));
            hResZ->Fill(z - it->second[2]);
            vResZ.push_back((double)TMath::Abs(z - it->second[2]));
            vContribs.push_back((double)contributors);
            vPurity.push_back((double)purity);
            umap.erase(it);
        }
    };

    vertInfo.Foreach(funzione2, {"eventId", "xCoord", "yCoord", "zCoord", "size", "purity"});
    auto canvasX = new TCanvas("resX", "resX", 800, 600);
    hResX->SetFillColor(kOrangeC);
    hResX->SetLineColor(kOrangeC);
    hResX->Draw();
    auto lResX = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResX->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hResX->GetEntries()), "C");
    lResX->AddEntry("resX", "residuals", "f");
    lResX->SetTextSize(0.022);
    lResX->SetBorderSize(1);
    lResX->Draw();
    auto canvasY = new TCanvas("resY", "resY", 800, 600);
    hResY->SetFillColor(kPurpleC);
    hResY->SetLineColor(kPurpleC);
    hResY->Draw();
    auto lResY = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResY->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hResY->GetEntries()), "C");
    lResY->AddEntry("resY", "residuals", "f");
    lResY->SetTextSize(0.022);
    lResY->SetBorderSize(1);
    lResY->Draw();
    auto canvasZ = new TCanvas("resZ", "resZ", 800, 600);
    hResZ->SetFillColor(kBlueC);
    hResZ->SetLineColor(kBlueC);
    hResZ->Draw();
    auto lResZ = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResZ->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hResZ->GetEntries()), "C");
    lResZ->AddEntry("resZ", "residuals", "f");
    lResZ->SetTextSize(0.022);
    lResZ->SetBorderSize(1);
    lResZ->Draw();
    canvasX->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResX.png", "r");
    canvasY->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResY.png", "r");
    canvasZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/vertexResZ.png", "r");

    // Residuals vs nContribs
    gStyle->SetOptStat(0);
    auto hContVsResx = new TH2F("hContVsResx", "Residuals X vs contributors; contributors; #DeltaX(cm) ", 150, 0, 7000, 150, 0.f, 0.014f);
    hContVsResx->FillN((int)vContribs.size(), vContribs.data(), vResX.data(), nullptr);
    auto hContVsResy = new TH2F("hContVsResy", "Residuals Y vs contributors; contributors; #DeltaY(cm) ", 150, 0, 7000, 150, 0.f, 0.014f);
    hContVsResy->FillN((int)vContribs.size(), vContribs.data(), vResY.data(), nullptr);
    auto hContVsResz = new TH2F("hContVsResz", "Residuals Z vs contributors; contributors; #DeltaZ(cm) ", 150, 0, 7000, 150, 0.f, 0.014f);
    hContVsResz->FillN((int)vContribs.size(), vContribs.data(), vResZ.data(), nullptr);
    auto canvasResidualsContX = new TCanvas("resVsContX", "resVsContX", 800, 600);
    canvasResidualsContX->SetLogy();
    hContVsResx->Draw();
    hContVsResx->SetMarkerColor(kOrangeC);
    hContVsResx->SetMarkerStyle(8);
    auto lResXvsCont = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResXvsCont->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hContVsResx->GetEntries()), "C");
    lResXvsCont->AddEntry("hContVsResx", "residuals", "p");
    lResXvsCont->SetTextSize(0.022);
    lResXvsCont->SetBorderSize(1);
    lResXvsCont->Draw();
    auto canvasResidualsContY = new TCanvas("resVsContY", "resVsContY", 800, 600);
    canvasResidualsContY->SetLogy();
    hContVsResy->Draw();
    hContVsResy->SetMarkerColor(kPurpleC);
    hContVsResy->SetMarkerStyle(8);
    auto lResYvsCont = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResYvsCont->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hContVsResy->GetEntries()), "C");
    lResYvsCont->AddEntry("hContVsResy", "residuals", "p");
    lResYvsCont->SetTextSize(0.022);
    lResYvsCont->SetBorderSize(1);
    lResYvsCont->Draw();
    auto canvasResidualsContZ = new TCanvas("resVsContZ", "resVsContZ", 800, 600);
    canvasResidualsContZ->SetLogy();
    hContVsResz->Draw();
    hContVsResz->SetMarkerColor(kBlueC);
    hContVsResz->SetMarkerStyle(8);
    auto lResZvsCont = new TLegend(0.60, 0.58, 0.86, 0.72);
    lResZvsCont->SetHeader(Form("Residuals: %d reconstructed vertices", (int)hContVsResz->GetEntries()), "C");
    lResZvsCont->AddEntry("hContVsResz", "residuals", "p");
    lResZvsCont->SetTextSize(0.022);
    lResZvsCont->SetBorderSize(1);
    lResZvsCont->Draw();
    canvasResidualsContX->SaveAs("/home/mconcas/cernbox/thesis_pictures/residualsVsContributorsX.png", "r");
    canvasResidualsContY->SaveAs("/home/mconcas/cernbox/thesis_pictures/residualsVsContributorsY.png", "r");
    canvasResidualsContZ->SaveAs("/home/mconcas/cernbox/thesis_pictures/residualsVsContributorsZ.png", "r");

    // Purity vs nContribs
    auto hPurityVsConts = new TH2F("hPurityVsConts", "Vertices purity vs contributors; contributors; #it{purity}", 150, 0, 7000, 150, 0.9f, 1.f);
    hPurityVsConts->FillN((int)vContribs.size(), vContribs.data(), vPurity.data(), nullptr);
    auto canvasPurityVsContribs = new TCanvas("purVsConts", "purVsConts", 800, 600);
    hPurityVsConts->SetMarkerColor(kGreenC);
    hPurityVsConts->SetMarkerStyle(23);
    hPurityVsConts->Draw();
    auto lPurityvsCont = new TLegend(0.60, 0.58, 0.86, 0.72);
    lPurityvsCont->SetHeader(Form("Purity: %d reconstructed vertices", (int)hPurityVsConts->GetEntries()), "C");
    lPurityvsCont->AddEntry("hPurityVsConts", "#it{purity} = #frac{validated lines}{used lines}", "p");
    lPurityvsCont->SetTextSize(0.022);
    lPurityvsCont->SetBorderSize(1);
    lPurityvsCont->Draw();
    canvasPurityVsContribs->SaveAs("/home/mconcas/cernbox/thesis_pictures/purityVsContributors.png", "r");
}

#endif