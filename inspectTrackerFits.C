#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TLeaf.h>
#include <THStack.h>
#include <TColor.h>

#include <iostream>
#include <array>
#endif

void inspectTrackerFits(const int lvl = 1)
{
    const float minchi{0};
    // const float maxchi{2e3};
    std::array<std::array<int, 7>, 3> maxchi;
    // maxchi[0] = {(int)2.5e4, (int)1.5e4, (int)6e3, 0, 0, 0, 0};
    // maxchi[1] = {(int)3e-5, (int)2e-3, (int)15, (int)1e2, (int)3.5e2, (int)3.5e2, (int)3e2};
    // maxchi[2] = {(int)1.5e3, (int)8e2, (int)5e2, (int)1.5e2, (int)1e2, (int)2e-6, (int)1e-9};
    maxchi[0] = {(int)5e4, (int)5e4, (int)5e4, (int)5e4, (int)5e4, (int)5e4, (int)5e4};
    maxchi[1] = {(int)1e3, (int)1e3, (int)5e1, (int)5e2, (int)1e3, (int)1e3, (int)1e3};
    maxchi[2] = {(int)5e2, (int)5e2, (int)1e3, (int)1e3, (int)1e3, (int)1e3, (int)1e3};
    std::array<float, 3> maxChiSmoother = {0.001, 0.01, 1000};
    const int nBins{300};
    const int cx{1600};
    const int cy{3200};
    const float minpt{0.};
    const float maxpt{5.};
    std::array<EColor, 7> stackColors = {kBlack, kYellow, kOrange, kRed, kBlue, kGreen, kViolet};
    bool isFake;
    auto inFile = TFile::Open("dbg_ITSTrackerCPU.root");
    std::array<TTree *, 3> fitTreesInfo;
    auto treeTrackParams = (TTree *)inFile->Get("TrackParams");
    auto treeSmoothingParams = (TTree *)inFile->Get("SmoothingParams");
    for (auto iTree{0}; iTree < 3; ++iTree)
    {
        fitTreesInfo[iTree] = (TTree *)inFile->Get(Form("Fit%dInfo", iTree));
    }

    std::array<std::array<TH1F *, 7>, 3> histChi2Valid;
    std::array<std::array<TH1F *, 7>, 3> histChi2Fake;
    std::array<TH1F *, 3> histSmoothChi2Valid;
    std::array<TH1F *, 3> histSmoothChi2Fake;
    std::array<TH1F *, 7> clusterFakeStacker;
    std::array<TCanvas *, 3> fitCanvas;
    std::array<TCanvas *, 3> smoothCanvas;
    auto stackClusters = new THStack("clustersStack", "Fake clusters;Layer");

    auto histClusterFakes = new TH1F("clustersFake", "Fake clusters;Layer", 7, -0.5, 6.5);
    auto histNFakeClusters = new TH1F("nFakeClusters", "Fake clusters per track;N_{fake}", 8, -0.5, 7.5);
    histClusterFakes->SetMinimum(0.0000000001);
    // histNFakeClusters->SetMinimum(0.0000000001);
    if (lvl > 1)
    {
        for (size_t iS{0}; iS < clusterFakeStacker.size(); ++iS)
        {
            clusterFakeStacker[iS] = new TH1F(Form("nFakeClusters%d", (int)iS), "Fake clusters;Layer", 7, -0.5, 6.5);
            clusterFakeStacker[iS]->SetLineColor(stackColors[iS]);
            clusterFakeStacker[iS]->SetFillColor(stackColors[iS]);
            stackClusters->Add(clusterFakeStacker[iS]);
        }

        for (auto iEnt{0}; iEnt < treeTrackParams->GetEntriesFast(); ++iEnt)
        {
            treeTrackParams->GetEntry(iEnt);
            isFake = treeTrackParams->GetLeaf("fake")->GetValue(0);
            if (isFake)
            {
                int nFakes = treeTrackParams->GetLeaf("nFakeClusters")->GetValue(0);
                histNFakeClusters->Fill(nFakes);
                for (auto iLayer{0}; iLayer < 7; ++iLayer)
                {
                    int cluLabel = treeTrackParams->GetLeaf(Form("clu%dLabel", iLayer))->GetValue(0);
                    if (!cluLabel)
                    {
                        clusterFakeStacker[nFakes - 1]->Fill(iLayer);
                    }
                }
            }
        }
        auto trackParamsCanvas = new TCanvas("clusterInfo", "clusterInfo", cx, cy / 2);
        trackParamsCanvas->Divide(2);
        trackParamsCanvas->cd(1);
        gPad->SetLogy();
        gPad->SetGridx();
        histNFakeClusters->Draw();
        trackParamsCanvas->cd(2);
        // gPad->SetLogy();
        stackClusters->Draw();
        gPad->Modified();
        gPad->Update();
        trackParamsCanvas->SaveAs("clustersCanvas.png");
    }

    if (lvl > 0)
    {
        for (auto iSmoothHist{0}; iSmoothHist < (int)histSmoothChi2Valid.size(); ++iSmoothHist)
        {
            histSmoothChi2Valid[iSmoothHist] = new TH1F(Form("histLayer%dvalid", iSmoothHist + 2), Form("Smoother #chi^{2} Layer %d;#chi^{2}_{Smoother}", iSmoothHist + 2), nBins, minchi, maxChiSmoother[iSmoothHist]);
            histSmoothChi2Valid[iSmoothHist]->SetLineColor(kBlue);
            histSmoothChi2Fake[iSmoothHist] = new TH1F(Form("histLayer%dfake", iSmoothHist + 2), Form("Smoother #chi^{2} Layer %d;#chi^{2}_{Smoother}", iSmoothHist + 2), nBins, minchi, maxChiSmoother[iSmoothHist]);
            histSmoothChi2Fake[iSmoothHist]->SetLineColor(kRed);
        }

        for (auto iEntry{0}; iEntry < treeSmoothingParams->GetEntriesFast(); ++iEntry)
        {
            treeSmoothingParams->GetEntry(iEntry);
            int layer = treeSmoothingParams->GetLeaf("layer")->GetValue(0);
            float schi2 = treeSmoothingParams->GetLeaf("schi2")->GetValue(0);
            bool fake = treeSmoothingParams->GetLeaf("fake")->GetValue(0);
            if (fake)
            {
                histSmoothChi2Fake[layer - 2]->Fill(schi2);
            }
            else
            {
                histSmoothChi2Valid[layer - 2]->Fill(schi2);
            }
        }

        for (auto iCanvas{0}; iCanvas < 3; ++iCanvas)
        {
            smoothCanvas[iCanvas] = new TCanvas(Form("canvasSmoothLayer%d", iCanvas + 2), Form("canvasSmoothLayer%d", iCanvas), cx, cy);
            smoothCanvas[iCanvas]->cd();
            gPad->SetLogy();
            histSmoothChi2Valid[iCanvas]->Scale(1. / histSmoothChi2Valid[iCanvas]->GetEntries());
            histSmoothChi2Fake[iCanvas]->Scale(1. / histSmoothChi2Fake[iCanvas]->GetEntries());
            histSmoothChi2Valid[iCanvas]->Draw("hist");
            histSmoothChi2Fake[iCanvas]->Draw("hist same");
            gPad->Modified();
            gPad->Update();

            smoothCanvas[iCanvas]->SaveAs(Form("SmootherChi2Layer%d", iCanvas + 2), "r");
        }
    }

    if (lvl > 2)
    {
        for (auto iFit{0}; iFit < 3; ++iFit)
        {
            for (auto iLayer{0}; iLayer < 7; ++iLayer)
            {

                histChi2Valid[iFit][iLayer] = new TH1F(Form("histSumChiValidLayer%d_%d", iLayer, iFit), Form("#chi^{2}_{predicted}, layer %d;#chi^{2}_{predicted}", iLayer), nBins, minchi, maxchi[iFit][iLayer]); //
                histChi2Valid[iFit][iLayer]->SetLineColor(kBlue);
                histChi2Fake[iFit][iLayer] = new TH1F(Form("histSumChiFakeLayer%d_%d", iLayer, iFit), Form("#chi^{2}_{predicted}, layer %d;#chi^{2}_{predicted}", iLayer), nBins, minchi, maxchi[iFit][iLayer]); //
                histChi2Fake[iFit][iLayer]->SetLineColor(kRed);
            }
        }

        float prChi2, inChi2;
        for (auto iFit{0}; iFit < 3; ++iFit)
        {
            for (auto iEnt{0}; iEnt < fitTreesInfo[iFit]->GetEntriesFast(); ++iEnt)
            {
                fitTreesInfo[iFit]->GetEntry(iEnt);
                isFake = fitTreesInfo[iFit]->GetLeaf("fake")->GetValue(0);
                for (auto iLayer{0}; iLayer < 7; ++iLayer)
                {
                    prChi2 = fitTreesInfo[iFit]->GetLeaf(Form("Layer%dchi2", iLayer))->GetValue(0);
                    inChi2 = fitTreesInfo[iFit]->GetLeaf(Form("Layer%dInChi2", iLayer))->GetValue(0);
                    if (!isFake)
                    {
                        histChi2Valid[iFit][iLayer]->Fill(prChi2);
                    }
                    else
                    {
                        histChi2Fake[iFit][iLayer]->Fill(prChi2);
                    }
                }
            }

            fitCanvas[iFit] = new TCanvas(Form("Chi2_fit%d", iFit), Form("Chi2_fit%d", iFit), cx, cy);
            fitCanvas[iFit]->Divide(2, 4);
            for (auto iHist{0}; iHist < 7; ++iHist)
            {
                fitCanvas[iFit]->cd(iHist + 1);
                gPad->SetLogy();
                // histChi2Valid[iFit][iHist]->Scale(1. / (histChi2Fake[iFit][iHist]->GetEntries() + histChi2Valid[iFit][iHist]->GetEntries()));
                // histChi2Fake[iFit][iHist]->Scale(1. / (histChi2Fake[iFit][iHist]->GetEntries() + histChi2Valid[iFit][iHist]->GetEntries()));
                histChi2Valid[iFit][iHist]->Draw("hist");
                histChi2Fake[iFit][iHist]->Draw("same hist");
            }
            fitCanvas[iFit]->SaveAs(Form("chiSquare_Fit%d.png", iFit), "r");
        }
    }
}
