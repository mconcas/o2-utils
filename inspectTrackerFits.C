#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TLeaf.h>

#include <iostream>
#include <array>
#endif

std::pair<int, int> getMainId(std::array<int, 7> &ids, std::array<int, 7> &eids)
{
    std::pair<int, int> maxOccurrenceValues = std::make_pair<int, int>(-999, -999);
    int count{0};

    for (int iCluster = 0; iCluster < 7; ++iCluster)
    {
        if (ids[iCluster] == maxOccurrenceValues.first && eids[iCluster] == maxOccurrenceValues.second)
        {
            ++count;
        }
        else
        {
            if (count != 0)
                // only in the first iteration count can be 0 at this point
                --count;
            if (count == 0)
            {
                maxOccurrenceValues.first = ids[iCluster];
                maxOccurrenceValues.second = eids[iCluster];
                count = 1;
            }
        }
    }
    return maxOccurrenceValues;
}

void inspectTrackerFits()
{
    const float minchi{0};
    const float maxchi{2e3};
    const int nBins{300};
    const int cx{1600};
    const int cy{1600};
    const float minpt{0.};
    const float maxpt{5.};

    auto inFile = TFile::Open("dbg_ITSTrackerCPU.root");

    auto treeFit1Info = (TTree *)inFile->Get("Fit1Info");
    auto treeFit2Info = (TTree *)inFile->Get("Fit2Info");
    auto treeFit3Info = (TTree *)inFile->Get("Fit3Info");

    // Layer 0
    auto histSumChiValidF2F3Layer0 = new TH1F("histSumChiValidF2F3Layer0", "Layer 0;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer0->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer0 = new TH1F("histSumChiFakeF2F3Layer0", "Layer 0;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer0->SetLineColor(kRed);

    // Layer 1
    auto histSumChiValidF2F3Layer1 = new TH1F("histSumChiValidF2F3Layer1", "Layer 1;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer1->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer1 = new TH1F("histSumChiFakeF2F3Layer1", "Layer 1;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer1->SetLineColor(kRed);

    // Layer 2
    auto histSumChiValidF2F3Layer2 = new TH1F("histSumChiValidF2F3Layer2", "Layer 2;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer2->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer2 = new TH1F("histSumChiFakeF2F3Layer2", "Layer 2;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer2->SetLineColor(kRed);

    // Layer 3
    auto histSumChiValidF2F3Layer3 = new TH1F("histSumChiValidF2F3Layer3", "Layer 3;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer3->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer3 = new TH1F("histSumChiFakeF2F3Layer3", "Layer 3;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer3->SetLineColor(kRed);

    // Layer 4
    auto histSumChiValidF2F3Layer4 = new TH1F("histSumChiValidF2F3Layer4", "Layer 4;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer4->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer4 = new TH1F("histSumChiFakeF2F3Layer4", "Layer 4;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer4->SetLineColor(kRed);

    // Layer 5
    auto histSumChiValidF2F3Layer5 = new TH1F("histSumChiValidF2F3Layer5", "Layer 5;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer5->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer5 = new TH1F("histSumChiFakeF2F3Layer5", "Layer 5;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer5->SetLineColor(kRed);

    // Layer 6
    auto histSumChiValidF2F3Layer6 = new TH1F("histSumChiValidF2F3Layer6", "Layer 6;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiValidF2F3Layer6->SetLineColor(kBlue);
    auto histSumChiFakeF2F3Layer6 = new TH1F("histSumChiFakeF2F3Layer6", "Layer 6;#chi^{2}_{F2+F3}", nBins, minchi, maxchi);
    histSumChiFakeF2F3Layer6->SetLineColor(kRed);

    // Fake clusters Layer
    auto histFakeIndices = new TH1I("histFakeIndices", "Fake clusters layers", 7, 0, 7);

    // Fake clusters number
    auto histFakeCounters = new TH1I("histFakeNumber", "Number of fake clusters", 7, 0, 7);

    std::cout << "Single treeFit2 has nEntries: " << treeFit2Info->GetEntriesFast() << std::endl;
    for (auto iEnt{0}; iEnt < treeFit2Info->GetEntriesFast(); ++iEnt)
    {
        treeFit2Info->GetEntry(iEnt);
        treeFit3Info->GetEntry(iEnt);

        float chiF2Layer0 = treeFit2Info->GetLeaf("Layer0chi2")->GetValue(0);
        float chiF3Layer0 = treeFit3Info->GetLeaf("Layer0chi2")->GetValue(0);
        float chiF2Layer1 = treeFit2Info->GetLeaf("Layer1chi2")->GetValue(0);
        float chiF3Layer1 = treeFit3Info->GetLeaf("Layer1chi2")->GetValue(0);
        float chiF2Layer2 = treeFit2Info->GetLeaf("Layer2chi2")->GetValue(0);
        float chiF3Layer2 = treeFit3Info->GetLeaf("Layer2chi2")->GetValue(0);
        float chiF2Layer3 = treeFit2Info->GetLeaf("Layer3chi2")->GetValue(0);
        float chiF3Layer3 = treeFit3Info->GetLeaf("Layer3chi2")->GetValue(0);
        float chiF2Layer4 = treeFit2Info->GetLeaf("Layer4chi2")->GetValue(0);
        float chiF3Layer4 = treeFit3Info->GetLeaf("Layer4chi2")->GetValue(0);
        float chiF2Layer5 = treeFit2Info->GetLeaf("Layer5chi2")->GetValue(0);
        float chiF3Layer5 = treeFit3Info->GetLeaf("Layer5chi2")->GetValue(0);
        float chiF2Layer6 = treeFit2Info->GetLeaf("Layer6chi2")->GetValue(0);
        float chiF3Layer6 = treeFit3Info->GetLeaf("Layer6chi2")->GetValue(0);

        std::array<int, 7> tIds = {(int)treeFit2Info->GetLeaf("Layer0tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer1tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer2tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer3tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer4tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer5tID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer6tID")->GetValue(0)};
        std::array<int, 7> eIds = {(int)treeFit2Info->GetLeaf("Layer0eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer1eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer2eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer3eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer4eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer5eID")->GetValue(0),
                                   (int)treeFit2Info->GetLeaf("Layer6eID")->GetValue(0)};

        bool isFakeF2 = treeFit2Info->GetLeaf("fake")->GetValue(0);

        if (isFakeF2)
        {
            histSumChiFakeF2F3Layer0->Fill(chiF2Layer0 + chiF3Layer0);
            histSumChiFakeF2F3Layer1->Fill(chiF2Layer1 + chiF3Layer1);
            histSumChiFakeF2F3Layer2->Fill(chiF2Layer2 + chiF3Layer2);
            histSumChiFakeF2F3Layer3->Fill(chiF2Layer3 + chiF3Layer3);
            histSumChiFakeF2F3Layer4->Fill(chiF2Layer4 + chiF3Layer4);
            histSumChiFakeF2F3Layer5->Fill(chiF2Layer5 + chiF3Layer5);
            histSumChiFakeF2F3Layer6->Fill(chiF2Layer6 + chiF3Layer6);

            int mainID = getMainId(tIds, eIds).first;
            int maineID = getMainId(tIds, eIds).second;
            int fakeCounter{0};
            for (auto iCluster{0}; iCluster < 7; ++iCluster)
            {
                if (tIds[iCluster] != mainID || eIds[iCluster] != maineID)
                {
                    histFakeIndices->AddBinContent(iCluster);
                    fakeCounter++;
                }
            }
            histFakeCounters->AddBinContent(fakeCounter);
        }
        else
        {
            histSumChiValidF2F3Layer0->Fill(chiF2Layer0 + chiF3Layer0);
            histSumChiValidF2F3Layer1->Fill(chiF2Layer1 + chiF3Layer1);
            histSumChiValidF2F3Layer2->Fill(chiF2Layer2 + chiF3Layer2);
            histSumChiValidF2F3Layer3->Fill(chiF2Layer3 + chiF3Layer3);
            histSumChiValidF2F3Layer4->Fill(chiF2Layer4 + chiF3Layer4);
            histSumChiValidF2F3Layer5->Fill(chiF2Layer5 + chiF3Layer5);
            histSumChiValidF2F3Layer6->Fill(chiF2Layer6 + chiF3Layer6);
        }
    }

    auto canvas = new TCanvas("chi2sums", "chi2sums", cy, cx);
    canvas->Divide(3, 3);
    canvas->cd(1);
    gPad->SetLogy();
    histSumChiValidF2F3Layer0->Draw();
    histSumChiFakeF2F3Layer0->Draw("same");
    canvas->cd(2);
    gPad->SetLogy();
    histSumChiValidF2F3Layer1->Draw();
    histSumChiFakeF2F3Layer1->Draw("same");
    canvas->cd(3);
    gPad->SetLogy();
    histSumChiValidF2F3Layer2->Draw();
    histSumChiFakeF2F3Layer2->Draw("same");
    canvas->cd(4);
    gPad->SetLogy();
    histSumChiValidF2F3Layer3->Draw();
    histSumChiFakeF2F3Layer3->Draw("same");
    canvas->cd(5);
    gPad->SetLogy();
    histSumChiValidF2F3Layer4->Draw();
    histSumChiFakeF2F3Layer4->Draw("same");
    canvas->cd(6);
    gPad->SetLogy();
    histSumChiValidF2F3Layer5->Draw();
    histSumChiFakeF2F3Layer5->Draw("same");
    canvas->cd(7);
    gPad->SetLogy();
    histSumChiValidF2F3Layer6->Draw();
    histSumChiFakeF2F3Layer6->Draw("same");
    canvas->cd(8);
    histFakeIndices->Draw();
    canvas->cd(9);
    gPad->SetLogy();
    histFakeCounters->Draw();

    canvas->SaveAs("fitInspectionChi2.pdf", "r");
}
