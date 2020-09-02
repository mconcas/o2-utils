#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#endif

void plotDeltas()
{
    const float dymin{-1.f};
    const float dymax{1.f};
    const float dzmin{-1.f};
    const float dzmax{1.f};
    const float minchi{0};
    const float maxchi{1e4};

    auto inFile = TFile::Open("dbg_ITSTrackerCPU.root");

    auto treeLayer4 = (TTree *)inFile->Get("layer4Deltas");
    auto fakeDYL4 = new TH1F("fakeDYL4", "Layer 4 Fake #DeltaY;#DeltaY(cm)", 300, dymin, dymax);
    fakeDYL4->SetLineColor(kRed);
    auto fakeDZL4 = new TH1F("fakeDZL4", "Layer 4 Fake #DeltaZ;#DeltaZ(cm)", 300, dymin, dymax);
    fakeDZL4->SetLineColor(kRed);
    auto trueDYL4 = new TH1F("trueDYL4", "Layer 4 True #DeltaY;#DeltaY(cm)", 300, dzmin, dzmax);
    trueDYL4->SetLineColor(kBlue);
    auto trueDZL4 = new TH1F("trueDZL4", "Layer 4 True #DeltaZ;#DeltaZ(cm)", 300, dzmin, dzmax);
    trueDZL4->SetLineColor(kBlue);
    auto fakechiout4 = new TH1F("fakechiout4", "Layer 4 Fake #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    fakechiout4->SetLineColor(kRed);
    auto truechiout4 = new TH1F("truechiout4", "Layer 4 True #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    truechiout4->SetLineColor(kBlue);
    auto fakechin4 = new TH1F("fakechin4", "Layer 4 Fake #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    fakechin4->SetLineColor(kRed);
    auto truechin4 = new TH1F("truechin4", "Layer 4 True #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    truechin4->SetLineColor(kBlue);

    auto fakeDYDZL4 = new TH2F("fakeDYDZL4", "Layer 4 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);
    auto trueDYDZL4 = new TH2F("trueDYDZL4", "Layer 4 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);

    treeLayer4->Draw("deltaY >> fakeDYL4", "fake", "goff");
    treeLayer4->Draw("deltaZ >> fakeDZL4", "fake", "goff");
    treeLayer4->Draw("deltaY >> trueDYL4", "!fake", "goff");
    treeLayer4->Draw("deltaZ >> trueDZL4", "!fake", "goff");
    treeLayer4->Draw("deltaY:deltaZ >> fakeDYDZL4", "fake", "goff");
    treeLayer4->Draw("deltaY:deltaZ >> trueDYDZL4", "!fake", "goff");
    treeLayer4->Draw("xout >> fakechiout4", "fake", "goff");
    treeLayer4->Draw("xout >> truechiout4", "!fake", "goff");
    treeLayer4->Draw("xin >> fakechin4", "fake", "goff");
    treeLayer4->Draw("xin >> truechin4", "!fake", "goff");

    auto treeLayer3 = (TTree *)inFile->Get("layer3Deltas");
    auto fakeDYL3 = new TH1F("fakeDYL3", "Layer 3 Fake #DeltaY;#DeltaY (cm)", 300, dymin, dymax);
    fakeDYL3->SetLineColor(kRed);
    auto fakeDZL3 = new TH1F("fakeDZL3", "Layer 3 Fake #DeltaZ;#DeltaZ (cm)", 300, dymin, dymax);
    fakeDZL3->SetLineColor(kRed);
    auto trueDYL3 = new TH1F("trueDYL3", "Layer 3 True #DeltaY;#DeltaY (cm)", 300, dzmin, dzmax);
    trueDYL3->SetLineColor(kBlue);
    auto trueDZL3 = new TH1F("trueDZL3", "Layer 3 True #DeltaZ;#DeltaZ (cm)", 300, dzmin, dzmax);
    trueDZL3->SetLineColor(kBlue);
    auto fakechiout3 = new TH1F("fakechiout3", "Layer 3 Fake #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    fakechiout3->SetLineColor(kRed);
    auto truechiout3 = new TH1F("truechiout3", "Layer 3 True #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    truechiout3->SetLineColor(kBlue);
    auto fakechin3 = new TH1F("fakechin3", "Layer 3 Fake #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    fakechin3->SetLineColor(kRed);
    auto truechin3 = new TH1F("truechin3", "Layer 3 True #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    truechin3->SetLineColor(kBlue);

    auto fakeDYDZL3 = new TH2F("fakeDYDZL3", "Layer 3 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);
    auto trueDYDZL3 = new TH2F("trueDYDZL3", "Layer 3 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);

    treeLayer3->Draw("deltaY >> fakeDYL3", "fake", "goff");
    treeLayer3->Draw("deltaZ >> fakeDZL3", "fake", "goff");
    treeLayer3->Draw("deltaY >> trueDYL3", "!fake", "goff");
    treeLayer3->Draw("deltaZ >> trueDZL3", "!fake", "goff");
    treeLayer3->Draw("deltaY:deltaZ >> fakeDYDZL3", "fake", "goff");
    treeLayer3->Draw("deltaY:deltaZ >> trueDYDZL3", "!fake", "goff");
    treeLayer3->Draw("xout >> fakechiout3", "fake", "goff");
    treeLayer3->Draw("xout >> truechiout3", "!fake", "goff");
    treeLayer3->Draw("xin >> fakechin3", "fake", "goff");
    treeLayer3->Draw("xin >> truechin3", "!fake", "goff");

    auto treeLayer2 = (TTree *)inFile->Get("layer2Deltas");
    auto fakeDYL2 = new TH1F("fakeDYL2", "Layer 2 Fake #DeltaY;#DeltaY (cm)", 300, dymin, dymax);
    fakeDYL2->SetLineColor(kRed);
    auto fakeDZL2 = new TH1F("fakeDZL2", "Layer 2 Fake #DeltaZ;#DeltaZ (cm)", 300, dymin, dymax);
    fakeDZL2->SetLineColor(kRed);
    auto trueDYL2 = new TH1F("trueDYL2", "Layer 2 True #DeltaY;#DeltaY (cm)", 300, dzmin, dzmax);
    trueDYL2->SetLineColor(kBlue);
    auto trueDZL2 = new TH1F("trueDZL2", "Layer 2 True #DeltaZ;#DeltaZ (cm)", 300, dzmin, dzmax);
    trueDZL2->SetLineColor(kBlue);
    auto fakechiout2 = new TH1F("fakechiout2", "Layer 2 Fake #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    fakechiout2->SetLineColor(kRed);
    auto truechiout2 = new TH1F("truechiout2", "Layer 2 True #chi^{2} Out;#chi^{2}", 300, minchi, maxchi);
    truechiout2->SetLineColor(kBlue);
    auto fakechin2 = new TH1F("fakechin2", "Layer 2 Fake #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    fakechin2->SetLineColor(kRed);
    auto truechin2 = new TH1F("truechin2", "Layer 2 True #chi^{2} In;#chi^{2}", 300, minchi, maxchi);
    truechin2->SetLineColor(kBlue);

    auto fakeDYDZL2 = new TH2F("fakeDYDZL2", "Layer 2 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);
    auto trueDYDZL2 = new TH2F("trueDYDZL2", "Layer 2 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", 300, dymin, dymax, 300, dzmin, dzmax);

    treeLayer2->Draw("deltaY >> fakeDYL2", "fake", "goff");
    treeLayer2->Draw("deltaZ >> fakeDZL2", "fake", "goff");
    treeLayer2->Draw("deltaY >> trueDYL2", "!fake", "goff");
    treeLayer2->Draw("deltaZ >> trueDZL2", "!fake", "goff");
    treeLayer2->Draw("deltaY:deltaZ >> fakeDYDZL2", "fake", "goff");
    treeLayer2->Draw("deltaY:deltaZ >> trueDYDZL2", "!fake", "goff");
    treeLayer2->Draw("xout >> fakechiout2", "fake", "goff");
    treeLayer2->Draw("xout >> truechiout2", "!fake", "goff");
    treeLayer2->Draw("xin >> fakechin2", "fake", "goff");
    treeLayer2->Draw("xin >> truechin2", "!fake", "goff");

    fakeDYL4->Scale(1. / fakeDYL4->GetEntries());
    fakeDZL4->Scale(1. / fakeDZL4->GetEntries());
    trueDYL4->Scale(1. / trueDYL4->GetEntries());
    trueDZL4->Scale(1. / trueDZL4->GetEntries());
    fakeDYL3->Scale(1. / fakeDYL3->GetEntries());
    fakeDZL3->Scale(1. / fakeDZL3->GetEntries());
    trueDYL3->Scale(1. / trueDYL3->GetEntries());
    trueDZL3->Scale(1. / trueDZL3->GetEntries());
    fakeDYL2->Scale(1. / fakeDYL2->GetEntries());
    fakeDZL2->Scale(1. / fakeDZL2->GetEntries());
    trueDYL2->Scale(1. / trueDYL2->GetEntries());
    trueDZL2->Scale(1. / trueDZL2->GetEntries());

    // Plot!
    auto canvasL4 = new TCanvas("canvasL4", "Layer 4", 1000, 1600);
    canvasL4->cd();
    canvasL4->Divide(2, 3);
    canvasL4->cd(1);
    gPad->SetLogy();
    trueDYL4->Draw("hist");
    fakeDYL4->Draw("same hist");
    canvasL4->cd(2);
    gPad->SetLogy();
    trueDZL4->Draw("hist");
    fakeDZL4->Draw("same hist");
    canvasL4->cd(3);
    trueDYDZL4->Draw("colz");
    canvasL4->cd(4);
    fakeDYDZL4->Draw("colz");
    canvasL4->cd(5);
    gPad->SetLogy();
    truechiout4->Draw();
    fakechiout4->Draw("same");
    canvasL4->cd(6);
    gPad->SetLogy();
    truechin4->Draw();
    fakechin4->Draw("same");

    auto canvasL3 = new TCanvas("canvasL3", "Layer 3", 1000, 1600);
    canvasL3->cd();
    canvasL3->Divide(2, 3);
    canvasL3->cd(1);
    gPad->SetLogy();
    trueDYL3->Draw("hist");
    fakeDYL3->Draw("same hist");
    canvasL3->cd(2);
    gPad->SetLogy();
    trueDZL3->Draw("hist");
    fakeDZL3->Draw("same hist");
    canvasL3->cd(3);
    trueDYDZL3->Draw("colz");
    canvasL3->cd(4);
    fakeDYDZL3->Draw("colz");
    canvasL3->cd(5);
    gPad->SetLogy();
    truechiout3->Draw();
    fakechiout3->Draw("same");
    canvasL3->cd(6);
    gPad->SetLogy();
    truechin3->Draw();
    fakechin3->Draw("same");

    auto canvasL2 = new TCanvas("canvasL2", "Layer 2", 1000, 1600);
    canvasL2->cd();
    canvasL2->Divide(2, 3);
    canvasL2->cd(1);
    gPad->SetLogy();
    trueDYL2->Draw("hist");
    fakeDYL2->Draw("same hist");
    canvasL2->cd(2);
    gPad->SetLogy();
    trueDZL2->Draw("hist");
    fakeDZL2->Draw("same hist");
    canvasL2->cd(3);
    trueDYDZL2->Draw("colz");
    canvasL2->cd(4);
    fakeDYDZL2->Draw("colz");
    canvasL2->cd(5);
    gPad->SetLogy();
    truechiout2->Draw();
    fakechiout2->Draw("same");
    canvasL2->cd(6);
    gPad->SetLogy();
    truechin2->Draw();
    fakechin2->Draw("same");

    canvasL4->SaveAs("layer4.png");
    canvasL3->SaveAs("layer3.png");
    canvasL2->SaveAs("layer2.png");
}