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
    const float maxchi{2e3};
    const int nBins{300};
    const int cx{1600};
    const int cy{1200};
    const float minpt{0.};
    const float maxpt{5.};

    auto inFile = TFile::Open("dbg_ITSTrackerCPU.root");

    auto treeLayer4 = (TTree *)inFile->Get("layer4Deltas");
    auto fakeDYL4 = new TH1F("fakeDYL4", "Layer 4 Fake #DeltaY;#DeltaY(cm)", nBins, dymin, dymax);
    fakeDYL4->SetLineColor(kRed);
    auto fakeDZL4 = new TH1F("fakeDZL4", "Layer 4 Fake #DeltaZ;#DeltaZ(cm)", nBins, dymin, dymax);
    fakeDZL4->SetLineColor(kRed);
    auto trueDYL4 = new TH1F("trueDYL4", "Layer 4 True #DeltaY;#DeltaY(cm)", nBins, dzmin, dzmax);
    trueDYL4->SetLineColor(kBlue);
    auto trueDZL4 = new TH1F("trueDZL4", "Layer 4 True #DeltaZ;#DeltaZ(cm)", nBins, dzmin, dzmax);
    trueDZL4->SetLineColor(kBlue);
    auto fakechiout4 = new TH1F("fakechiout4", "Layer 4 Fake;#chi^{2}_{Out}", nBins, minchi, maxchi);
    fakechiout4->SetLineColor(kRed);
    auto truechiout4 = new TH1F("truechiout4", "Layer 4 True;#chi^{2}_{Out}", nBins, minchi, maxchi);
    truechiout4->SetLineColor(kBlue);
    auto fakechin4 = new TH1F("fakechin4", "Layer 4 Fake;#chi^{2}_{In};", nBins, minchi, maxchi / 3);
    fakechin4->SetLineColor(kRed);
    auto truechin4 = new TH1F("truechin4", "Layer 4 True;#chi^{2}_{In};", nBins, minchi, maxchi / 3);
    truechin4->SetLineColor(kBlue);
    auto fakechi4 = new TH1F("fakechi4", "Layer 4 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    fakechi4->SetLineColor(kRed);
    auto truechi4 = new TH1F("truechi4", "Layer 4 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    truechi4->SetLineColor(kBlue);
    auto validInTrackMatchingCluster4 = new TH1F("validInTrackMatchingCluster4",
                                                 "Layer 4 valid In track Matching cluster;#chi^{2}_{In}",
                                                 nBins, minchi, maxchi);
    validInTrackMatchingCluster4->SetLineColor(kBlue);
    auto validInTrackMismatchingCluster4 = new TH1F("validInTrackMismatchingCluster4",
                                                    "Layer 4 valid In track Mismatching cluster;#chi^{2}_{In}",
                                                    nBins, minchi, maxchi);
    validInTrackMismatchingCluster4->SetLineColor(kRed);
    auto validOutTrackMatchingCluster4 = new TH1F("validOutTrackMatchingCluster4",
                                                  "Layer 4 valid Out track Matching cluster;#chi^{2}_{Out}",
                                                  nBins, minchi, maxchi);
    validOutTrackMatchingCluster4->SetLineColor(kBlue);
    auto validOutTrackMismatchingCluster4 = new TH1F("validOutTrackMismatchingCluster4",
                                                     "Layer 4 valid Out track Mismatching cluster;#chi^{2}_{Out}",
                                                     nBins, minchi, maxchi);
    validOutTrackMismatchingCluster4->SetLineColor(kRed);

    auto fakeDYDZL4 = new TH2F("fakeDYDZL4", "Layer 4 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);
    auto trueDYDZL4 = new TH2F("trueDYDZL4", "Layer 4 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);

    auto fakechinchiout4 = new TH2F("fakechinchiout4", "Layer 4 Fake;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    auto truechinchiout4 = new TH2F("truechinchiout4", "Layer 4 True;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    truechinchiout4->SetLineColor(kBlue);
    fakechinchiout4->SetLineColor(kRed);

    treeLayer4->Draw("deltaY >> fakeDYL4", "fake", "goff");
    treeLayer4->Draw("deltaZ >> fakeDZL4", "fake", "goff");
    treeLayer4->Draw("deltaY >> trueDYL4", "!fake", "goff");
    treeLayer4->Draw("deltaZ >> trueDZL4", "!fake", "goff");
    treeLayer4->Draw("deltaY:deltaZ >> fakeDYDZL4", "fake", "goff");
    treeLayer4->Draw("deltaY:deltaZ >> trueDYDZL4", "!fake", "goff");
    treeLayer4->Draw("xin:xout >> fakechinchiout4", "fake", "goff");
    treeLayer4->Draw("xin:xout >> truechinchiout4", "!fake", "goff");
    treeLayer4->Draw("xout >> fakechiout4", "fake", "goff");
    treeLayer4->Draw("xout >> truechiout4", "!fake", "goff");
    treeLayer4->Draw("xin >> fakechin4", "fake", "goff");
    treeLayer4->Draw("xin >> truechin4", "!fake", "goff");
    treeLayer4->Draw("xin+xout >> fakechi4", "fake", "goff");
    treeLayer4->Draw("xin+xout >> truechi4", "!fake", "goff");
    treeLayer4->Draw("xin >> validInTrackMatchingCluster4", "compTrIn && compClIn", "goff");
    treeLayer4->Draw("xin >> validInTrackMismatchingCluster4", "compTrIn && !compClIn", "goff");
    treeLayer4->Draw("xout >> validOutTrackMatchingCluster4", "compTrOut && compClOut", "goff");
    treeLayer4->Draw("xout >> validOutTrackMismatchingCluster4", "compTrOut && !compClOut", "goff");

    auto treeLayer3 = (TTree *)inFile->Get("layer3Deltas");
    auto fakeDYL3 = new TH1F("fakeDYL3", "Layer 3 Fake #DeltaY;#DeltaY (cm)", nBins, dymin, dymax);
    fakeDYL3->SetLineColor(kRed);
    auto fakeDZL3 = new TH1F("fakeDZL3", "Layer 3 Fake #DeltaZ;#DeltaZ (cm)", nBins, dymin, dymax);
    fakeDZL3->SetLineColor(kRed);
    auto trueDYL3 = new TH1F("trueDYL3", "Layer 3 True #DeltaY;#DeltaY (cm)", nBins, dzmin, dzmax);
    trueDYL3->SetLineColor(kBlue);
    auto trueDZL3 = new TH1F("trueDZL3", "Layer 3 True #DeltaZ;#DeltaZ (cm)", nBins, dzmin, dzmax);
    trueDZL3->SetLineColor(kBlue);
    auto fakechiout3 = new TH1F("fakechiout3", "Layer 3 Fake;#chi^{2}_{Out}", nBins, minchi, maxchi);
    fakechiout3->SetLineColor(kRed);
    auto truechiout3 = new TH1F("truechiout3", "Layer 3 True;#chi^{2}_{Out}", nBins, minchi, maxchi);
    truechiout3->SetLineColor(kBlue);
    auto fakechin3 = new TH1F("fakechin3", "Layer 3 Fake;#chi^{2}_{In}", nBins, minchi, maxchi / 3);
    fakechin3->SetLineColor(kRed);
    auto truechin3 = new TH1F("truechin3", "Layer 3 True;#chi^{2}_{In}", nBins, minchi, maxchi / 3);
    truechin3->SetLineColor(kBlue);
    auto fakechi3 = new TH1F("fakechi3", "Layer 3 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    fakechi3->SetLineColor(kRed);
    auto truechi3 = new TH1F("truechi3", "Layer 3 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    auto validInTrackMatchingCluster3 = new TH1F("validInTrackMatchingCluster3",
                                                 "Layer 3 valid In track Matching cluster;#chi^{2}_{In}",
                                                 nBins, minchi, maxchi);
    validInTrackMatchingCluster3->SetLineColor(kBlue);
    auto validInTrackMismatchingCluster3 = new TH1F("validInTrackMismatchingCluster3",
                                                    "Layer 3 valid In track Mismatching cluster;#chi^{2}_{In}",
                                                    nBins, minchi, maxchi);
    validInTrackMismatchingCluster3->SetLineColor(kRed);
    auto validOutTrackMatchingCluster3 = new TH1F("validOutTrackMatchingCluster3",
                                                  "Layer 3 valid Out track Matching cluster;#chi^{2}_{Out}",
                                                  nBins, minchi, maxchi);
    validOutTrackMatchingCluster3->SetLineColor(kBlue);
    auto validOutTrackMismatchingCluster3 = new TH1F("validOutTrackMismatchingCluster3",
                                                     "Layer 3 valid Out track Mismatching cluster;#chi^{2}_{Out}",
                                                     nBins, minchi, maxchi);
    validOutTrackMismatchingCluster3->SetLineColor(kRed);
    truechi3->SetLineColor(kBlue);

    auto fakeDYDZL3 = new TH2F("fakeDYDZL3", "Layer 3 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);
    auto trueDYDZL3 = new TH2F("trueDYDZL3", "Layer 3 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);

    auto fakechinchiout3 = new TH2F("fakechinchiout3", "Layer 3 Fake;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    auto truechinchiout3 = new TH2F("truechinchiout3", "Layer 3 True;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    truechinchiout3->SetLineColor(kBlue);
    fakechinchiout3->SetLineColor(kRed);

    treeLayer3->Draw("deltaY >> fakeDYL3", "fake", "goff");
    treeLayer3->Draw("deltaZ >> fakeDZL3", "fake", "goff");
    treeLayer3->Draw("deltaY >> trueDYL3", "!fake", "goff");
    treeLayer3->Draw("deltaZ >> trueDZL3", "!fake", "goff");
    treeLayer3->Draw("deltaY:deltaZ >> fakeDYDZL3", "fake", "goff");
    treeLayer3->Draw("deltaY:deltaZ >> trueDYDZL3", "!fake", "goff");
    treeLayer3->Draw("xin:xout >> fakechinchiout3", "fake", "goff");
    treeLayer3->Draw("xin:xout >> truechinchiout3", "!fake", "goff");
    treeLayer3->Draw("xout >> fakechiout3", "fake", "goff");
    treeLayer3->Draw("xout >> truechiout3", "!fake", "goff");
    treeLayer3->Draw("xin >> fakechin3", "fake", "goff");
    treeLayer3->Draw("xin >> truechin3", "!fake", "goff");
    treeLayer3->Draw("xin+xout >> fakechi3", "fake", "goff");
    treeLayer3->Draw("xin+xout >> truechi3", "!fake", "goff");
    treeLayer3->Draw("xin >> validInTrackMatchingCluster3", "compTrIn && compClIn", "goff");
    treeLayer3->Draw("xin >> validInTrackMismatchingCluster3", "compTrIn && !compClIn", "goff");
    treeLayer3->Draw("xout >> validOutTrackMatchingCluster3", "compTrOut && compClOut", "goff");
    treeLayer3->Draw("xout >> validOutTrackMismatchingCluster3", "compTrOut && !compClOut", "goff");

    auto treeLayer2 = (TTree *)inFile->Get("layer2Deltas");
    auto fakeDYL2 = new TH1F("fakeDYL2", "Layer 2 Fake #DeltaY;#DeltaY (cm)", nBins, dymin, dymax);
    fakeDYL2->SetLineColor(kRed);
    auto fakeDZL2 = new TH1F("fakeDZL2", "Layer 2 Fake #DeltaZ;#DeltaZ (cm)", nBins, dymin, dymax);
    fakeDZL2->SetLineColor(kRed);
    auto trueDYL2 = new TH1F("trueDYL2", "Layer 2 True #DeltaY;#DeltaY (cm)", nBins, dzmin, dzmax);
    trueDYL2->SetLineColor(kBlue);
    auto trueDZL2 = new TH1F("trueDZL2", "Layer 2 True #DeltaZ;#DeltaZ (cm)", nBins, dzmin, dzmax);
    trueDZL2->SetLineColor(kBlue);
    auto fakechiout2 = new TH1F("fakechiout2", "Layer 2 Fake;#chi^{2}_{Out}", nBins, minchi, maxchi);
    fakechiout2->SetLineColor(kRed);
    auto truechiout2 = new TH1F("truechiout2", "Layer 2 True;#chi^{2}_{Out}", nBins, minchi, maxchi);
    truechiout2->SetLineColor(kBlue);
    auto fakechin2 = new TH1F("fakechin2", "Layer 2 Fake;#chi^{2}_{In}", nBins, minchi, maxchi / 3);
    fakechin2->SetLineColor(kRed);
    auto truechin2 = new TH1F("truechin2", "Layer 2 True;#chi^{2}_{In}", nBins, minchi, maxchi / 3);
    truechin2->SetLineColor(kBlue);
    auto fakechi2 = new TH1F("fakechi2", "Layer 2 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    fakechi2->SetLineColor(kRed);
    auto truechi2 = new TH1F("truechi2", "Layer 2 True;#chi^{2}_{In}+#chi^{2}_{Out};", nBins, minchi, maxchi);
    truechi2->SetLineColor(kBlue);
    auto validInTrackMatchingCluster2 = new TH1F("validInTrackMatchingCluster2",
                                                 "Layer 2 valid In track Matching cluster;#chi^{2}_{In}",
                                                 nBins, minchi, maxchi);
    validInTrackMatchingCluster2->SetLineColor(kBlue);
    auto validInTrackMismatchingCluster2 = new TH1F("validInTrackMismatchingCluster2",
                                                    "Layer 2 valid In track Mismatching cluster;#chi^{2}_{In}",
                                                    nBins, minchi, maxchi);
    validInTrackMismatchingCluster2->SetLineColor(kRed);
    auto validOutTrackMatchingCluster2 = new TH1F("validOutTrackMatchingCluster2",
                                                  "Layer 2 valid Out track Matching cluster;#chi^{2}_{Out}",
                                                  nBins, minchi, maxchi);
    validOutTrackMatchingCluster2->SetLineColor(kBlue);
    auto validOutTrackMismatchingCluster2 = new TH1F("validOutTrackMismatchingCluster2",
                                                     "Layer 2 valid Out track Mismatching cluster;#chi^{2}_{Out}",
                                                     nBins, minchi, maxchi);
    validOutTrackMismatchingCluster2->SetLineColor(kRed);
    truechi2->SetLineColor(kBlue);

    auto fakeDYDZL2 = new TH2F("fakeDYDZL2", "Layer 2 Fake #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);
    auto trueDYDZL2 = new TH2F("trueDYDZL2", "Layer 2 True #DeltaY#DeltaZ;#DeltaY (cm);#DeltaZ (cm)", nBins, dymin, dymax, nBins, dzmin, dzmax);

    auto fakechinchiout2 = new TH2F("fakechinchiout2", "Layer 2 Fake;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    auto truechinchiout2 = new TH2F("truechinchiout2", "Layer 2 True;#chi^{2}_{in};#chi^{2}_{Out}", nBins, minchi / 100, maxchi / 100, nBins, minchi / 100, maxchi / 100);
    truechinchiout2->SetLineColor(kBlue);
    fakechinchiout2->SetLineColor(kRed);

    treeLayer2->Draw("deltaY >> fakeDYL2", "fake", "goff");
    treeLayer2->Draw("deltaZ >> fakeDZL2", "fake", "goff");
    treeLayer2->Draw("deltaY >> trueDYL2", "!fake", "goff");
    treeLayer2->Draw("deltaZ >> trueDZL2", "!fake", "goff");
    treeLayer2->Draw("deltaY:deltaZ >> fakeDYDZL2", "fake", "goff");
    treeLayer2->Draw("deltaY:deltaZ >> trueDYDZL2", "!fake", "goff");
    treeLayer2->Draw("xin:xout >> fakechinchiout2", "fake", "goff");
    treeLayer2->Draw("xin:xout >> truechinchiout2", "!fake", "goff");
    treeLayer2->Draw("xout >> fakechiout2", "fake", "goff");
    treeLayer2->Draw("xout >> truechiout2", "!fake", "goff");
    treeLayer2->Draw("xin >> fakechin2", "fake", "goff");
    treeLayer2->Draw("xin >> truechin2", "!fake", "goff");
    treeLayer2->Draw("xin+xout >> fakechi2", "fake", "goff");
    treeLayer2->Draw("xin+xout >> truechi2", "!fake", "goff");
    treeLayer2->Draw("xin >> validInTrackMatchingCluster2", "compTrIn && compClIn", "goff");
    treeLayer2->Draw("xin >> validInTrackMismatchingCluster2", "compTrIn && !compClIn", "goff");
    treeLayer2->Draw("xout >> validOutTrackMatchingCluster2", "compTrOut && compClOut", "goff");
    treeLayer2->Draw("xout >> validOutTrackMismatchingCluster2", "compTrOut && !compClOut", "goff");

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

    fakechiout4->Scale(1. / fakechiout4->GetEntries());
    truechiout4->Scale(1. / truechiout4->GetEntries());
    fakechin4->Scale(1. / fakechin4->GetEntries());
    truechin4->Scale(1. / truechin4->GetEntries());
    fakechiout3->Scale(1. / fakechiout3->GetEntries());
    truechiout3->Scale(1. / truechiout3->GetEntries());
    fakechin3->Scale(1. / fakechin3->GetEntries());
    truechin3->Scale(1. / truechin3->GetEntries());
    fakechiout2->Scale(1. / fakechiout2->GetEntries());
    truechiout2->Scale(1. / truechiout2->GetEntries());
    fakechin2->Scale(1. / fakechin2->GetEntries());
    truechin2->Scale(1. / truechin2->GetEntries());

    fakechi4->Scale(1. / fakechi4->GetEntries());
    truechi4->Scale(1. / fakechi4->GetEntries());
    fakechi3->Scale(1. / fakechi3->GetEntries());
    truechi3->Scale(1. / truechi3->GetEntries());
    fakechi2->Scale(1. / truechi2->GetEntries());
    truechi2->Scale(1. / truechi2->GetEntries());

    // ------
    auto trackTree = (TTree *)inFile->Get("TrackParams");
    auto validTracksChi2 = new TH1F("validTracksChi2", "Valid final tracks #chi^{2};#chi^{2}",
                                    nBins, minchi, maxchi);
    validTracksChi2->SetLineColor(kBlue);
    auto fakeTracksChi2 = new TH1F("fakeTracksChi2", "Fake final tracks #chi^{2};#chi^{2}",
                                   nBins, minchi, maxchi);
    fakeTracksChi2->SetLineColor(kRed);

    auto validTracksPt = new TH1F("validTracksPt", "Valid final tracks #it{p}_{T};#it{p}_{T}",
                                  nBins, minpt, maxpt);
    validTracksPt->SetLineColor(kBlue);
    auto fakeTracksPt = new TH1F("fakeTracksPt", "Fake final tracks #it{p}_{T};#it{p}_{T}",
                                 nBins, minpt, maxpt);
    fakeTracksPt->SetLineColor(kRed);

    trackTree->Draw("chi2>>validTracksChi2", "!fake", "goff");
    trackTree->Draw("chi2>>fakeTracksChi2", "fake", "goff");
    trackTree->Draw("pt>>validTracksPt", "!fake", "goff");
    trackTree->Draw("pt>>fakeTracksPt", "fake", "goff");

    validTracksChi2->Scale(1. / (fakeTracksChi2->GetEntries() + fakeTracksChi2->GetEntries()));
    fakeTracksChi2->Scale(1. / (fakeTracksChi2->GetEntries() + fakeTracksChi2->GetEntries()));
    validTracksPt->Scale(1. / (fakeTracksPt->GetEntries() + validTracksPt->GetEntries()));
    fakeTracksPt->Scale(1. / (fakeTracksPt->GetEntries() + validTracksPt->GetEntries()));

    // Plot!
    auto canvasL4 = new TCanvas("canvasL4", "Layer 4", cx, cy);
    auto canvasL4_chi2 = new TCanvas("canvasL4_chi2", "Layer 4", cx, cy);
    canvasL4->cd();
    canvasL4->Divide(2, 2);
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

    canvasL4_chi2->Divide(2, 2);
    canvasL4_chi2->cd(1);
    gPad->SetLogy();
    truechiout4->Draw("hist");
    fakechiout4->Draw("same hist");
    canvasL4_chi2->cd(2);
    gPad->SetLogy();
    truechin4->Draw("hist");
    fakechin4->Draw("same hist");
    canvasL4_chi2->cd(3);
    fakechinchiout4->Draw("colz");
    canvasL4_chi2->cd(4);
    truechinchiout4->Draw("colz");

    auto canvasL3 = new TCanvas("canvasL3", "Layer 3", cx, cy);
    auto canvasL3_chi2 = new TCanvas("canvasL3_chi2", "Layer 3", cx, cy);
    canvasL3->cd();
    canvasL3->Divide(2, 2);
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

    canvasL3_chi2->Divide(2, 2);
    canvasL3_chi2->cd(1);
    gPad->SetLogy();
    truechiout3->Draw("hist");
    fakechiout3->Draw("same hist");
    canvasL3_chi2->cd(2);
    gPad->SetLogy();
    truechin3->Draw("hist");
    fakechin3->Draw("same hist");
    canvasL3_chi2->cd(3);
    fakechinchiout3->Draw("colz");
    canvasL3_chi2->cd(4);
    truechinchiout3->Draw("colz");

    auto canvasL2 = new TCanvas("canvasL2", "Layer 2", cx, cy);
    auto canvasL2_chi2 = new TCanvas("canvasL2_chi2", "Layer 2", cx, cy);
    canvasL2->cd();
    canvasL2->Divide(2, 2);
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

    canvasL2_chi2->Divide(2, 2);
    canvasL2_chi2->cd(1);
    gPad->SetLogy();
    truechiout2->Draw("hist");
    fakechiout2->Draw("same hist");
    canvasL2_chi2->cd(2);
    gPad->SetLogy();
    truechin2->Draw("hist");
    fakechin2->Draw("same hist");
    canvasL2_chi2->cd(3);
    fakechinchiout2->Draw("colz");
    canvasL2_chi2->cd(4);
    truechinchiout2->Draw("colz");

    auto canvasSumChi2 = new TCanvas("canvasSumChi2", "Chi2 Sum", cx, cy / 2);
    canvasSumChi2->Divide(3);
    canvasSumChi2->cd(1);
    gPad->SetLogy();
    fakechi4->Draw("hist");
    truechi4->Draw("same hist");
    canvasSumChi2->cd(2);
    gPad->SetLogy();
    fakechi3->Draw("hist");
    truechi3->Draw("same hist");
    canvasSumChi2->cd(3);
    gPad->SetLogy();
    fakechi2->Draw("hist");
    truechi2->Draw("same hist");

    auto canvasL4Matches_chi2 = new TCanvas("canvasL4Matches_chi2", "Layer 4", cx, cy / 2);
    canvasL4Matches_chi2->cd();
    canvasL4Matches_chi2->Divide(2);
    canvasL4Matches_chi2->cd(1);
    gPad->SetLogy();
    validInTrackMatchingCluster4->Draw("hist");
    validInTrackMismatchingCluster4->Draw("same hist");
    canvasL4Matches_chi2->cd(2);
    gPad->SetLogy();
    validOutTrackMatchingCluster4->Draw("hist");
    validOutTrackMismatchingCluster4->Draw("same hist");

    auto canvasL3Matches_chi2 = new TCanvas("canvasL3Matches_chi2", "Layer 3", cx, cy / 2);
    canvasL3Matches_chi2->cd();
    canvasL3Matches_chi2->Divide(2);
    canvasL3Matches_chi2->cd(1);
    gPad->SetLogy();
    validInTrackMatchingCluster3->Draw("hist");
    validInTrackMismatchingCluster3->Draw("same hist");
    canvasL3Matches_chi2->cd(2);
    gPad->SetLogy();
    validOutTrackMatchingCluster3->Draw("hist");
    validOutTrackMismatchingCluster3->Draw("same hist");

    auto canvasL2Matches_chi2 = new TCanvas("canvasL2Matches_chi2", "Layer 2", cx, cy / 2);
    canvasL2Matches_chi2->cd();
    canvasL2Matches_chi2->Divide(2);
    canvasL2Matches_chi2->cd(1);
    gPad->SetLogy();
    validInTrackMatchingCluster2->Draw("hist");
    validInTrackMismatchingCluster2->Draw("same hist");
    canvasL2Matches_chi2->cd(2);
    gPad->SetLogy();
    validOutTrackMatchingCluster2->Draw("hist");
    validOutTrackMismatchingCluster2->Draw("same hist");

    auto canvasTracks = new TCanvas("canvasTracks", "Tracks", cx, cy);
    canvasTracks->cd();
    canvasTracks->Divide(2);
    canvasTracks->cd(1);
    gPad->SetLogy();
    validTracksChi2->Draw("hist");
    fakeTracksChi2->Draw("same hist");
    canvasTracks->cd(2);
    // gPad->SetLogy();
    validTracksPt->Draw("hist");
    fakeTracksPt->Draw("same hist");

    canvasL4->SaveAs("layer4dxdy.png");
    canvasL3->SaveAs("layer3dxdy.png");
    canvasL2->SaveAs("layer2dxdy.png");
    canvasL4_chi2->SaveAs("layer4_chi2.png");
    canvasL3_chi2->SaveAs("layer3_chi2.png");
    canvasL2_chi2->SaveAs("layer2_chi2.png");
    canvasSumChi2->SaveAs("sumChi2.png");
    canvasL4Matches_chi2->SaveAs("layer4_matchesChi2.png");
    canvasL3Matches_chi2->SaveAs("layer3_matchesChi2.png");
    canvasL2Matches_chi2->SaveAs("layer2_matchesChi2.png");
    canvasTracks->SaveAs("trackParams.png");
}