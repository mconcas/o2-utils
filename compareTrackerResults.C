#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TRatioPlot.h>
#endif

const float dymin{-5.f};
const float dymax{5.f};
const float dzmin{-5.f};
const float dzmax{5.f};
const float minchi{0};
const float maxchi{2e3};
const int nBins{300};
const int cx{2000};
const int cy{1500};
const float minpt{0.};
const float maxpt{5.};

void compareTrackerResults(const std::string inFileNormal = "nosmooth/dbg_ITSTrackerCPU.root", const std::string inFileSmoothed = "smooth/dbg_ITSTrackerCPU.root")
{
    auto noSmoothedFile = TFile::Open(inFileNormal.data());
    auto trackTree = (TTree *)noSmoothedFile->Get("TrackParams");
    auto validTracksChi2 = new TH1F("validTracksChi2", "Valid final tracks #chi^{2};#chi^{2}",
                                    nBins, minchi, maxchi);
    validTracksChi2->SetLineColor(kBlack);
    auto fakeTracksChi2 = new TH1F("fakeTracksChi2", "Fake final tracks #chi^{2};#chi^{2}",
                                   nBins, minchi, maxchi);
    fakeTracksChi2->SetLineColor(kBlack);

    auto validTracksPt = new TH1F("validTracksPt", "Valid final tracks #it{p}_{T};#it{p}_{T}",
                                  nBins, minpt, maxpt);
    validTracksPt->SetLineColor(kBlack);
    auto fakeTracksPt = new TH1F("fakeTracksPt", "Fake final tracks #it{p}_{T};#it{p}_{T}",
                                 nBins, minpt, maxpt);
    fakeTracksPt->SetLineColor(kBlack);

    trackTree->Draw("chi2>>validTracksChi2", "!fake", "goff");
    trackTree->Draw("chi2>>fakeTracksChi2", "fake", "goff");
    trackTree->Draw("pt>>validTracksPt", "!fake", "goff");
    trackTree->Draw("pt>>fakeTracksPt", "fake", "goff");

    // validTracksChi2->Scale(1. / (fakeTracksChi2->GetEntries() + fakeTracksChi2->GetEntries()));
    // fakeTracksChi2->Scale(1. / (fakeTracksChi2->GetEntries() + fakeTracksChi2->GetEntries()));
    // validTracksPt->Scale(1. / (fakeTracksPt->GetEntries() + validTracksPt->GetEntries()));
    // fakeTracksPt->Scale(1. / (fakeTracksPt->GetEntries() + validTracksPt->GetEntries()));

    auto smoothedFile = TFile::Open(inFileSmoothed.data());
    auto trackTreeSmoothed = (TTree *)smoothedFile->Get("TrackParams");
    auto validTracksChi2Smoothed = new TH1F("validTracksChi2Smoothed", "Valid final tracks #chi^{2};#chi^{2}",
                                            nBins, minchi, maxchi);
    validTracksChi2Smoothed->SetLineColor(kBlue);
    auto fakeTracksChi2Smoothed = new TH1F("fakeTracksChi2Smoothed", "Fake final tracks #chi^{2};#chi^{2}",
                                           nBins, minchi, maxchi);
    fakeTracksChi2Smoothed->SetLineColor(kRed);

    auto validTracksPtSmoothed = new TH1F("validTracksPtSmoothed", "Valid final tracks #it{p}_{T};#it{p}_{T}",
                                          nBins, minpt, maxpt);
    validTracksPtSmoothed->SetLineColor(kBlue);
    auto fakeTracksPtSmoothed = new TH1F("fakeTracksPtSmoothed", "Fake final tracks #it{p}_{T};#it{p}_{T}",
                                         nBins, minpt, maxpt);
    fakeTracksPtSmoothed->SetLineColor(kRed);

    trackTreeSmoothed->Draw("chi2>>validTracksChi2Smoothed", "!fake", "goff");
    trackTreeSmoothed->Draw("chi2>>fakeTracksChi2Smoothed", "fake", "goff");
    trackTreeSmoothed->Draw("pt>>validTracksPtSmoothed", "!fake", "goff");
    trackTreeSmoothed->Draw("pt>>fakeTracksPtSmoothed", "fake", "goff");

    // validTracksChi2Smoothed->Scale(1. / (fakeTracksChi2Smoothed->GetEntries() + fakeTracksChi2Smoothed->GetEntries()));
    // fakeTracksChi2Smoothed->Scale(1. / (fakeTracksChi2Smoothed->GetEntries() + fakeTracksChi2Smoothed->GetEntries()));
    // validTracksPtSmoothed->Scale(1. / (fakeTracksPtSmoothed->GetEntries() + validTracksPtSmoothed->GetEntries()));
    // fakeTracksPtSmoothed->Scale(1. / (fakeTracksPtSmoothed->GetEntries() + validTracksPtSmoothed->GetEntries()));

    auto canvasTracks = new TCanvas("canvasTracks", "Tracks", cx, cy);
    canvasTracks->cd();
    canvasTracks->Divide(2);
    auto ratioFake = new TRatioPlot(fakeTracksPtSmoothed, fakeTracksPt);
    canvasTracks->cd(1);
    ratioFake->Draw();
    ratioFake->GetLowerRefGraph()->SetMinimum(0.f);
    ratioFake->GetLowerRefGraph()->SetMaximum(10.f);
    canvasTracks->cd(2);
    auto ratioTrue = new TRatioPlot(validTracksPtSmoothed, validTracksPt);
    ratioTrue->Draw();
    ratioTrue->GetLowerRefGraph()->SetMinimum(0.f);
    ratioTrue->GetLowerRefGraph()->SetMaximum(10.f);
    canvasTracks->Update();

    canvasTracks->SaveAs("compareTracker.png");
}
