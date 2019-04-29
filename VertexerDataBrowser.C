#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphPolar.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLeaf.h"
#include "TNtuple.h"
#include "TPaveStats.h"

void VertexerDataBrowser(
    const std::string fileNameResults, /*, const std::string fileNameData*/
    const std::string fileNameCompareTo = "") {
  TFile vertexerData(fileNameResults.data());
  // Comparison block
  char compare{!fileNameCompareTo.empty()};
  TFile *vertexerComparisonData{NULL};

  if (compare) {
    vertexerComparisonData = TFile::Open(fileNameCompareTo.data());
  }

  auto canvasClusters = new TCanvas("Clusters", "Clusters data", 1000, 1000);
  auto canvasTracklets = new TCanvas("Tracklets", "Tracklets data", 1000, 1000);
  auto canvasLines = new TCanvas("Lines", "Lines data", 1000, 1000);
  auto canvasPolar = new TCanvas("Polar", "Polar", 1000, 1000);
  auto canvasDebugTanL = new TCanvas("dbTanL", "dbTanL", 1000, 1000);

  TNtuple *lines = (TNtuple *)vertexerData.Get("Tracklets");
  TNtuple *comb01 = (TNtuple *)vertexerData.Get("comb01");
  TNtuple *comb12 = (TNtuple *)vertexerData.Get("comb12");
  TNtuple *phi01 = (TNtuple *)vertexerData.Get("clus_phi01");
  TNtuple *phi12 = (TNtuple *)vertexerData.Get("clus_phi12");
  TNtuple *deltaTanLambda = (TNtuple *)vertexerData.Get("dtl");
  TNtuple *centroids = (TNtuple *)vertexerData.Get("centroids");
  TNtuple *linesData = (TNtuple *)vertexerData.Get("ld");

  TNtuple *linesRef{NULL};
  TNtuple *comb01Ref{NULL};
  TNtuple *comb12Ref{NULL};
  TNtuple *phi01Ref{NULL};
  TNtuple *phi12Ref{NULL};
  TNtuple *deltaTanLambdaRef{NULL};
  TNtuple *centroidsRef{NULL};
  TNtuple *linesDataRef{NULL};

  // Comparison data
  if (compare) {
    linesRef = (TNtuple *)vertexerComparisonData->Get("Tracklets");
    comb01Ref = (TNtuple *)vertexerComparisonData->Get("comb01");
    comb12Ref = (TNtuple *)vertexerComparisonData->Get("comb12");
    phi01Ref = (TNtuple *)vertexerComparisonData->Get("clus_phi01");
    phi12Ref = (TNtuple *)vertexerComparisonData->Get("clus_phi12");
    deltaTanLambdaRef = (TNtuple *)vertexerComparisonData->Get("dtl");
    centroidsRef = (TNtuple *)vertexerComparisonData->Get("centroids");
    linesDataRef = (TNtuple *)vertexerComparisonData->Get("ld");
  }

  TH1F *histComb01Phi =
      new TH1F("hComb01phi", "Comb01 Phi angle", 200, -3.15f, 3.15f);
  TH1F *histComb12Phi =
      new TH1F("hComb12phi", "Comb12 Phi angle", 200, -3.15f, 3.15f);
  TH1F *histClus0Phi =
      new TH1F("hCl0phi", "Clusters Layer 0 Phi", 200, 0.f, 6.3f);
  TH1F *histClus1Phi =
      new TH1F("hCl1phi", "Clusters Layer 1 Phi", 200, 0.f, 6.3f);
  TH1F *histClus2Phi =
      new TH1F("hCl2phi", "Clusters Layer 2 Phi", 200, 0.f, 6.3f);
  TH1F *histComb01TanLambda =
      new TH1F("hComb01tlambda", "Comb01 tanLambda", 100, -65.f, 65.f);
  TH1F *histComb12TanLambda =
      new TH1F("hComb12tlambda", "Comb12 tanLambda", 100, -65.f, 65.f);
  TH2F *histDelta01Phi1 = new TH2F("hDelta01Phi1", "DeltaPhi 01 vs Phi1", 100,
                                   0.f, 6.29f, 200, -1.f, 1.f);
  TH2F *histDelta12Phi1 = new TH2F("hDelta12Phi1", "DeltaPhi 12 vs Phi1", 100,
                                   0.f, 6.29f, 200, -1.f, 1.f);
  TH2F *histPhi0VsPhi1 =
      new TH2F("hPhi0Phi1", "phi0 vs phi1", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH2F *histPhi1VsPhi2 =
      new TH2F("hPhi1Phi2", "phi1 vs phi2", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH1F *histPhi01 = new TH1F("hPhi01", "Phi0 - Phi1", 200, -1.f, 1.f);
  TH1F *histPhi12 = new TH1F("hPhi12", "Phi1 - Phi2", 200, -1.f, 1.f);
  TH1F *histOrigX = new TH1F("hOrigX", "Origin X", 200, -3, 3);
  TH1F *histOrigY = new TH1F("hOrigY", "Origin Y", 200, -3, 3);
  TH1F *histOrigZ = new TH1F("hOrigZ", "Origin Z", 200, -25, 25);
  TH1F *histCosDir1 =
      new TH1F("hCosDir1", "Cosine director 1", 100, -1.2f, 1.2f);
  TH1F *histCosDir2 =
      new TH1F("hCosDir2", "Cosine director 2", 100, -1.2f, 1.2f);
  TH1F *histCosDir3 =
      new TH1F("hCosDir3", "Cosine director 3", 100, -1.2f, 1.2f);
  TH2F *histClusTanLambda = new TH2F("hClusTanLambda", "Clusters on 3 Layers",
                                     300, -16.33f, 16.33f, 300, 0.f, 5.f);

  histComb01Phi->SetMinimum(0);
  histComb12Phi->SetMinimum(0);
  histClus0Phi->SetMinimum(0);
  histClus1Phi->SetMinimum(0);
  histClus2Phi->SetMinimum(0);

  TH1F *histComb01Phi_ref{NULL};
  TH1F *histComb12Phi_ref{NULL};
  TH1F *histComb01Phi_diff{NULL};
  TH1F *histComb12Phi_diff{NULL};

  TH1F *histClus0Phi_ref{NULL};
  TH1F *histClus1Phi_ref{NULL};
  TH1F *histClus2Phi_ref{NULL};

  TH1F *histClus0Phi_diff{NULL};
  TH1F *histClus1Phi_diff{NULL};
  TH1F *histClus2Phi_diff{NULL};

  TH1F *histComb01TanLambda_ref{NULL};
  TH1F *histComb12TanLambda_ref{NULL};
  TH1F *histComb01TanLambda_diff{NULL};
  TH1F *histComb12TanLambda_diff{NULL};

  TH1F *histPhi01_ref{NULL};
  TH1F *histPhi12_ref{NULL};
  TH1F *histPhi01_diff{NULL};
  TH1F *histPhi12_diff{NULL};

  TH2F *histDelta01Phi1_ref{NULL};
  TH2F *histDelta12Phi1_ref{NULL};

  TH2F *histPhi0VsPhi1_ref{NULL};
  TH2F *histPhi1VsPhi2_ref{NULL};

  TH1F *histOrigX_ref{NULL};
  TH1F *histOrigY_ref{NULL};
  TH1F *histOrigZ_ref{NULL};
  TH1F *histCosDir1_ref{NULL};
  TH1F *histCosDir2_ref{NULL};
  TH1F *histCosDir3_ref{NULL};

  TH1F *histOrigX_diff{NULL};
  TH1F *histOrigY_diff{NULL};
  TH1F *histOrigZ_diff{NULL};
  TH1F *histCosDir1_diff{NULL};
  TH1F *histCosDir2_diff{NULL};
  TH1F *histCosDir3_diff{NULL};

  TH2F *histClusTanLambda_ref{NULL};

  // Comparison data
  if (compare) {
    histComb01Phi_ref =
        new TH1F("hComb01phiRef", "Comb01 Phi angle", 200, -3.15f, 3.15f);
    histComb12Phi_ref =
        new TH1F("hComb12phiRef", "Comb12 Phi angle", 200, -3.15f, 3.15f);
    histClus0Phi_ref =
        new TH1F("hCl0phiRef", "Clusters Layer 0 Phi", 200, 0.f, 6.3f);
    histClus1Phi_ref =
        new TH1F("hCl1phiRef", "Clusters Layer 1 Phi", 200, 0.f, 6.3f);
    histClus2Phi_ref =
        new TH1F("hCl2phiRef", "Clusters Layer 2 Phi", 200, 0.f, 6.3f);
    histComb01TanLambda_ref =
        new TH1F("hComb01tlambdaRef", "Comb01 tanLambda", 100, -65.f, 65.f);
    histComb12TanLambda_ref =
        new TH1F("hComb12tlambdaRef", "Comb12 tanLambda", 100, -65.f, 65.f);
    histDelta01Phi1_ref = new TH2F("hDelta01Phi1Ref", "DeltaPhi 01 vs Phi1",
                                   100, 0.f, 6.29f, 200, -1.f, 1.f);
    histDelta12Phi1_ref = new TH2F("hDelta12Phi1Ref", "DeltaPhi 12 vs Phi1",
                                   100, 0.f, 6.29f, 200, -1.f, 1.f);
    histPhi0VsPhi1_ref = new TH2F("hPhi0Phi1Ref", "phi0 vs phi1", 200, 0.f,
                                  6.29f, 200, 0.f, 6.29f);
    histPhi1VsPhi2_ref = new TH2F("hPhi1Phi2Ref", "phi1 vs phi2", 200, 0.f,
                                  6.29f, 200, 0.f, 6.29f);
    histPhi01_ref = new TH1F("hPhi01Ref", "Phi0 - Phi1", 200, -1.f, 1.f);
    histPhi12_ref = new TH1F("hPhi12Ref", "Phi1 - Phi2", 200, -1.f, 1.f);
    histOrigX_ref = new TH1F("hOrigXRef", "Origin X", 200, -3, 3);
    histOrigY_ref = new TH1F("hOrigYRef", "Origin Y", 200, -3, 3);
    histOrigZ_ref = new TH1F("hOrigZRef", "Origin Z", 200, -25, 25);
    histCosDir1_ref =
        new TH1F("hCosDir1Ref", "Cosine director 1", 100, -1.2f, 1.2f);
    histCosDir2_ref =
        new TH1F("hCosDir2Ref", "Cosine director 2", 100, -1.2f, 1.2f);
    histCosDir3_ref =
        new TH1F("hCosDir3Ref", "Cosine director 3", 100, -1.2f, 1.2f);
    histClusTanLambda_ref =
        new TH2F("hClusTanLambdaRef", "Clusters on 3 Layers", 300, -16.33f,
                 16.33f, 300, 0.f, 5.f);
  }

  canvasTracklets->Divide(2, 2);
  canvasClusters->Divide(2, 3);
  canvasLines->Divide(3, 2);
  canvasPolar->Divide(2, 2);

  comb01->Project("hComb01phi", "phi");
  comb12->Project("hComb12phi", "phi");
  comb01->Project("hComb01tlambda", "tanLambda");
  comb12->Project("hComb12tlambda", "tanLambda");

  phi01->Project("hDelta01Phi1", "phi0-phi1:phi1");
  phi12->Project("hDelta12Phi1", "phi1-phi2:phi1");
  phi01->Project("hPhi01", "phi0-phi1");
  phi12->Project("hPhi12", "phi1-phi2");
  phi01->Project("hPhi0Phi1", "phi0:phi1");
  phi12->Project("hPhi1Phi2", "phi1:phi2");

  phi01->Project("hCl0phi", "phi0");
  phi01->Project("hCl1phi", "phi1");
  phi12->Project("hCl2phi", "phi2");

  lines->Project("hOrigX", "oX");
  lines->Project("hOrigY", "oY");
  lines->Project("hOrigZ", "oZ");
  lines->Project("hCosDir1", "c1");
  lines->Project("hCosDir2", "c2");
  lines->Project("hCosDir3", "c3");

  // Comparison data
  if (compare) {
    comb01Ref->Project("hComb01phiRef", "phi");
    comb12Ref->Project("hComb12phiRef", "phi");
    comb01Ref->Project("hComb01tlambdaRef", "tanLambda");
    comb12Ref->Project("hComb12tlambdaRef", "tanLambda");

    phi01Ref->Project("hDelta01Phi1Ref", "phi0-phi1:phi1");
    phi12Ref->Project("hDelta12Phi1Ref", "phi1-phi2:phi1");
    phi01Ref->Project("hPhi01Ref", "phi0-phi1");
    phi12Ref->Project("hPhi12Ref", "phi1-phi2");
    phi01Ref->Project("hPhi0Phi1Ref", "phi0:phi1");
    phi12Ref->Project("hPhi1Phi2Ref", "phi1:phi2");

    phi01Ref->Project("hCl0phiRef", "phi0");
    phi01Ref->Project("hCl1phiRef", "phi1");
    phi12Ref->Project("hCl2phiRef", "phi2");

    linesRef->Project("hOrigXRef", "oX");
    linesRef->Project("hOrigYRef", "oY");
    linesRef->Project("hOrigZRef", "oZ");
    linesRef->Project("hCosDir1Ref", "c1");
    linesRef->Project("hCosDir2Ref", "c2");
    linesRef->Project("hCosDir3Ref", "c3");
  }

  for (auto i{0}; i < deltaTanLambda->GetEntries(); ++i) {
    deltaTanLambda->GetEntry(i);
    // std::cout<<deltaTanLambda->GetLeaf("deltatanlambda")->GetValue()<<std::endl;
    if (deltaTanLambda->GetLeaf("deltatanlambda")->GetValue() > 0.0f) {
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c0z")->GetValue(),
                              deltaTanLambda->GetLeaf("c0r")->GetValue());
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c1z")->GetValue(),
                              deltaTanLambda->GetLeaf("c1r")->GetValue());
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c2z")->GetValue(),
                              deltaTanLambda->GetLeaf("c2r")->GetValue());
      // break;
    }
  }

  // Comparison data
  if (compare) {
    for (auto i{0}; i < deltaTanLambdaRef->GetEntries(); ++i) {
      deltaTanLambdaRef->GetEntry(i);
      // std::cout<<deltaTanLambda->GetLeaf("deltatanlambda")->GetValue()<<std::endl;
      if (deltaTanLambdaRef->GetLeaf("deltatanlambda")->GetValue() > 0.0f) {
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c0z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c0r")->GetValue());
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c1z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c1r")->GetValue());
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c2z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c2r")->GetValue());
        // break;
      }
    }
  }

  const int npoints = 200;

  double r01[npoints];
  double phiArr01[npoints];

  double r12[npoints];
  double phiArr12[npoints];

  for (int ipt = 1; ipt < npoints + 1; ipt++) {
    r01[ipt - 1] = histComb01Phi->GetXaxis()->GetBinCenter(ipt);
    phiArr01[ipt - 1] = histComb01Phi->GetBinContent(ipt);
    r12[ipt - 1] = histComb12Phi->GetXaxis()->GetBinCenter(ipt);
    phiArr12[ipt - 1] = histComb12Phi->GetBinContent(ipt);
  }

  // Comparison data
  double r01Ref[npoints];
  double phiArr01Ref[npoints];

  double r12Ref[npoints];
  double phiArr12Ref[npoints];

  if (compare) {
    for (int ipt = 1; ipt < npoints + 1; ipt++) {
      r01Ref[ipt - 1] = histComb01Phi_ref->GetXaxis()->GetBinCenter(ipt);
      phiArr01Ref[ipt - 1] = histComb01Phi_ref->GetBinContent(ipt);
      r12Ref[ipt - 1] = histComb12Phi_ref->GetXaxis()->GetBinCenter(ipt);
      phiArr12Ref[ipt - 1] = histComb12Phi_ref->GetBinContent(ipt);
    }
  }

  TGraphPolar *grP01 = new TGraphPolar(npoints, r01, phiArr01);
  TGraphPolar *grP12 = new TGraphPolar(npoints, r12, phiArr12);
  grP01->SetLineWidth(2);
  grP01->SetLineColor(kBlue + 2);
  grP01->SetMinRadial(0);
  grP12->SetLineWidth(2);
  grP12->SetLineColor(kRed - 3);
  grP12->SetMinRadial(0);
  canvasPolar->cd(4);
  grP01->Draw("LEP");
  grP12->Draw("LEP");
  canvasPolar->Update();
  grP01->GetPolargram()->SetRadialLabelSize(0.025);
  grP12->GetPolargram()->SetRadialLabelSize(0.025);

  histClus0Phi->SetDirectory(0);
  histClus0Phi->SetLineColor(kOrange + 8);
  canvasPolar->cd(1);
  histClus0Phi->Draw();
  histClus1Phi->SetDirectory(0);
  histClus1Phi->SetLineColor(kOrange + 8);
  canvasPolar->cd(2);
  histClus1Phi->Draw();
  histClus2Phi->SetDirectory(0);
  histClus2Phi->SetLineColor(kOrange + 8);
  canvasPolar->cd(3);
  histClus2Phi->Draw();

  // Comparison data
  if (compare) {
    TGraphPolar *grP01Ref = new TGraphPolar(npoints, r01Ref, phiArr01Ref);
    TGraphPolar *grP12Ref = new TGraphPolar(npoints, r12Ref, phiArr12Ref);
    grP01Ref->SetLineWidth(2);
    grP01Ref->SetLineColor(kAzure + 6);
    grP01Ref->SetMinRadial(0);
    grP12Ref->SetLineWidth(2);
    grP12Ref->SetLineColor(kMagenta + 1);
    grP12Ref->SetMinRadial(0);
    canvasPolar->cd(4);
    grP01Ref->Draw("LEP");
    grP12Ref->Draw("LEP");
    canvasPolar->Update();
    grP01Ref->GetPolargram()->SetRadialLabelSize(0.025);
    grP12Ref->GetPolargram()->SetRadialLabelSize(0.025);

    histClus0Phi_ref->SetDirectory(0);
    canvasPolar->cd(1);
    histClus0Phi_ref->Draw("same");
    histClus1Phi_ref->SetDirectory(0);
    canvasPolar->cd(2);
    histClus1Phi_ref->Draw("same");
    histClus2Phi_ref->SetDirectory(0);
    canvasPolar->cd(3);
    histClus2Phi_ref->Draw("same");

    // Diff histos
    histClus0Phi_diff = (TH1F *)histClus0Phi_ref->Clone("hCl0phiDiff");
    histClus0Phi_diff->Add(histClus0Phi, -1);
    histClus0Phi_diff->SetDirectory(0);
    histClus0Phi_diff->SetFillStyle(3001);
    histClus0Phi_diff->SetFillColor(kBlue);
    canvasPolar->cd(1);
    // TPaveStats* ps0 = (TPaveStats*) histClus0Phi_diff->FindObject("stats");
    histClus0Phi_diff->Draw("sames");
    histClus1Phi_diff = (TH1F *)histClus1Phi_ref->Clone("hCl0phiDiff");
    histClus1Phi_diff->Add(histClus1Phi, -1);
    histClus1Phi_diff->SetDirectory(0);
    histClus1Phi_diff->SetFillStyle(3001);
    histClus1Phi_diff->SetFillColor(kBlue);
    canvasPolar->cd(2);
    // TPaveStats* ps1 = (TPaveStats*) histClus1Phi_diff->FindObject("stats");
    histClus1Phi_diff->Draw("sames");
    histClus2Phi_diff = (TH1F *)histClus2Phi_ref->Clone("hCl0phiDiff");
    histClus2Phi_diff->Add(histClus2Phi, -1);
    histClus2Phi_diff->SetDirectory(0);
    histClus2Phi_diff->SetFillStyle(3001);
    histClus2Phi_diff->SetFillColor(kBlue);
    canvasPolar->cd(3);
    // TPaveStats* ps2 = (TPaveStats*) histClus2Phi_diff->FindObject("stats");
    histClus2Phi_diff->Draw("sames");
  }

  // Lines
  histComb01Phi->SetDirectory(0);
  histComb01Phi->SetLineColor(kOrange + 8);
  canvasTracklets->cd(1);
  histComb01Phi->Draw();
  histComb12Phi->SetDirectory(0);
  histComb12Phi->SetLineColor(kOrange + 8);
  canvasTracklets->cd(2);
  histComb12Phi->Draw();
  histComb01TanLambda->SetDirectory(0);
  histComb01TanLambda->SetLineColor(kOrange + 8);
  canvasTracklets->cd(3)->SetLogy();
  histComb01TanLambda->Draw();
  histComb12TanLambda->SetDirectory(0);
  histComb12TanLambda->SetLineColor(kOrange + 8);
  canvasTracklets->cd(4)->SetLogy();
  histComb12TanLambda->Draw();

  if (compare) {
    histComb01Phi_ref->SetDirectory(0);
    histComb01Phi_ref->SetLineColor(kBlue + 2);
    canvasTracklets->cd(1);
    histComb01Phi_ref->Draw("sames");
    histComb12Phi_ref->SetDirectory(0);
    histComb12Phi_ref->SetLineColor(kBlue + 2);
    canvasTracklets->cd(2);
    histComb12Phi_ref->Draw("sames");
    histComb01TanLambda_ref->SetDirectory(0);
    histComb01TanLambda_ref->SetLineColor(kBlue + 2);
    canvasTracklets->cd(3)->SetLogy();
    histComb01TanLambda_ref->Draw("sames");
    histComb12TanLambda_ref->SetDirectory(0);
    histComb12TanLambda_ref->SetLineColor(kBlue + 2);
    canvasTracklets->cd(4)->SetLogy();
    histComb12TanLambda_ref->Draw("sames");

    histComb01Phi_diff = (TH1F *)histComb01Phi_ref->Clone("hComb01phiDiff");
    histComb01Phi_diff->Add(histComb01Phi, -1);
    histComb01Phi_diff->SetDirectory(0);
    histComb01Phi_diff->SetFillStyle(3001);
    histComb01Phi_diff->SetFillColor(kBlue);
    canvasTracklets->cd(1);
    histComb01Phi_diff->Draw("sames");
    histComb12Phi_diff = (TH1F *)histComb12Phi_ref->Clone("hComb12phiDiff");
    histComb12Phi_diff->Add(histComb12Phi, -1);
    histComb12Phi_diff->SetDirectory(0);
    histComb12Phi_diff->SetFillStyle(3001);
    histComb12Phi_diff->SetFillColor(kBlue);
    canvasTracklets->cd(2);
    histComb12Phi_diff->Draw("sames");
    histComb01TanLambda_diff =
        (TH1F *)histComb01TanLambda_ref->Clone("hComb01tlambdaDiff");
    histComb01TanLambda_diff->Add(histComb01TanLambda, -1);
    histComb01TanLambda_diff->SetDirectory(0);
    histComb01TanLambda_diff->SetFillStyle(3001);
    histComb01TanLambda_diff->SetFillColor(kBlue);
    canvasTracklets->cd(3);
    histComb01TanLambda_diff->Draw("sames");
    histComb12TanLambda_diff =
        (TH1F *)histComb12TanLambda_ref->Clone("hComb12tlambdaDiff");
    histComb12TanLambda_diff->Add(histComb12TanLambda, -1);
    histComb12TanLambda_diff->SetDirectory(0);
    histComb12TanLambda_diff->SetFillStyle(3001);
    histComb12TanLambda_diff->SetFillColor(kBlue);
    canvasTracklets->cd(4);
    histComb12TanLambda_diff->Draw("sames");
  }

  // Phi Clusters block
  histDelta01Phi1->SetDirectory(0);
  canvasClusters->cd(1);
  histDelta01Phi1->Draw("colz");
  histDelta12Phi1->SetDirectory(0);
  canvasClusters->cd(2);
  histDelta12Phi1->Draw("colz");
  histPhi0VsPhi1->SetDirectory(0);
  canvasClusters->cd(3);
  histPhi0VsPhi1->Draw("colz");
  histPhi1VsPhi2->SetDirectory(0);
  canvasClusters->cd(4);
  histPhi1VsPhi2->Draw("colz");
  histPhi01->SetDirectory(0);
  histPhi01->SetLineColor(kOrange + 8);
  histPhi01->SetMinimum(0.0001);
  canvasClusters->cd(5)->SetLogy();
  histPhi01->Draw();
  histPhi12->SetDirectory(0);
  histPhi12->SetLineColor(kOrange + 8);
  histPhi12->SetMinimum(0.0001);
  canvasClusters->cd(6)->SetLogy();
  histPhi12->Draw();

  // Comparison data
  histPhi01_ref->SetDirectory(0);
  canvasClusters->cd(5);
  histPhi01_ref->Draw("sames");
  histPhi12_ref->SetDirectory(0);
  canvasClusters->cd(6);
  histPhi12_ref->Draw("sames");

  histPhi01_diff = (TH1F *)histPhi01_ref->Clone("hPhi01Diff");
  histPhi01_diff->Add(histPhi01, -1);
  histPhi01_diff->SetDirectory(0);
  histPhi01_diff->SetFillStyle(3001);
  histPhi01_diff->SetFillColor(kBlue);
  canvasClusters->cd(5);
  histPhi01_diff->Draw("sames");

  histPhi12_diff = (TH1F *)histPhi12_ref->Clone("hPhi12Diff");
  histPhi12_diff->Add(histPhi12, -1);
  histPhi12_diff->SetDirectory(0);
  histPhi12_diff->SetFillStyle(3001);
  histPhi12_diff->SetFillColor(kBlue);
  canvasClusters->cd(6);
  histPhi12_diff->Draw("sames");

  // Lines block
  histOrigX->SetDirectory(0);
  canvasLines->cd(1)->SetLogy();
  histOrigX->SetLineColor(kOrange + 8);
  histOrigX->Draw();
  histOrigY->SetDirectory(0);
  canvasLines->cd(2)->SetLogy();
  histOrigY->SetLineColor(kOrange + 8);
  histOrigY->Draw();
  histOrigZ->SetDirectory(0);
  canvasLines->cd(3)->SetLogy();
  histOrigZ->SetLineColor(kOrange + 8);
  histOrigZ->Draw();
  histCosDir1->SetDirectory(0);
  canvasLines->cd(4)->SetLogy();
  histCosDir1->SetLineColor(kOrange + 8);
  histCosDir1->Draw();
  histCosDir2->SetDirectory(0);
  canvasLines->cd(5)->SetLogy();
  histCosDir2->SetLineColor(kOrange + 8);
  histCosDir2->Draw();
  histCosDir3->SetDirectory(0);
  canvasLines->cd(6)->SetLogy();
  histCosDir3->SetLineColor(kOrange + 8);
  histCosDir3->Draw();

  histOrigX_ref->SetDirectory(0);
  canvasLines->cd(1);
  histOrigX_ref->Draw("same");
  histOrigY_ref->SetDirectory(0);
  canvasLines->cd(2);
  histOrigY_ref->Draw("same");
  histOrigZ_ref->SetDirectory(0);
  canvasLines->cd(3);
  histOrigZ_ref->Draw("same");
  histCosDir1_ref->SetDirectory(0);
  canvasLines->cd(4);
  histCosDir1_ref->Draw("same");
  histCosDir2_ref->SetDirectory(0);
  canvasLines->cd(5);
  histCosDir2_ref->Draw("same");
  histCosDir3_ref->SetDirectory(0);
  canvasLines->cd(6);
  histCosDir3_ref->Draw("same");

  histOrigX_diff = (TH1F *)histOrigX_ref->Clone("hOrigXDiff");
  histOrigX_diff->Add(histOrigX, -1);
  histOrigX_diff->SetDirectory(0);
  histOrigX_diff->SetFillStyle(3001);
  histOrigX_diff->SetFillColor(kBlue);
  canvasLines->cd(1);
  histOrigX_diff->Draw("sames");
  histOrigY_diff = (TH1F *)histOrigY_ref->Clone("hOrigYDiff");
  histOrigY_diff->Add(histOrigY, -1);
  histOrigY_diff->SetDirectory(0);
  histOrigY_diff->SetFillStyle(3001);
  histOrigY_diff->SetFillColor(kBlue);
  canvasLines->cd(2);
  histOrigY_diff->Draw("sames");
  histOrigZ_diff = (TH1F *)histOrigZ_ref->Clone("hOrigZDiff");
  histOrigZ_diff->Add(histOrigZ, -1);
  histOrigZ_diff->SetDirectory(0);
  histOrigZ_diff->SetFillStyle(3001);
  histOrigZ_diff->SetFillColor(kBlue);
  canvasLines->cd(3);
  histOrigZ_diff->Draw("sames");
  histCosDir1_diff = (TH1F *)histCosDir1_ref->Clone("hCosDir1Diff");
  histCosDir1_diff->Add(histCosDir1, -1);
  histCosDir1_diff->SetDirectory(0);
  histCosDir1_diff->SetFillStyle(3001);
  histCosDir1_diff->SetFillColor(kBlue);
  canvasLines->cd(4);
  histCosDir1_diff->Draw("sames");
  histCosDir2_diff = (TH1F *)histCosDir2_ref->Clone("hCosDir2Diff");
  histCosDir2_diff->Add(histCosDir2, -1);
  histCosDir2_diff->SetDirectory(0);
  histCosDir2_diff->SetFillStyle(3001);
  histCosDir2_diff->SetFillColor(kBlue);
  canvasLines->cd(5);
  histCosDir2_diff->Draw("sames");
  histCosDir3_diff = (TH1F *)histCosDir3_ref->Clone("hCosDir3Diff");
  histCosDir3_diff->Add(histCosDir3, -1);
  histCosDir3_diff->SetDirectory(0);
  histCosDir3_diff->SetFillStyle(3001);
  histCosDir3_diff->SetFillColor(kBlue);
  canvasLines->cd(6);
  histCosDir3_diff->Draw("sames");

  // Debug clusters block
  histClusTanLambda->SetDirectory(0);
  canvasDebugTanL->cd();
  histClusTanLambda->Draw("colz");
}
#endif
