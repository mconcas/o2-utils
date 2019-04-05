#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <string>

#include "TCanvas.h"
#include "TGraphPolar.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

void VertexerDataBrowser(const std::string fileNameResults /*, const std::string fileNameData*/)
{
    TFile vertexerData(fileNameResults.data());
    // TFile simData(fileNameData.data());

    auto canvasClusters = new TCanvas("Clusters", "Clusters data", 1000, 1000);
    auto canvasTracklets = new TCanvas("Tracklets", "Tracklets data", 1000, 1000);
    auto canvasLines = new TCanvas("Lines", "Lines data", 1000, 1000);
    auto canvasPolar = new TCanvas("Polar", "Polar", 1000, 1000);

    TNtuple *lines = (TNtuple *)vertexerData.Get("Tracklets");

    TNtuple *comb01 = (TNtuple *)vertexerData.Get("comb01");
    TNtuple *comb12 = (TNtuple *)vertexerData.Get("comb12");

    TNtuple *phi01 = (TNtuple *)vertexerData.Get("clus_phi01");
    TNtuple *phi12 = (TNtuple *)vertexerData.Get("clus_phi12");

    TH1F *histComb01Phi = new TH1F("hComb01phi", "Comb01 Phi angle", 200, -3.15f, 3.15f);
    histComb01Phi->SetMinimum(0);
    TH1F *histComb12Phi = new TH1F("hComb12phi", "Comb12 Phi angle", 200, -3.15f, 3.15f);
    histComb12Phi->SetMinimum(0);

    TH1F *histClus0Phi = new TH1F("hCl0phi", "Clusters Layer 0 Phi", 200, 0.f, 6.3f);
    histClus0Phi->SetMinimum(0);
    TH1F *histClus1Phi = new TH1F("hCl1phi", "Clusters Layer 1 Phi", 200, 0.f, 6.3f);
    histClus1Phi->SetMinimum(0);
    TH1F *histClus2Phi = new TH1F("hCl2phi", "Clusters Layer 2 Phi", 200, 0.f, 6.3f);
    histClus2Phi->SetMinimum(0);

    TH1F *histComb01TanLambda = new TH1F("hComb01tlambda", "Comb01 tanLambda", 100, -100.f, 100.f);
    TH1F *histComb12TanLambda = new TH1F("hComb12tlambda", "Comb12 tanLambda", 100, -100.f, 100.f);

    TH2F *histDelta01Phi1 = new TH2F("hDelta01Phi1", "DeltaPhi 01 vs Phi1", 100, 0.f, 6.29f, 200, -1.f, 1.f);
    TH2F *histDelta12Phi1 = new TH2F("hDelta12Phi1", "DeltaPhi 12 vs Phi1", 100, 0.f, 6.29f, 200, -1.f, 1.f);

    TH2F *histPhi0VsPhi1 = new TH2F("hPhi0Phi1", "phi0 vs phi1", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
    TH2F *histPhi1VsPhi2 = new TH2F("hPhi1Phi2", "phi1 vs phi2", 200, 0.f, 6.29f, 200, 0.f, 6.29f);

    TH1F *histPhi01 = new TH1F("hPhi01", "Phi0 - Phi1", 200, -1.f, 1.f);
    TH1F *histPhi12 = new TH1F("hPhi12", "Phi1 - Phi2", 200, -1.f, 1.f);

    TH1F *histOrigX = new TH1F("hOrigX", "Origin X", 200, -3, 3);
    TH1F *histOrigY = new TH1F("hOrigY", "Origin Y", 200, -3, 3);
    TH1F *histOrigZ = new TH1F("hOrigZ", "Origin Z", 200, -25, 25);

    TH1F *histCosDir1 = new TH1F("hCosDir1", "Cosine director 1", 100, -1.2f, 1.2f);
    TH1F *histCosDir2 = new TH1F("hCosDir2", "Cosine director 2", 100, -1.2f, 1.2f);
    TH1F *histCosDir3 = new TH1F("hCosDir3", "Cosine director 3", 100, -1.2f, 1.2f);

    canvasTracklets->Divide(2, 2);
    canvasClusters->Divide(2, 5);
    canvasLines->Divide(3, 2);

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

    const int npoints = 200;

    double r01[npoints];
    double phiArr01[npoints];

    double r12[npoints];
    double phiArr12[npoints];

    for (int ipt = 1; ipt < npoints + 1; ipt++)
    {
        r01[ipt - 1] = histComb01Phi->GetXaxis()->GetBinCenter(ipt);
        phiArr01[ipt - 1] = histComb01Phi->GetBinContent(ipt);
        r12[ipt - 1] = histComb12Phi->GetXaxis()->GetBinCenter(ipt);
        phiArr12[ipt - 1] = histComb12Phi->GetBinContent(ipt);
    }

    TGraphPolar *grP01 = new TGraphPolar(npoints, r01, phiArr01);
    TGraphPolar *grP12 = new TGraphPolar(npoints, r12, phiArr12);
    grP01->SetLineWidth(2);
    grP01->SetLineColor(3);
    grP01->SetMinRadial(0);
    grP12->SetLineWidth(2);
    grP12->SetLineColor(4);
    grP12->SetMinRadial(0);
    canvasPolar->cd();
    grP01->Draw("LEP");
    grP12->Draw("LEP");
    canvasPolar->Update();
    grP01->GetPolargram()->SetRadialLabelSize(0.025);

    histComb01Phi->SetDirectory(0);
    canvasTracklets->cd(1);
    histComb01Phi->Draw();
    histComb12Phi->SetDirectory(0);
    canvasTracklets->cd(2);
    histComb12Phi->Draw();
    histComb01TanLambda->SetDirectory(0);
    canvasTracklets->cd(3);
    histComb01TanLambda->Draw();
    histComb12TanLambda->SetDirectory(0);
    canvasTracklets->cd(4);
    histComb12TanLambda->Draw();

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
    canvasClusters->cd(5);
    histPhi01->Draw();
    histPhi12->SetDirectory(0);
    canvasClusters->cd(6);
    histPhi12->Draw();
    histClus0Phi->SetDirectory(0);
    canvasClusters->cd(7);
    histClus0Phi->Draw();
    histClus1Phi->SetDirectory(0);
    canvasClusters->cd(8);
    histClus1Phi->Draw();
    histClus2Phi->SetDirectory(0);
    canvasClusters->cd(9);
    histClus2Phi->Draw();

    // Lines block
    histOrigX->SetDirectory(0);
    canvasLines->cd(1);
    histOrigX->Draw();
    histOrigY->SetDirectory(0);
    canvasLines->cd(2);
    histOrigY->Draw();
    histOrigZ->SetDirectory(0);
    canvasLines->cd(3);
    histOrigZ->Draw();
    histCosDir1->SetDirectory(0);
    canvasLines->cd(4);
    histCosDir1->Draw();
    histCosDir2->SetDirectory(0);
    canvasLines->cd(5);
    histCosDir2->Draw();
    histCosDir3->SetDirectory(0);
    canvasLines->cd(6);
    histCosDir3->Draw();

    // canvasTracklets->Draw();
    // comb01->Project("hComb01phi", "phi");
    // histComb01Phi->SetDirectory(0);
    // histComb01Phi->Draw();
    // comb01->Draw("phi:tanLambda");
}
#endif
