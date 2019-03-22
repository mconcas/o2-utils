#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

void VertexerDataBrowser(
    const std::string fileName =
        "/home/alidock/alice/data/000vtx-alien/vertexer_gpu_data.root")
{
  TFile vertexerData(fileName.data());
  auto canvasTracklets = new TCanvas("Tracklets", "Tracklets data", 1000, 400);
  auto canvasClusters = new TCanvas("Clusters", "Clusters data", 1000, 400);

  TNtuple *comb01 = (TNtuple *)vertexerData.Get("comb01");
  TNtuple *comb12 = (TNtuple *)vertexerData.Get("comb12");
  TNtuple *phi01 = (TNtuple *)vertexerData.Get("clus_phi01");
  TNtuple *phi12 = (TNtuple *)vertexerData.Get("clus_phi12");

  TH1F *histComb01Phi =
      new TH1F("hComb01phi", "Comb01 Phi angle", 100, -3.15f, 3.15f);
  TH1F *histComb12Phi =
      new TH1F("hComb12phi", "Comb12 Phi angle", 100, -3.15f, 3.15f);
  TH1F *histComb01TanLambda =
      new TH1F("hComb01tlambda", "Comb01 tanLambda", 100, -100.f, 100.f);
  TH1F *histComb12TanLambda =
      new TH1F("hComb12tlambda", "Comb12 tanLambda", 100, -100.f, 100.f);

  TH2F *histDelta01Phi1 = new TH2F("hDelta01Phi1", "DeltaPhi 01 vs Phi1", 100,
                                   0.f, 6.29f, 200, -0.5f, 0.5f);
  TH2F *histDelta12Phi1 = new TH2F("hDelta12Phi1", "DeltaPhi 12 vs Phi1", 100,
                                   0.f, 6.29f, 200, -0.5f, 0.5f);
  TH2F *histPhi0VsPhi1 =
      new TH2F("hPhi0Phi1", "phi0 vs phi1", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH2F *histPhi1VsPhi2 =
      new TH2F("hPhi1Phi2", "phi1 vs phi2", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH1F *histPhi01 = new TH1F("hPhi01", "Phi0 - Phi1", 200, -0.5f, 0.5f);
  TH1F *histPhi12 = new TH1F("hPhi12", "Phi1 - Phi2", 200, -0.5f, 0.5f);

  canvasTracklets->Divide(2, 2);
  canvasClusters->Divide(2, 3);

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
  // canvasTracklets->Draw();
  // comb01->Project("hComb01phi", "phi");
  // histComb01Phi->SetDirectory(0);
  // histComb01Phi->Draw();
  // comb01->Draw("phi:tanLambda");
}
#endif
