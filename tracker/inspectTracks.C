#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include "DataFormatsITS/TrackITS.h"

#endif

using namespace std;

void inspectTracks(string trackFileName = "o2trac_its.root")
{
  using namespace o2::its;

  ROOT::EnableImplicitMT();                          // Enable ROOT's implicit multi-threading
  ROOT::RDataFrame d("o2sim", trackFileName.data()); // Interface to TTree and TChain
  TH1F *chiHist = new TH1F("h", "Tracks #chi^2;#chi^2;counts", 500, 0, 300);
  d.Foreach([&chiHist](vector<TrackITS> &vTracks) { for (TrackITS& t : vTracks) { chiHist->Fill(t.getChi2()); } }, {"ITSTrack"});

  auto canvasChi = new TCanvas("c", "chisquare");
  canvasChi->cd();
  canvasChi->SetLogy();

  chiHist->Draw();
}
