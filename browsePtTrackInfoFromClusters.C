// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include <gsl/gsl>

#include "DataFormatsITSMFT/Cluster.h"
#include "ITStracking/Cluster.h"

#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/ROframe.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"

void browsePtTrackInfoFromClusters(
    const std::string path = "./",
    const std::string inputClustersITS = "o2clus_its.root",
    const std::string simfilename = "o2sim.root",
    const std::string paramfilename = "O2geometry.root",
    const std::string outputfile = "label2Track.root") {
  // Setup Runtime DB
  TFile paramFile((path + paramfilename).data());
  paramFile.Get("FAIRGeom");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(
      o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                          o2::TransformType::L2G)); // request cached transforms

  // ROframe data
  TChain itsClustersROF("ITSClustersROF");
  TChain itsClusters("o2sim");
  itsClusters.AddFile((path + inputClustersITS).data());
  itsClustersROF.AddFile((path + inputClustersITS).data());
  std::vector<o2::itsmft::Cluster> *clusters = nullptr;
  std::vector<o2::itsmft::ROFRecord> *rofs = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels = nullptr;
  itsClusters.SetBranchAddress("ITSCluster", &clusters);
  itsClustersROF.SetBranchAddress("ITSClustersROF", &rofs);
  itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);
  itsClustersROF.GetEntry(0);
  itsClusters.GetEntry(0);

  // Simulation data
  TChain mcHeaderTree("o2sim");
  mcHeaderTree.AddFile((path + simfilename).data());
  std::vector<o2::MCTrack> *tracks = nullptr;
  mcHeaderTree.SetBranchAddress("MCTrack", &tracks);

  // Output data
  TFile *outfile = new TFile(outputfile.data(), "recreate");
  TTree outTree("Labels2Tracks", "labels2tracks");
  std::vector<o2::MCCompLabel> *labelsToSave = new std::vector<o2::MCCompLabel>;
  std::vector<o2::MCTrack> *tracksToSave = new std::vector<o2::MCTrack>;
  outTree.Branch("MCLabels", &labelsToSave);
  outTree.Branch("Tracks", &tracksToSave);

  o2::its::GeometryTGeo *geom = o2::its::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));

  int ROCounter{0};

  for (o2::itsmft::ROFRecord &rofRec : *rofs) {
    itsClustersROF.GetEntry(rofRec.getROFEntry().getEvent());
    int first = static_cast<int>(rofRec.getROFEntry().getIndex());
    int number = static_cast<int>(rofRec.getNROFEntries());
    std::vector<o2::MCCompLabel> roLabels;
    std::vector<o2::MCTrack> roTracks;
    std::cout << "ROframe: " << ROCounter++ << std::endl;
    for (int iCluster{first}; iCluster < first + number; ++iCluster) {
      auto &c = (*clusters)[iCluster];
      if (0 == geom->getLayer(c.getSensorID())) {
        auto &lbl = *labels->getLabels(first + iCluster).begin();
        if (lbl.isValid() && !lbl.isFake() && lbl.getSourceID() == 0) {
          if (mcHeaderTree.GetReadEntry() != lbl.getEventID()) {
            // to avoid repetitive reading of the same event
            mcHeaderTree.GetEntry(lbl.getEventID());
          }
          if (!(lbl.getTrackID() == 2147483647)) {
            const auto mcTr = (*tracks)[lbl.getTrackID()];
            std::cout << "\tid: " << lbl.getTrackID()
                      << " PDG: " << mcTr.GetPdgCode()
                      << " pT: " << mcTr.GetPt() << "\n";
            roLabels.push_back(lbl);
            roTracks.push_back(mcTr);

          } else {
            std::cout << "\t -> " << lbl.getTrackID()
                      << " evt ID: " << lbl.getEventID() << std::endl;
          }
        }
      }
    }
    labelsToSave->swap(roLabels);
    tracksToSave->swap(roTracks);
    outTree.Fill();
    if (ROCounter > 10)
      break;
  }
  outTree.Write();
  outfile->Close();
  
}
#endif