#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <iomanip>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"
// DEBUG
#include "ITStracking/ClusterLines.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/Cluster.h"

static bool compare(const o2::itsmft::MC2ROFRecord &a, const o2::itsmft::MC2ROFRecord &b)
{
  return (a.maxROF < b.maxROF);
}

void browseMCHeader()
{
  // Get event header
  TChain mcHeaderTree("o2sim");
  TChain itsClustersMC2ROF("ITSClustersMC2ROF");
  mcHeaderTree.AddFile("o2sim.root");
  itsClustersMC2ROF.AddFile("o2clus_its.root");

  o2::dataformats::MCEventHeader *mcHeader = nullptr;
  std::vector<o2::MCTrack> *mcTracks = nullptr;
  std::vector<o2::itsmft::MC2ROFRecord> *MC2ROFRecordsAccumPtr = nullptr;

  mcHeaderTree.GetBranch("MCEventHeader.");
  mcHeaderTree.SetBranchAddress("MCEventHeader.", &mcHeader);
  mcHeaderTree.SetBranchAddress("MCTrack", &mcTracks);

  itsClustersMC2ROF.GetBranch("ITSClustersMC2ROF");
  itsClustersMC2ROF.SetBranchAddress("ITSClustersMC2ROF", &MC2ROFRecordsAccumPtr);
  itsClustersMC2ROF.GetEntry(0);

  TFile *outfile = new TFile("ROFInfo.root", "recreate");
  TTree outTree("ROFInfo", "ROFInfo");
  std::vector<int> *roframe2evtsToSave = new std::vector<int>;
  std::vector<float> *primVert2evtsToSaveX = new std::vector<float>;
  std::vector<float> *primVert2evtsToSaveY= new std::vector<float>;
  std::vector<float> *primVert2evtsToSaveZ = new std::vector<float>;

  outTree.Branch("ROF2Evts", &roframe2evtsToSave);
  outTree.Branch("PrimaryVerticesX", &primVert2evtsToSaveX);
  outTree.Branch("PrimaryVerticesY", &primVert2evtsToSaveY);
  outTree.Branch("PrimaryVerticesZ", &primVert2evtsToSaveZ);

  std::vector<std::vector<int>> roframe2evts;
  std::vector<std::array<float, 3>> primaryVertices; // n elements = n events
  std::vector<std::vector<float>> vertices2roframesX; // n elements = n roframes
  std::vector<std::vector<float>> vertices2roframesY; // n elements = n roframes
  std::vector<std::vector<float>> vertices2roframesZ; // n elements = n roframes

  auto MC2ROFRecordsAccum = *MC2ROFRecordsAccumPtr;
  const int nRoframes = 1 + (*std::max_element(MC2ROFRecordsAccum.begin(), MC2ROFRecordsAccum.end(), compare)).maxROF;
  const int mcEntries = mcHeaderTree.GetEntries();

  roframe2evts.resize(nRoframes);
  primaryVertices.reserve(mcEntries);   // n elements = n events
  vertices2roframesX.resize(nRoframes); // n elements = n roframes
  vertices2roframesY.resize(nRoframes); // n elements = n roframes
  vertices2roframesZ.resize(nRoframes); // n elements = n roframes

  for (auto i{0}; i < mcEntries; ++i)
  {
    mcHeaderTree.GetEntry(i);
    for (auto &t : *mcTracks)
    {
      if (t.getMotherTrackId() == -1)
      {
        std::cout<<"MC entry: "<<i<<std::endl;
        std::cout << "Nominal position for the vertex: " << t.GetStartVertexCoordinatesX() << " " << t.GetStartVertexCoordinatesY() << " " << t.GetStartVertexCoordinatesZ() << std::endl;

        primaryVertices.emplace_back(std::array<float,3>{static_cast<float>(t.GetStartVertexCoordinatesX()), static_cast<float>(t.GetStartVertexCoordinatesY()), static_cast<float>(t.GetStartVertexCoordinatesZ())});
        break;
      }
    }
  }

  for (size_t i{0}; i < MC2ROFRecordsAccum.size(); ++i) // size is num_events
  {
    for (auto j{MC2ROFRecordsAccum[i].minROF}; j <= MC2ROFRecordsAccum[i].maxROF; ++j)
    {
      roframe2evts[j].push_back(i);
      vertices2roframesX[j].push_back(primaryVertices[i][0]);
      vertices2roframesY[j].push_back(primaryVertices[i][1]);
      vertices2roframesZ[j].push_back(primaryVertices[i][2]);
    }
    // std::cout << "Event " << std::setw(3) << i << ":\n\tstarts at: " << std::setw(3) << MC2ROFRecordsAccum[i].minROF << " ends at: " << std::setw(3) << MC2ROFRecordsAccum[i].maxROF << std::endl;
  }

  // std::cout << "Dumping :\n";
  for (size_t i{0}; i < roframe2evts.size(); ++i)
  {
    std::cout << "ROF " << i << ">\t";
    // std::copy(roframe2evts[i].begin(), roframe2evts[i].end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    for (size_t j{0}; j < vertices2roframesX[i].size(); ++j)
    {
      std::cout << "\tevent: " << roframe2evts[i][j] << "\t vertices: " << vertices2roframesX[i][j] << " " << vertices2roframesY[i][j] << " " << vertices2roframesZ[i][j] << "\n";
    }
    roframe2evtsToSave->swap(roframe2evts[i]);
    primVert2evtsToSaveX->swap(vertices2roframesX[i]);
    primVert2evtsToSaveY->swap(vertices2roframesY[i]);
    primVert2evtsToSaveZ->swap(vertices2roframesZ[i]);
    outTree.Fill();
  }

  // Output data

  // std::vector<std::vector<int>> *roframe2evtsToSave = new std::vector<std::vector<int>>;
  // std::vector<std::array<float, 3>> *evtVerticesToSave = new std::vector<std::array<float, 3>>;

  // roframe2evtsToSave->swap(roframe2evts);
  // evtVerticesToSave->swap(evtVertices);

  outTree.Write();
  outfile->Close();
}

#endif
