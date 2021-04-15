#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TFile.h>
#endif

void dumpParticle(const int evtId, const int trackId, const string path = "./")
{
    auto kinFile = TFile::Open((path + "o2sim_Kine.root").data());
    auto simtree = (TTree *)kinFile->Get("o2sim");
    auto tracksBranch = simtree->GetBranch("TrackRefs");
    for (auto iTRef{0}; iTref < tracksBranch->GetEntriesFast(); ++iTRef)
    {
        tracksBranch->GetEntry(iTref);
        auto structure = 
    }
}