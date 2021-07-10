#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"
#include <iostream>


int main(int argc, char* argv[]) {
//  HM
//  std::vector<float> mTppBins = { 1.08, 1.32, 1.65, 4.5 };
//  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
 std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
  double norm1 = 0.;
  double norm2 = 6.;
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(norm1,norm2);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  // p-antiL
  DreamDist* pAp = DreamFile->GetPairDistributions(0, 1, "");
  DreamCF* CFpApDef = CATSinput->ObtainCFSystBBar(1, "pApDef", pAp, nullptr);//rebin = 1
  // p-L
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(1, "ppDef", pp, ApAp);//rebin = 1
  const int pairCountsDefault = CFpApDef->GetFemtoPairs(0, 0.2);
  const int pairCountsDefault_pp = CFppDef->GetFemtoPairs(0, 0.2);

    TH1F* DefRebin = CFpApDef->FindCorrelationFunction("hCk_FixShiftedpApDef_0");
    TH1F* DefRebinMeV = CFpApDef->FindCorrelationFunction(
        "hCk_FixShiftedpApDefMeV_0");
    TH1F* DefReweighted = CFpApDef->FindCorrelationFunction(
        "hCk_ReweightedpApDef_0");
    TH1F* DefReweightedMeV = CFpApDef->FindCorrelationFunction(
        "hCk_ReweightedpApDefMeV_0");

    TH1F* DefRebin_pp = CFppDef->FindCorrelationFunction("hCk_FixShiftedppDef_0");
    TH1F* DefRebinMeV_pp = CFppDef->FindCorrelationFunction(
        "hCk_FixShiftedppDefMeV_0");
    TH1F* DefReweighted_pp = CFppDef->FindCorrelationFunction(
        "hCk_ReweightedppDef_0");
    TH1F* DefReweightedMeV_pp = CFppDef->FindCorrelationFunction(
        "hCk_ReweightedppDefMeV_0");

    DreamFile->SetQuite();
    DreamKayTee* mTppDists;
    DreamKayTee* mTpApDists;

    DreamFile->ReadmTHistos(filename, prefix, "8");
    mTpApDists = DreamFile->GetmTPairDistributionsBBar(0, 1);
    mTppDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

    mTpApDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTpApDists->SetKayTeeBins(mTppBins);
    mTpApDists->SetNormalization(norm1, norm2);
    mTpApDists->SetRebin( { 1 });
    mTpApDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpApDists->ObtainTheCorrelationFunctionBBar(
        gSystem->pwd(), prefix, TString::Format("pApDef"));


    mTppDists->SetSEMEReweightingRatio(DefRebin_pp, DefReweighted_pp, DefRebinMeV_pp,
                                       DefReweightedMeV_pp);
    mTppDists->SetKayTeeBins(mTppBins);
    mTppDists->SetNormalization(norm1,norm2);
    mTppDists->SetRebin( { 1 });
    mTppDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTppDists->ObtainTheCorrelationFunction(
        gSystem->pwd(), prefix, TString::Format("ppDef"));

  CFpApDef->WriteOutput("CFpAp.root");
  CFppDef->WriteOutput("CFpp.root");

}
