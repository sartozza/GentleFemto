// #include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"
#include <iostream>

int main(int argc, char* argv[]) {
//  HM
//  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
//  std::vector<float> mTppBins = { 1.08, 1.32, 1.65, 4.5 };
 std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  double norm1 = 0.;
  double norm2 = 6.;

  CATSinput->SetNormalization(norm1, norm2);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  // p-antiL
  DreamDist* LALCommon = DreamFile->GetPairDistributionsCommon(0, 1, "");
  DreamDist* LALNonCommon = DreamFile->GetPairDistributionsNonCommon(0, 1, "");

  DreamCF* CFLALDefCommon = CATSinput->ObtainCFSystBBar(1, "pApDefCommon", LALCommon, nullptr);//rebin = 5
  DreamCF* CFLALDefNonCommon = CATSinput->ObtainCFSystBBar(1, "pApDefNonCommon", LALNonCommon, nullptr);//rebin = 5

  const int pairCountsDefaultCommon = CFLALDefCommon->GetFemtoPairs(0, 0.2);
  const int pairCountsDefaultNonCommon = CFLALDefNonCommon->GetFemtoPairs(0, 0.2);

    TH1F* DefRebinCommon = CFLALDefCommon->FindCorrelationFunction("hCk_FixShiftedpApDefCommon_0");
    TH1F* DefRebinMeVCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_FixShiftedpApDefCommonMeV_0");
    TH1F* DefReweightedCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefCommon_0");
    TH1F* DefReweightedMeVCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefCommonMeV_0");

    TH1F* DefRebinNonCommon = CFLALDefNonCommon->FindCorrelationFunction("hCk_FixShiftedpApDefNonCommon_0");
    TH1F* DefRebinMeVNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_FixShiftedpApDefNonCommonMeV_0");
    TH1F* DefReweightedNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefNonCommon_0");
    TH1F* DefReweightedMeVNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefNonCommonMeV_0");

    DreamFile->SetQuite();
    DreamKayTee* mTLALDistsCommon;
    DreamKayTee* mTLALDistsNonCommon;


    DreamFile->ReadmTHistosAncestors(filename, prefix, "8");
    mTLALDistsCommon = DreamFile->GetmTPairDistributionsCommon(0, 1);
    mTLALDistsNonCommon = DreamFile->GetmTPairDistributionsNonCommon(0, 1);

    mTLALDistsCommon->SetSEMEReweightingRatio(DefRebinCommon, DefReweightedCommon, DefRebinMeVCommon,
                                       DefReweightedMeVCommon);
    mTLALDistsCommon->SetKayTeeBins(mTppBins);
    mTLALDistsCommon->SetNormalization(norm1, norm2);
    mTLALDistsCommon->SetRebin( { 1 });
    mTLALDistsCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLALDistsCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("pApDefCommon"), TString::Format("Common"));

    mTLALDistsNonCommon->SetSEMEReweightingRatio(DefRebinNonCommon, DefReweightedNonCommon, DefRebinMeVNonCommon,
                                       DefReweightedMeVNonCommon);
    mTLALDistsNonCommon->SetKayTeeBins(mTppBins);
    mTLALDistsNonCommon->SetNormalization(norm1, norm2);
    mTLALDistsNonCommon->SetRebin( { 1 });
    mTLALDistsNonCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLALDistsNonCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("pApDefNonCommon"), TString::Format("NonCommon"));

  CFLALDefCommon->WriteOutput("CFpApCommon.root");
  CFLALDefNonCommon->WriteOutput("CFpApNonCommon.root");

}
