#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCut.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"

using namespace TMVA;

void TMVA_ggHClassification()
{
  // -------------------------------------------------------
  // Configuration
  // -------------------------------------------------------
  TString dir = ""; 
  TString outName = dir + "TMVA_ggH.root";

  // Input files: add as many as you want
  std::vector<TString> sigFiles = {
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-12_2024.root"
  };

  std::vector<TString> bkgFiles = {
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-800to1000_2024.root"
    // dir + "ggHAnalysis_qcd_QCD_HT1000to1500_2024.root",
    // ...
  };

  const TString treeName = "RecoAfterCuts";

  // -------------------------------------------------------
  // Output
  // -------------------------------------------------------
  TFile *outputFile = TFile::Open(outName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "ERROR: cannot create " << outName << "\n";
    return;
  }

  TMVA::Tools::Instance();

  TMVA::Factory *factory = new TMVA::Factory(
    "MVAnalysis",
    outputFile,
    "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification"
  );

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset_ggH");

  // -------------------------------------------------------
  // Variables: MUST match your RecoAfterCuts branches
  // -------------------------------------------------------
  dataloader->AddVariable("m_2b",        'F');
  dataloader->AddVariable("pt_2b",       'F');
  dataloader->AddVariable("eta_2b",      'F');
  dataloader->AddVariable("met",         'F');
  dataloader->AddVariable("ht",          'F');
  dataloader->AddVariable("nbjets",      'F');

  dataloader->AddVariable("dphi_j1j2",   'F');
  dataloader->AddVariable("dR_j1j2",     'F');
  dataloader->AddVariable("abs_dm_j1j2", 'F');

  dataloader->AddVariable("dphi_met_bb", 'F');
  dataloader->AddVariable("dphi_met_j1", 'F');
  dataloader->AddVariable("min_dphi_met_jet", 'F');

  dataloader->AddVariable("pt_j1",  'F');
  dataloader->AddVariable("eta_j1", 'F');
  dataloader->AddVariable("pt_j2",  'F');
  dataloader->AddVariable("eta_j2", 'F');

  // -------------------------------------------------------
  // Use per-event weight branch (this is the correct way)
  // -------------------------------------------------------
  dataloader->SetSignalWeightExpression("wgt");
  dataloader->SetBackgroundWeightExpression("wgt");

  // -------------------------------------------------------
  // Load signal trees
  // -------------------------------------------------------
  std::vector<TFile*> keepOpen;
  keepOpen.reserve(sigFiles.size() + bkgFiles.size());

  for (const auto &fn : sigFiles) {
    TFile *f = TFile::Open(fn, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: cannot open signal file " << fn << "\n";
      return;
    }
    keepOpen.push_back(f);

    TTree *t = (TTree*)f->Get(treeName);
    if (!t) {
      std::cerr << "ERROR: cannot find tree " << treeName << " in " << fn << "\n";
      return;
    }

    dataloader->AddSignalTree(t, 1.0);
  }

  // -------------------------------------------------------
  // Load background trees (multiple bins, multiple files)
  // -------------------------------------------------------
  for (const auto &fn : bkgFiles) {
    TFile *f = TFile::Open(fn, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: cannot open background file " << fn << "\n";
      return;
    }
    keepOpen.push_back(f);

    TTree *t = (TTree*)f->Get(treeName);
    if (!t) {
      std::cerr << "ERROR: cannot find tree " << treeName << " in " << fn << "\n";
      return;
    }

    dataloader->AddBackgroundTree(t, 1.0);
  }

  // -------------------------------------------------------
  // Preselection: remove invalid sentinel values if any
  // (You set min_dphi_met_jet = -1 when invalid.)
  // -------------------------------------------------------
  TCut preselectionCut = "min_dphi_met_jet >= 0 && met>=0 && ht>=0";

  // Train/test split
  dataloader->PrepareTrainingAndTestTree(
    preselectionCut, preselectionCut,
    "SplitMode=Random:NormMode=NumEvents:!V"
    // You can cap if you want:
    // ":nTrain_Signal=50000:nTrain_Background=50000:nTest_Signal=20000:nTest_Background=20000"
  );

  // -------------------------------------------------------
  // Methods
  // -------------------------------------------------------
  factory->BookMethod(
    dataloader,
    TMVA::Types::kBDT,
    "BDT",
    "!H:!V:"
    "NTrees=800:"
    "MinNodeSize=2.5%:"
    "MaxDepth=3:"
    "BoostType=AdaBoost:"
    "AdaBoostBeta=0.5:"
    "UseBaggedBoost:"
    "BaggedSampleFraction=0.5:"
    "SeparationType=GiniIndex:"
    "nCuts=20"
  );


  // Run
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();
  delete factory;
  delete dataloader;

  std::cout << "TMVA finished. Output: " << outName << "\n";
  // gROOT->ProcessLine(Form(".x $ROOTSYS/tutorials/tmva/TMVAGui.C(\"%s\")", outName.Data()));

}
