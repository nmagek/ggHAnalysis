#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>
#include "TChain.h"
#include <cmath>     // για floor
#include <cstdio>    // για Form (συνήθως το έχεις ήδη μέσω TString, αλλά βάλε το για σιγουριά)
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
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-12_2024.root",
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-15_2024.root",
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-20_2024.root",
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-25_2024.root",
    dir + "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-30_2024.root"
  };

  std::vector<TString> bkgFiles = {
    // dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-600to800_2024.root",
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-800to1000_2024.root",
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1000-1200_2024.root",
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1200-1500_2024.root",
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1500-2000_2024.root",
    dir + "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-2000toinf_2024.root"
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
  //dataloader->AddVariable("eta_2b",      'F');
  dataloader->AddVariable("met",         'F');
  dataloader->AddVariable("ht",          'F');
  //dataloader->AddVariable("nbjets",      'F');

  dataloader->AddVariable("dphi_j1j2",   'F');
  dataloader->AddVariable("dR_j1j2",     'F');
  dataloader->AddVariable("abs_dm_j1j2", 'F');

  dataloader->AddVariable("dphi_met_bb", 'F');
  dataloader->AddVariable("dphi_met_j1", 'F');
  dataloader->AddVariable("min_dphi_met_jet", 'F');

  dataloader->AddVariable("pt_j1",  'F');
  //dataloader->AddVariable("eta_j1", 'F');
  dataloader->AddVariable("pt_j2",  'F');
  //dataloader->AddVariable("eta_j2", 'F');

  // -------------------------------------------------------
  // Use per-event weight branch 
  // -------------------------------------------------------
  //dataloader->SetSignalWeightExpression("wgt");
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

    // --- count events AFTER preselection to set 60/40 split ---
    TChain sigCh(treeName);
    for (const auto &fn : sigFiles) sigCh.Add(fn);

    TChain bkgCh(treeName);
    for (const auto &fn : bkgFiles) bkgCh.Add(fn);

    // same preselection as you use in TMVA
    TCut preselectionCut = "min_dphi_met_jet >= 0 && met>=0 && ht>=0";

    Long64_t nSigTot = sigCh.GetEntries(preselectionCut);
    Long64_t nBkgTot = bkgCh.GetEntries(preselectionCut);

    // 60/40 split
    Long64_t nTrainSig = (Long64_t) std::floor(0.7 * (double)nSigTot);
    Long64_t nTestSig  = nSigTot - nTrainSig;

    Long64_t nTrainBkg = (Long64_t) std::floor(0.7 * (double)nBkgTot);
    Long64_t nTestBkg  = nBkgTot - nTrainBkg;

    std::cout << "[Split] Sig total=" << nSigTot
	      << " train=" << nTrainSig << " test=" << nTestSig << "\n";
    std::cout << "[Split] Bkg total=" << nBkgTot
	      << " train=" << nTrainBkg << " test=" << nTestBkg << "\n";

    // Now prepare with explicit numbers
    dataloader->PrepareTrainingAndTestTree(
					   preselectionCut, preselectionCut,
					   Form("nTrain_Signal=%lld:nTest_Signal=%lld:"
						"nTrain_Background=%lld:nTest_Background=%lld:"
						"SplitMode=Random:NormMode=NumEvents:!V",
						nTrainSig, nTestSig, nTrainBkg, nTestBkg)
					   );
    
    dataloader->AddBackgroundTree(t, 1.0);
  }

  // -------------------------------------------------------
  // Methods
  // -------------------------------------------------------
  factory->BookMethod(dataloader, TMVA::Types::kBDT,"BDT","!H:!V:NTrees=150:BoostType=Grad:Shrinkage=0.1:MaxDepth=2");


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
