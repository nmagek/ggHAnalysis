# ggHAnalysis
# Higgs Event Analysis Pipeline (ROOT/C++)

End-to-end event-analysis pipeline for Higgs production studies using ROOT.
The analysis performs dataset normalization using cross section and integrated luminosity, applies reconstruction-level object selections (leptons/jets), jet–lepton cleaning, b-tag requirements, MET and Higgs mass-window cuts, and writes both histograms and a structured output TTree for downstream studies.

## Features
- Event weighting: `w = (σ [pb] × L [fb⁻¹] × 10³) / N_events`
- GEN-level truth matching (PDG traversal + mother indices) for signal samples
- RECO-level selection:
  - loose lepton veto
  - jet selection + jet–lepton cleaning (ΔR < 0.4)
  - b-tag jet requirements
  - MET cut
  - Higgs candidate from 2 b-jets + mass window
- Cutflow report (efficiencies + weighted yields)
- Outputs:
  - ROOT histograms (kinematics before/after cuts)
  - `RecoAfterCuts` TTree with selected variables

## Requirements
- ROOT (tested with ROOT 6.x)
- C++ compiler available in ROOT environment

## ML Classification (TMVA)
This repository also includes a TMVA-based classifier to separate signal/background using high-level kinematic variables.


## How to run
From a ROOT session:

```cpp
root -l
.L ggHAnalysis.C+
.x run.C(1)   // signal
.x run.C(2)   // QCD background


### Run TMVA
root -l
.L TMVA_ggHClassification.C
TMVA_ggHClassification.C()
TMVA::TMVAGui("TMVA_ggH.root");
