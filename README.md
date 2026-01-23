# Higgs Event Analysis & ML Classification (ROOT / C++ / TMVA)

End-to-end data analysis and machine learning pipeline for Higgs production studies, built with ROOT/C++ and TMVA.

This project implements a full workflow:
- Event processing and normalization using cross section × luminosity
- Physics object selection (jets, leptons, b-tagging, MET, mass window cuts)
- Feature engineering of high-level kinematic variables
- Creation of an ML-ready dataset (`RecoAfterCuts` TTree)
- Training and evaluation of a TMVA BDT classifier (signal vs background)

---

## Pipeline Overview

1. **Data processing (ggHAnalysis.C)**  
   - Reads ROOT ntuples  
   - Applies event selection and data-quality checks  
   - Computes derived variables (m_bb, MET, HT, angular variables, etc.)  
   - Writes structured dataset (`RecoAfterCuts` tree)

2. **Modeling (TMVA_ggHClassification.C)**  
   - Trains a BDT classifier using engineered features  
   - Uses per-event weights  
   - Produces evaluation outputs (ROC, classifier response, etc.)

3. **Visualization (DrawggHAnalysis.C)**  
   - Generates plots of key distributions before/after selection

---

## Technologies used
- C++
- ROOT Framework
- TMVA (BDT classification)
- Scientific data analysis
- Feature engineering
- Git/GitHub

---

## How to run

```cpp
root -l
.L ggHAnalysis.C+
.x run.C
