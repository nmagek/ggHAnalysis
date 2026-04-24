// run.C
// Usage inside ROOT:
//   root [0] .L ggHAnalysis.C+
//   root [1] .x run.C(1)   // signalΜ12
//   root [2] .x run.C(2)   // signalM15
//   root [3] .x run.C(3)   // signalM20
//   root [4] .x run.C(4)   // signalM25
//   root [5] .x run.C(5)   // signal30
//   root [6] .x run.C(6)   // QCDbin600-800
//   root [7] .x run.C(7)   // QCDbin800-1000
//   root [8] .x run.C(8)   // QCDbin1000-1200
//   root [9] .x run.C(9)   // QCDbin1200-1500
//   root [10] .x run.C(10)   // QCDbin1500-2000
//   root [11] .x run.C(11)  //QCDbin2000-inf
//   root [12] .x run.C(12)  //TTto4Q

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>


void run(int sample = 1)
{
    TString infile;
    double sigma_pb = 1.0;        // set per-sample
    double lumi_fb  = 108.96;     // your integrated luminosity

    // -------------------------------
    // Select input file + sigma
    // -------------------------------
    if (sample == 1) {
        infile   = "GluGluH-01J_HToAATo4B_M-12_2024.root";
        sigma_pb = 52.23;  
        std::cout << "[run] Selected SIGNAL sample\n";
    }
    else if (sample == 2){
      infile = "GluGluH-01J_HToAATo4B_M-15_2024.root";
      sigma_pb = 52.23;
      std::cout << "[run] Selected Signal sample\n";
    }
    else if (sample == 3){
      infile = "GluGluH-01J_HToAATo4B_M-20_2024.root";
      sigma_pb = 52.23;
      std::cout << "[run] Selected Signal sample\n";
    }
    else if (sample == 4){
      infile = "GluGluH-01J_HToAATo4B_M-25_2024.root";
      sigma_pb = 52.23;
      std::cout << "[run] Selected Signal sample\n";
    }
    else if (sample == 5){
      infile = "GluGluH-01J_HToAATo4B_M-30_2024.root";
      sigma_pb = 52.23;
      std::cout << "[run] Selected Signal sample\n";
    }
    else if (sample == 6) {
        infile   = "QCD-4Jets_Bin-HT-600to800_2024.root";
        sigma_pb = 13540;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    
    else if (sample == 7) {
        infile   = "QCD-4Jets_Bin-HT-800to1000_2024.root";
        sigma_pb = 3033;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else if (sample == 8) {
        infile   = "QCD-4Jets_Bin-HT-1000-1200_2024.root";
        sigma_pb = 883.7;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else if (sample == 9) {
        infile   = "QCD-4Jets_Bin-HT-1200-1500_2024.root";
        sigma_pb = 383.5;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else if (sample == 10) {
        infile   = "QCD-4Jets_Bin-HT-1500-2000_2024.root";
        sigma_pb = 125.2;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else if (sample == 11) {
        infile   = "QCD-4Jets_Bin-HT-2000toinf_2024.root";
        sigma_pb = 26.49;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else if (sample == 12) {
      infile = "TTto4Q_2024.root";
      sigma_pb = 419.69;
      std::cout << "[run] Selected TTto4Q backround sample";
    }
    else {
        std::cout << "[run] ERROR: unknown option " << sample << "\n";
        std::cout << "       Use run(1) for signal or run(2) for QCD\n";
        return;
    }

    std::cout << "[run] Input file: " << infile << "\n";
    std::cout << "[run] Using sigma [pb] = " << sigma_pb << ", lumi [fb^-1] = " << lumi_fb << "\n";

    // -------------------------------
    // Build output file name from input file name
    // -------------------------------
    TString base = gSystem->BaseName(infile);   // e.g. GluGluH-...root
    base.ReplaceAll(".root", "");               // strip extension

    TString category = "unknown";
    if (base.BeginsWith("GluGluH")) category = "signal";
    else if (base.BeginsWith("QCD")) category = "qcd";

    TString outfile = Form("ggHAnalysis_%s_%s.root", category.Data(), base.Data());
    std::cout << "[run] Output file: " << outfile << "\n";

    
    // -------------------------------
    // Open input ROOT file
    // -------------------------------
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) {
        std::cout << "[run] ERROR: cannot open file " << infile << "\n";
        return;
    }

    // -------------------------------
    // Get Events tree
    // -------------------------------
    TTree *tree = nullptr;
    f->GetObject("Events", tree);
    if (!tree) {
        std::cout << "[run] ERROR: TTree \"Events\" not found in " << infile << "\n";
        f->Close();
        return;
    }

    // -------------------------------
    // Run analysis (weighted with sigma)
    // -------------------------------
    ggHAnalysis ana(tree);
    ana.Loop(sigma_pb, lumi_fb, outfile.Data());

    std::cout << "[run] Analysis completed successfully.\n";

    f->Close();
}
