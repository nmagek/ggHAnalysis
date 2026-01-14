// run.C
// Usage inside ROOT:
//   root [0] .L ggHAnalysis.C+
//   root [1] .x run.C(1)   // signal
//   root [2] .x run.C(2)   // QCD

#include <TFile.h>
#include <TTree.h>
#include <iostream>

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
        sigma_pb = 48.58;  
        std::cout << "[run] Selected SIGNAL sample\n";
    }
    else if (sample == 2) {
        infile   = "QCD-4Jets_Bin-HT-800to1000_2024.root";
        sigma_pb = 3033;  
        std::cout << "[run] Selected QCD BACKGROUND sample\n";
    }
    else {
        std::cout << "[run] ERROR: unknown option " << sample << "\n";
        std::cout << "       Use run(1) for signal or run(2) for QCD\n";
        return;
    }

    std::cout << "[run] Input file: " << infile << "\n";
    std::cout << "[run] Using sigma [pb] = " << sigma_pb << ", lumi [fb^-1] = " << lumi_fb << "\n";

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
    ana.Loop(sigma_pb, lumi_fb);

    std::cout << "[run] Analysis completed successfully.\n";

    f->Close();
}
