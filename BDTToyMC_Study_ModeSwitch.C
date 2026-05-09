// BDTToyMC_Study_ModeSwitch.C
//
// Structured Bayesian toy-MC limit code.
// - Posterior plots have x-axis in Nexp_sig, including negative values.
// - Posterior is computed only in the physical region Nexp_sig >= 0.
// - Negative Nexp_sig region is shown with zero posterior.
//
// Usage:
//   root -l
//   .L BDTToyMC_Study_ModeSwitch.C+
//   BDTToyMC_Study_ModeSwitch();

#ifndef BDTTOYMC_STUDY_MODESWITCH_C
#define BDTTOYMC_STUDY_MODESWITCH_C

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cmath>

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooBinning.h"
#include "RooFit.h"

using namespace RooFit;

// ============================================================
// Options
// ============================================================
struct AnalysisOptions {
  std::vector<int> masses = {12, 15, 20, 25, 30};

  const char* histName = "h_BDT";

  // Use the same variable BDT binning style as the second code.
  int nbins = 5;
  double var_bins[6] = {-1.0, -0.70, -0.40, -0.08, 0.52, 1.0};

  // Toy and scan settings.
  int Ntoys_limits = 2001;      // use 1001 for final production
  int Ntoys_pulls  = 1000;

  int nNexpSigScan = 1000;
  int nNexpBkgScan = 400;

  double CL = 0.95;

  // Posterior scan range in Nexp_sig.
  // Negative values are allowed, but points that make any expected bin
  // non-positive are skipped before calling RooFit.
  double NexpSigMinPhysicalFactor = -1.0;
  double NexpSigMaxPhysicalFactor = 1.5;

  // Plot range in Nexp_sig.
  double NexpSigPlotMinFactor = -1.0;
  double NexpSigPlotMaxFactor = 1.5;

  double NexpBkgMinFactor = 0.5;
  double NexpBkgMaxFactor = 1.5;

  // Background Gaussian constraint width, currently constructed but not used
  // in the default posterior/fits unless you switch model_0/model_1 to total_model_0/1.
  double rel_bkg = 0.15;

  int seedBase = 12345;

  bool makeOneFitPlots = true;
  bool runPullStudy = false;
};

AnalysisOptions opt;

// ============================================================
// Input files from the first code
// ============================================================
struct Sample {
  std::string label;
  std::string fileName;
};

std::vector<Sample> backgrounds = {
  {"QCD 600-800",   "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-600to800_2024.root"},
  {"QCD 800-1000",  "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-800to1000_2024.root"},
  {"QCD 1000-1200", "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1000-1200_2024.root"},
  {"QCD 1200-1500", "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1200-1500_2024.root"},
  {"QCD 1500-2000", "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1500-2000_2024.root"},
  {"QCD 2000toinf", "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-2000toinf_2024.root"}
};

std::map<int, std::string> signalFiles = {
  {12, "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-12_2024.root"},
  {15, "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-15_2024.root"},
  {20, "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-20_2024.root"},
  {25, "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-25_2024.root"},
  {30, "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-30_2024.root"}
};

// Replace these with the real nominal cross sections in pb.
std::map<int, double> sigmaNom_pb = {
  {12, 1.0},
  {15, 1.0},
  {20, 1.0},
  {25, 1.0},
  {30, 1.0}
};

// ============================================================
// Small helpers
// ============================================================
double safeSigma(double nominal, double rel) {
  const double s = std::fabs(nominal) * rel;
  return (s > 1e-6 ? s : 1.0);
}

void FixNegativeBins(TH1* h) {
  if (!h) return;

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    if (h->GetBinContent(ib) < 0.0) {
      h->SetBinContent(ib, 0.0);
      h->SetBinError(ib, 0.0);
    }
  }
}

double Yield(const TH1* h) {
  if (!h) return 0.0;
  return h->Integral(1, h->GetNbinsX());
}

double Quantile(std::vector<double> v, double q) {
  if (v.empty()) return -1.0;

  std::sort(v.begin(), v.end());

  if (q <= 0.0) return v.front();
  if (q >= 1.0) return v.back();

  const double pos = q * (v.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(pos));
  const size_t hi = static_cast<size_t>(std::ceil(pos));

  if (lo == hi) return v[lo];

  const double frac = pos - lo;
  return v[lo] + frac * (v[hi] - v[lo]);
}

bool ExpectedBinsPositive(const TH1* hSig,
                          const TH1* hBkg,
                          double NexpSig,
                          double NexpBkg) {
  if (!hSig || !hBkg) return false;

  const double sigNorm = hSig->Integral(1, hSig->GetNbinsX());
  const double bkgNorm = hBkg->Integral(1, hBkg->GetNbinsX());

  if (!(sigNorm > 0.0) || !(bkgNorm > 0.0)) return false;

  const int nb = hSig->GetNbinsX();

  for (int ib = 1; ib <= nb; ++ib) {
    const double sFrac = hSig->GetBinContent(ib) / sigNorm;
    const double bFrac = hBkg->GetBinContent(ib) / bkgNorm;

    const double expected = NexpSig * sFrac + NexpBkg * bFrac;

    // RooAbsPdf::getLogVal requires the total p.d.f. to be strictly positive
    // where data can be evaluated. A tiny tolerance avoids numerical noise.
    if (expected <= 1e-12) return false;
  }

  return true;
}

// ============================================================
// Containers
// ============================================================
struct Histograms {
  TH1* h_sig_raw = nullptr;
  TH1* h_bkg_raw = nullptr;

  TH1* h_sig_r = nullptr;
  TH1* h_bkg_r = nullptr;

  double Nsig = 0.0;
  double Nbkg = 0.0;

  void cleanup() {
    delete h_sig_raw;
    delete h_bkg_raw;
    delete h_sig_r;
    delete h_bkg_r;

    h_sig_raw = nullptr;
    h_bkg_raw = nullptr;
    h_sig_r = nullptr;
    h_bkg_r = nullptr;
  }
};

struct AnalysisModel {
  RooRealVar* output_BDT = nullptr;
  RooBinning* customBinning = nullptr;

  RooDataHist* dh_sig = nullptr;
  RooDataHist* dh_bkg = nullptr;

  RooHistPdf* pdf_sig = nullptr;
  RooHistPdf* pdf_bkg = nullptr;

  RooRealVar* Nexp_sig = nullptr;
  RooRealVar* Nexp_bkg = nullptr;

  RooConstVar* bkg_nom = nullptr;
  RooConstVar* bkg_sig = nullptr;
  RooGaussian* c_bkg = nullptr;

  RooAddPdf* model_0 = nullptr;
  RooAddPdf* model_1 = nullptr;

  RooProdPdf* total_model_0 = nullptr;
  RooProdPdf* total_model_1 = nullptr;

  void cleanup() {
    delete total_model_1;
    delete total_model_0;

    delete model_1;
    delete model_0;

    delete c_bkg;
    delete bkg_sig;
    delete bkg_nom;

    delete Nexp_bkg;
    delete Nexp_sig;

    delete pdf_bkg;
    delete pdf_sig;

    delete dh_bkg;
    delete dh_sig;

    delete customBinning;
    delete output_BDT;

    total_model_1 = nullptr;
    total_model_0 = nullptr;

    model_1 = nullptr;
    model_0 = nullptr;

    c_bkg = nullptr;
    bkg_sig = nullptr;
    bkg_nom = nullptr;

    Nexp_bkg = nullptr;
    Nexp_sig = nullptr;

    pdf_bkg = nullptr;
    pdf_sig = nullptr;

    dh_bkg = nullptr;
    dh_sig = nullptr;

    customBinning = nullptr;
    output_BDT = nullptr;
  }
};

struct LimitSummary {
  std::vector<double> mass;

  std::vector<double> NexpSig_obs;
  std::vector<double> NexpSig_exp_m2;
  std::vector<double> NexpSig_exp_m1;
  std::vector<double> NexpSig_exp_med;
  std::vector<double> NexpSig_exp_p1;
  std::vector<double> NexpSig_exp_p2;

  std::vector<double> mu_obs;
  std::vector<double> mu_exp_m2;
  std::vector<double> mu_exp_m1;
  std::vector<double> mu_exp_med;
  std::vector<double> mu_exp_p1;
  std::vector<double> mu_exp_p2;

  std::vector<double> sigma_obs;
  std::vector<double> sigma_exp_m2;
  std::vector<double> sigma_exp_m1;
  std::vector<double> sigma_exp_med;
  std::vector<double> sigma_exp_p1;
  std::vector<double> sigma_exp_p2;

  std::vector<double> NsigTheory;
  std::vector<double> NbkgTheory;
};

// ============================================================
// Histogram loading and rebinning
// ============================================================
TH1* GetClonedHist(TFile* f, const char* hname, const char* newname) {
  if (!f) return nullptr;

  TH1* h = dynamic_cast<TH1*>(f->Get(hname));
  if (!h) return nullptr;

  TH1* c = dynamic_cast<TH1*>(h->Clone(newname));
  if (!c) return nullptr;

  c->SetDirectory(nullptr);
  return c;
}

TH1* BuildSingleRawTemplate(const char* fileName,
                            const char* histName,
                            const char* outName) {
  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] cannot open file " << fileName << "\n";
    if (f) delete f;
    return nullptr;
  }

  TH1* h = GetClonedHist(f, histName, outName);

  f->Close();
  delete f;

  if (!h) {
    std::cerr << "[ERROR] histogram " << histName
              << " not found in " << fileName << "\n";
    return nullptr;
  }

  FixNegativeBins(h);
  return h;
}

TH1* BuildSummedRawTemplate(const std::vector<Sample>& samples,
                            const char* histName,
                            const char* outName) {
  TH1* hSum = nullptr;

  for (size_t i = 0; i < samples.size(); ++i) {
    const auto& s = samples[i];

    TFile* f = TFile::Open(s.fileName.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "[ERROR] cannot open file " << s.fileName << "\n";
      if (f) delete f;
      delete hSum;
      return nullptr;
    }

    TH1* h = GetClonedHist(f, histName, Form("%s_comp_%zu", outName, i));

    f->Close();
    delete f;

    if (!h) {
      std::cerr << "[ERROR] histogram " << histName
                << " not found in " << s.fileName << "\n";
      delete hSum;
      return nullptr;
    }

    FixNegativeBins(h);

    if (!hSum) {
      hSum = dynamic_cast<TH1*>(h->Clone(outName));
      hSum->SetDirectory(nullptr);
      hSum->Reset("ICES");
    }

    hSum->Add(h);
    delete h;
  }

  return hSum;
}

bool BuildInputHistograms(int mass, Histograms& hs) {
  hs.h_bkg_raw = BuildSummedRawTemplate(backgrounds,
                                        opt.histName,
                                        Form("h_bkg_raw_m%d", mass));

  if (!hs.h_bkg_raw) return false;

  auto it = signalFiles.find(mass);
  if (it == signalFiles.end()) {
    std::cerr << "[ERROR] missing signal file for mass " << mass << "\n";
    return false;
  }

  hs.h_sig_raw = BuildSingleRawTemplate(it->second.c_str(),
                                        opt.histName,
                                        Form("h_sig_raw_m%d", mass));

  if (!hs.h_sig_raw) return false;

  hs.h_bkg_r = hs.h_bkg_raw->Rebin(opt.nbins,
                                   Form("h_bkg_r_m%d", mass),
                                   opt.var_bins);

  hs.h_sig_r = hs.h_sig_raw->Rebin(opt.nbins,
                                   Form("h_sig_r_m%d", mass),
                                   opt.var_bins);

  if (!hs.h_bkg_r || !hs.h_sig_r) {
    std::cerr << "[ERROR] variable rebin failed for mass " << mass << "\n";
    return false;
  }

  hs.h_bkg_r->SetDirectory(nullptr);
  hs.h_sig_r->SetDirectory(nullptr);

  FixNegativeBins(hs.h_bkg_r);
  FixNegativeBins(hs.h_sig_r);

  hs.Nsig = Yield(hs.h_sig_r);
  hs.Nbkg = Yield(hs.h_bkg_r);

  if (!(hs.Nsig > 0.0) || !(hs.Nbkg > 0.0)) {
    std::cerr << "[ERROR] non-positive yields for mass " << mass
              << " Nsig=" << hs.Nsig
              << " Nbkg=" << hs.Nbkg << "\n";
    return false;
  }

  return true;
}

// ============================================================
// 1) Construct RooFit model, in the style you requested
// ============================================================
bool ConstructModel(const Histograms& hs, AnalysisModel& fm) {
  fm.output_BDT = new RooRealVar("output_BDT", "BDT score", -1.0, 1.0);

  fm.customBinning = new RooBinning(opt.nbins, opt.var_bins);
  fm.output_BDT->setBinning(*fm.customBinning, "customBinning");

  fm.dh_sig = new RooDataHist("sig",
                              "sig",
                              RooArgList(*fm.output_BDT),
                              Import(*hs.h_sig_r));

  fm.dh_bkg = new RooDataHist("bkg",
                              "bkg",
                              RooArgList(*fm.output_BDT),
                              Import(*hs.h_bkg_r));

  fm.pdf_sig = new RooHistPdf("sig_pdf",
                              "sig_pdf",
                              RooArgSet(*fm.output_BDT),
                              *fm.dh_sig);

  fm.pdf_bkg = new RooHistPdf("bkg_pdf",
                              "bkg_pdf",
                              RooArgSet(*fm.output_BDT),
                              *fm.dh_bkg);

  fm.Nexp_sig = new RooRealVar("Nexp_sig",
                               "Expected signal events",
                               hs.Nsig,
                               -30.0 * std::max(1.0, hs.Nsig),
                                30.0 * std::max(1.0, hs.Nsig));

  fm.Nexp_bkg = new RooRealVar("Nexp_bkg",
                               "Expected combined background",
                               hs.Nbkg,
                               -0.5 * hs.Nbkg,
                                2.0 * hs.Nbkg);

  fm.bkg_nom = new RooConstVar("bkg_nom", "bkg_nom", hs.Nbkg);
  fm.bkg_sig = new RooConstVar("bkg_sig", "bkg_sig", safeSigma(hs.Nbkg, opt.rel_bkg));

  fm.c_bkg = new RooGaussian("c_bkg",
                             "constraint bkg",
                             *fm.Nexp_bkg,
                             *fm.bkg_nom,
                             *fm.bkg_sig);

  fm.model_0 = new RooAddPdf("model_0",
                             "Background-only (extended)",
                             RooArgList(*fm.pdf_bkg),
                             RooArgList(*fm.Nexp_bkg));

  fm.model_1 = new RooAddPdf("model_1",
                             "Signal+Background (extended)",
                             RooArgList(*fm.pdf_sig, *fm.pdf_bkg),
                             RooArgList(*fm.Nexp_sig, *fm.Nexp_bkg));

  // Optional constrained versions.
  // Switch to total_model_0/1 in fits and posterior if you want the Gaussian
  // background constraint included explicitly in the likelihood.
  fm.total_model_0 = new RooProdPdf("total_model_0",
                                    "b-only with constraints",
                                    RooArgList(*fm.model_0, *fm.c_bkg));

  fm.total_model_1 = new RooProdPdf("total_model_1",
                                    "s+b with constraints",
                                    RooArgList(*fm.model_1, *fm.c_bkg));

  return true;
}

// ============================================================
// 2) Generate pseudodata
// ============================================================
RooDataHist* GenerateToyBinned(AnalysisModel& fm,
                               TRandom3& rng,
                               double nBkgTrue,
                               const char* name) {
  const int nEvents = rng.PoissonD(nBkgTrue);

  fm.Nexp_sig->setVal(0.0);
  fm.Nexp_bkg->setVal(nBkgTrue);

  std::unique_ptr<RooDataSet> toyUnbinned(
    fm.model_0->generate(*fm.output_BDT,
                         nEvents,
                         RooFit::Name(Form("%s_unb", name)))
  );

  if (!toyUnbinned) return nullptr;

  return new RooDataHist(name,
                         name,
                         RooArgList(*fm.output_BDT),
                         *toyUnbinned);
}

// ============================================================
// Fit diagnostics
// ============================================================
RooFitResult* FitBOnlyModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.model_0->fitTo(data,
                           Save(),
                           Extended(kTRUE),
                           SumW2Error(kFALSE),
                           RooFit::MaxCalls(10000),
                           Strategy(1),
                           RooFit::Optimize(kTRUE),
                           PrintLevel(-1));
}

RooFitResult* FitSBModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.model_1->fitTo(data,
                           Save(),
                           Extended(kTRUE),
                           SumW2Error(kFALSE),
                           RooFit::MaxCalls(20000),
                           Strategy(2),
                           RooFit::Optimize(kTRUE),
                           PrintLevel(-1));
}

void PrintOneFitSummary(int mass,
                        const Histograms& hs,
                        AnalysisModel& fm,
                        RooFitResult* fitres) {
  std::cout << "\n============================================================\n";
  std::cout << " ONE-FIT SUMMARY  m_A = " << mass << " GeV\n";
  std::cout << "============================================================\n";

  if (!fitres) {
    std::cout << "[ERROR] null fit result\n";
    return;
  }

  std::cout << "status   = " << fitres->status()  << "\n";
  std::cout << "covQual  = " << fitres->covQual() << "\n";
  std::cout << "EDM      = " << fitres->edm()     << "\n";
  std::cout << "minNll   = " << fitres->minNll()  << "\n";

  std::cout << "\nNominal values:\n";
  std::cout << "NsigTheory = " << hs.Nsig << "\n";
  std::cout << "NbkgTheory = " << hs.Nbkg << "\n";

  std::cout << "\nFit result:\n";
  std::cout << "Nexp_bkg = " << fm.Nexp_bkg->getVal()
            << " +/- " << fm.Nexp_bkg->getError() << "\n";

  std::cout << "\nFinal parameters:\n";
  fitres->floatParsFinal().Print("v");

  std::cout << "\nCorrelation matrix:\n";
  fitres->correlationMatrix().Print();

  std::cout << "============================================================\n";
}

void SaveCorrelationPlot(RooFitResult* fitres,
                         const TString& outname,
                         const char* title) {
  if (!fitres) return;

  TH2* hcorr = fitres->correlationHist();
  if (!hcorr) return;

  TCanvas* c = new TCanvas(Form("c_corr_%s", outname.Data()), title, 900, 800);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);

  hcorr->SetTitle(title);
  hcorr->GetZaxis()->SetRangeUser(-1.0, 1.0);
  hcorr->LabelsOption("v", "X");
  hcorr->SetMarkerSize(1.5);

  gStyle->SetPaintTextFormat(".2f");
  hcorr->Draw("COLZ TEXT");

  c->SaveAs(outname);

  delete c;
}

void PlotOneToyFit(int mass,
                   AnalysisModel& fm,
                   RooDataHist& toyData) {
  TCanvas* c = new TCanvas(Form("c_onefit_m%d", mass),
                           "Pseudodata + fit",
                           900,
                           800);

  RooPlot* fr = fm.output_BDT->frame();

  toyData.plotOn(fr, Name("toydata"));
  fm.model_0->plotOn(fr, Name("fullfit"), LineColor(kBlue));
  fm.model_0->plotOn(fr,
                     Components(fm.pdf_bkg->GetName()),
                     LineColor(kRed + 1),
                     LineStyle(kSolid),
                     Name("bkgcomp"));

  fr->SetTitle(Form("One pseudodata fit, m_{A}=%d GeV", mass));
  fr->GetXaxis()->SetTitle("BDT score");
  fr->GetYaxis()->SetTitle("Events");

  c->cd();
  c->SetLogy(0);
  fr->SetMinimum(0.0);
  fr->Draw();

  TLegend* leg = new TLegend(0.56, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fr->findObject("toydata"), "Pseudodata", "pe");
  leg->AddEntry(fr->findObject("fullfit"), "B-only fit", "l");
  leg->AddEntry(fr->findObject("bkgcomp"), "QCD background", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, Form("m_{A}=%d GeV", mass));

  c->SaveAs(Form("fits_and_limit_plots/step3_onefit_bonly_bdt_ma%d.pdf", mass));

  c->SetLogy(1);
  fr->SetMinimum(0.1);
  fr->Draw();
  leg->Draw();
  lat.DrawLatex(0.15, 0.92, Form("m_{A}=%d GeV", mass));

  c->SaveAs(Form("fits_and_limit_plots/step3_onefit_bonly_bdt_ma%d_log.pdf", mass));

  delete leg;
  delete fr;
  delete c;
}

// ============================================================
// Optional pull study
// ============================================================
void RunToyPullStudy(int mass,
                     const Histograms& hs,
                     AnalysisModel& fm,
                     TRandom3& rng) {
  TH1D* h_pull_bkg = new TH1D(Form("h_pull_bkg_m%d", mass),
                              "Background pull;Pull;Toys",
                              60,
                              -5,
                              5);

  TH1D* h_fit_bkg = new TH1D(Form("h_fit_bkg_m%d", mass),
                             "Fitted Nexp_bkg;Nexp_bkg;Toys",
                             60,
                             0.0,
                             2.0 * hs.Nbkg);

  TH1D* h_err_bkg = new TH1D(Form("h_err_bkg_m%d", mass),
                             "Fit error on Nexp_bkg;#sigma(Nexp_bkg);Toys",
                             60,
                             0.0,
                             hs.Nbkg);

  int nFail = 0;

  for (int it = 0; it < opt.Ntoys_pulls; ++it) {
    std::unique_ptr<RooDataHist> toy(
      GenerateToyBinned(fm, rng, hs.Nbkg, Form("toy_pull_m%d_%d", mass, it))
    );

    if (!toy) continue;

    fm.Nexp_bkg->setVal(hs.Nbkg);

    std::unique_ptr<RooFitResult> fitres(
      FitBOnlyModel(fm, *toy)
    );

    if (!fitres || fitres->status() != 0 || fitres->covQual() < 2) {
      ++nFail;
      continue;
    }

    const double val = fm.Nexp_bkg->getVal();
    const double err = fm.Nexp_bkg->getError();

    if (err > 0.0) h_pull_bkg->Fill((hs.Nbkg - val) / err);

    h_fit_bkg->Fill(val);
    h_err_bkg->Fill(err);
  }

  std::cout << "[PULLS] mass=" << mass
            << " failed fits = " << nFail
            << " / " << opt.Ntoys_pulls << "\n";

  std::cout << "[PULLS] mean=" << h_pull_bkg->GetMean()
            << " RMS=" << h_pull_bkg->GetRMS() << "\n";

  TCanvas* c1 = new TCanvas(Form("c_pull_bkg_m%d", mass),
                            "Background pull",
                            800,
                            600);
  h_pull_bkg->Fit("gaus", "Q", "");
  h_pull_bkg->Draw();
  c1->SaveAs(Form("fits_and_limit_plots/step4_pull_distribution_bkg_ma%d.pdf", mass));
  delete c1;

  TCanvas* c2 = new TCanvas(Form("c_fit_bkg_m%d", mass),
                            "Fitted background",
                            800,
                            600);
  h_fit_bkg->Draw();
  c2->SaveAs(Form("fits_and_limit_plots/step4_fitval_bkg_ma%d.pdf", mass));
  delete c2;

  TCanvas* c3 = new TCanvas(Form("c_err_bkg_m%d", mass),
                            "Background fit error",
                            800,
                            600);
  h_err_bkg->Draw();
  c3->SaveAs(Form("fits_and_limit_plots/step4_fiterr_bkg_ma%d.pdf", mass));
  delete c3;

  delete h_pull_bkg;
  delete h_fit_bkg;
  delete h_err_bkg;
}

// ============================================================
// Deterministic Bayesian posterior in Nexp_sig,
// marginalized over Nexp_bkg.
// Nexp_sig is allowed to be negative in the scan and in the posterior plot.
// ============================================================
bool ComputePosteriorNexpSigMarginalized(AnalysisModel& fm,
                                         RooDataHist& data,
                                         const Histograms& hs,
                                         double NexpSigMinPhysical,
                                         double NexpSigMaxPhysical,
                                         double NexpSigPlotMin,
                                         double NexpSigPlotMax,
                                         double NexpBkgMin,
                                         double NexpBkgMax,
                                         TH1D*& hPosterior,
                                         double& NexpSig95) {
  hPosterior = nullptr;
  NexpSig95 = -1.0;

  if (opt.nNexpSigScan < 2 || opt.nNexpBkgScan < 2) return false;
  if (NexpSigMaxPhysical <= NexpSigMinPhysical) return false;
  if (NexpSigPlotMax <= NexpSigPlotMin) return false;
  if (NexpBkgMax <= NexpBkgMin) return false;

  std::unique_ptr<RooAbsReal> nll(
    fm.model_1->createNLL(data, RooFit::Extended(true))
  );

  if (!nll) return false;

  std::vector<double> sigVals(opt.nNexpSigScan, 0.0);
  std::vector<double> postSig(opt.nNexpSigScan, 0.0);

  const double dSig =
    (NexpSigMaxPhysical - NexpSigMinPhysical) / (opt.nNexpSigScan - 1);

  const double dBkg =
    (NexpBkgMax - NexpBkgMin) / (opt.nNexpBkgScan - 1);

  for (int i = 0; i < opt.nNexpSigScan; ++i) {
    sigVals[i] = NexpSigMinPhysical + i * dSig;
  }

  double minNll = 1e100;

  // Pass 1: global minimum in physical signal region.
  for (int iSig = 0; iSig < opt.nNexpSigScan; ++iSig) {
    fm.Nexp_sig->setVal(sigVals[iSig]);
    fm.Nexp_sig->setConstant(true);

    for (int iB = 0; iB < opt.nNexpBkgScan; ++iB) {
      const double nb = NexpBkgMin + iB * dBkg;

      if (!ExpectedBinsPositive(hs.h_sig_r, hs.h_bkg_r, sigVals[iSig], nb)) {
        continue;
      }

      fm.Nexp_bkg->setVal(nb);

      const double val = nll->getVal();

      if (std::isfinite(val) && val < minNll) {
        minNll = val;
      }
    }
  }

  if (!std::isfinite(minNll)) {
    fm.Nexp_sig->setConstant(false);
    return false;
  }

  // Pass 2: marginalize over Nexp_bkg.
  for (int iSig = 0; iSig < opt.nNexpSigScan; ++iSig) {
    fm.Nexp_sig->setVal(sigVals[iSig]);
    fm.Nexp_sig->setConstant(true);

    double sumB = 0.0;

    for (int iB = 0; iB < opt.nNexpBkgScan; ++iB) {
      const double nb = NexpBkgMin + iB * dBkg;

      if (!ExpectedBinsPositive(hs.h_sig_r, hs.h_bkg_r, sigVals[iSig], nb)) {
        continue;
      }

      fm.Nexp_bkg->setVal(nb);

      const double val = nll->getVal();
      if (!std::isfinite(val)) continue;

      const double w = std::exp(-(val - minNll));
      if (std::isfinite(w)) sumB += w;
    }

    postSig[iSig] = sumB;
  }

  fm.Nexp_sig->setConstant(false);

  const double norm =
    std::accumulate(postSig.begin(), postSig.end(), 0.0) * dSig;

  if (!(norm > 0.0)) return false;

  hPosterior = new TH1D(Form("hPosteriorNexpSig_%s", data.GetName()),
                        ";N_{exp}^{sig};p(N_{exp}^{sig} | data)",
                        500,
                        NexpSigPlotMin,
                        NexpSigPlotMax);

  for (int ib = 1; ib <= hPosterior->GetNbinsX(); ++ib) {
    hPosterior->SetBinContent(ib, 0.0);
  }

  double cumulative = 0.0;
  bool found95 = false;

  for (int i = 0; i < opt.nNexpSigScan; ++i) {
    const double p = postSig[i] / norm;

    const int bin = hPosterior->FindBin(sigVals[i]);
    hPosterior->SetBinContent(bin, p);

    cumulative += p * dSig;

    if (!found95 && cumulative >= opt.CL) {
      NexpSig95 = sigVals[i];
      found95 = true;
    }
  }

  if (!found95) NexpSig95 = NexpSigMaxPhysical;

  return true;
}

void SavePosteriorPlotNexpSig(int mass,
                              TH1D* h,
                              double NexpSig95,
                              const TString& outname) {
  if (!h) return;

  TCanvas* c = new TCanvas(Form("c_post_NexpSig_m%d", mass),
                           "Posterior",
                           900,
                           700);

  c->SetTicks(1, 1);

  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(3);
  h->SetTitle("");
  h->Draw("HIST");

  const double ymax = 1.15 * h->GetMaximum();

  TLine* zeroLine = new TLine(0.0, 0.0, 0.0, ymax);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(3);
  zeroLine->SetLineWidth(2);
  zeroLine->Draw("SAME");

  TLine* line95 = new TLine(NexpSig95, 0.0, NexpSig95, ymax);
  line95->SetLineColor(kRed + 1);
  line95->SetLineStyle(2);
  line95->SetLineWidth(3);
  line95->Draw("SAME");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.16, 0.92, Form("Posterior for m_{A}=%d GeV", mass));
  lat.DrawLatex(0.16, 0.87, Form("95%% upper limit: N_{exp}^{sig}=%.4g", NexpSig95));

  TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h, "Posterior", "l");
  leg->AddEntry(zeroLine, "N_{exp}^{sig}=0", "l");
  leg->AddEntry(line95, "95% upper limit", "l");
  leg->Draw();

  c->SaveAs(outname);

  delete leg;
  delete line95;
  delete zeroLine;
  delete c;
}

void SavePosteriorPlotMu(int mass,
                         TH1D* h,
                         double mu95,
                         const TString& outname) {
  if (!h) return;

  TCanvas* c = new TCanvas(Form("c_post_mu_m%d", mass),
                           "Posterior in mu",
                           900,
                           700);

  c->SetTicks(1, 1);

  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(3);
  h->SetTitle("");
  h->Draw("HIST");

  const double ymax = 1.15 * h->GetMaximum();

  TLine* zeroLine = new TLine(0.0, 0.0, 0.0, ymax);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(3);
  zeroLine->SetLineWidth(2);
  zeroLine->Draw("SAME");

  TLine* line95 = new TLine(mu95, 0.0, mu95, ymax);
  line95->SetLineColor(kRed + 1);
  line95->SetLineStyle(2);
  line95->SetLineWidth(3);
  line95->Draw("SAME");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.16, 0.92, Form("Posterior for m_{A}=%d GeV", mass));
  lat.DrawLatex(0.16, 0.87, Form("95%% upper limit: #mu=%.4g", mu95));

  TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h, "Posterior", "l");
  leg->AddEntry(zeroLine, "#mu=0", "l");
  leg->AddEntry(line95, "95% upper limit", "l");
  leg->Draw();

  c->SaveAs(outname);

  delete leg;
  delete line95;
  delete zeroLine;
  delete c;
}

TH1D* BuildPosteriorMuFromNexpSig(const TH1D* hNexp,
                                  double NsigTheory,
                                  int mass,
                                  const char* tag) {
  if (!hNexp || !(NsigTheory > 0.0)) return nullptr;

  const double muMinPlot = hNexp->GetXaxis()->GetXmin() / NsigTheory;
  const double muMaxPlot = hNexp->GetXaxis()->GetXmax() / NsigTheory;

  TH1D* hMu = new TH1D(Form("hPosteriorMu_%s_m%d", tag, mass),
                       ";#mu;p(#mu | data)",
                       hNexp->GetNbinsX(),
                       muMinPlot,
                       muMaxPlot);

  for (int ib = 1; ib <= hNexp->GetNbinsX(); ++ib) {
    const double Ncenter = hNexp->GetBinCenter(ib);
    const double muCenter = Ncenter / NsigTheory;

    // If hNexp stores p(N), then p(mu) = p(N) * dN/dmu = p(N) * NsigTheory.
    // The shape and cumulative limits are identical after this linear transformation.
    const double pMu = hNexp->GetBinContent(ib) * NsigTheory;

    const int mbin = hMu->FindBin(muCenter);
    hMu->SetBinContent(mbin, pMu);
  }

  return hMu;
}

// ============================================================
// Final Brazil-style plots
// ============================================================
void MakeFinalLimitPlot(const std::vector<double>& mass,
                        const std::vector<double>& obs,
                        const std::vector<double>& med,
                        const std::vector<double>& m1,
                        const std::vector<double>& p1,
                        const std::vector<double>& m2,
                        const std::vector<double>& p2,
                        const char* ytitle,
                        const char* outbase) {
  const int n = static_cast<int>(mass.size());
  if (n <= 0) return;

  std::vector<double> x(n), y(n), yobs(n);
  std::vector<double> exl(n, 0.0), exh(n, 0.0);
  std::vector<double> eyl1(n), eyh1(n), eyl2(n), eyh2(n);

  for (int i = 0; i < n; ++i) {
    x[i] = mass[i];
    y[i] = med[i];
    yobs[i] = obs[i];

    eyl1[i] = med[i] - m1[i];
    eyh1[i] = p1[i] - med[i];

    eyl2[i] = med[i] - m2[i];
    eyh2[i] = p2[i] - med[i];
  }

  TCanvas* c = new TCanvas(Form("c_%s", outbase),
                           "Expected Bayesian limits",
                           900,
                           700);

  TGraphAsymmErrors* g2 =
    new TGraphAsymmErrors(n, &x[0], &y[0],
                          &exl[0], &exh[0],
                          &eyl2[0], &eyh2[0]);

  TGraphAsymmErrors* g1 =
    new TGraphAsymmErrors(n, &x[0], &y[0],
                          &exl[0], &exh[0],
                          &eyl1[0], &eyh1[0]);

  TGraph* gmed = new TGraph(n, &x[0], &y[0]);
  TGraph* gobs = new TGraph(n, &x[0], &yobs[0]);

  g2->SetFillColor(kYellow);
  g2->SetLineColor(kYellow);
  g2->SetTitle(Form(";m_{A} [GeV];%s", ytitle));

  g1->SetFillColor(kGreen + 1);
  g1->SetLineColor(kGreen + 1);

  gmed->SetLineColor(kBlack);
  gmed->SetLineWidth(2);
  gmed->SetLineStyle(2);
  gmed->SetMarkerStyle(20);

  gobs->SetLineColor(kRed + 1);
  gobs->SetLineWidth(2);
  gobs->SetMarkerStyle(24);

  g2->SetMinimum(0.0);
  g2->Draw("A3");
  g1->Draw("3 SAME");
  gmed->Draw("LP SAME");
  gobs->Draw("LP SAME");

  TLegend* leg = new TLegend(0.55, 0.64, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gmed, "Median expected", "lp");
  leg->AddEntry(gobs, "Pseudo-observed", "lp");
  leg->AddEntry(g1, "#pm1#sigma", "f");
  leg->AddEntry(g2, "#pm2#sigma", "f");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");

  c->SaveAs(Form("fits_and_limit_plots/%s.pdf", outbase));

  c->Clear();
  c->SetLogy();

  g2->SetMinimum(1e-4);
  g2->Draw("A3");
  g1->Draw("3 SAME");
  gmed->Draw("LP SAME");
  gobs->Draw("LP SAME");
  leg->Draw();
  lat.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");

  c->SaveAs(Form("fits_and_limit_plots/%s_log.pdf", outbase));

  delete leg;
  delete gobs;
  delete gmed;
  delete g1;
  delete g2;
  delete c;
}

// ============================================================
// Main
// ============================================================
void BDTToyMC_Study_ModeSwitch() {
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  gSystem->mkdir("fits_and_limit_plots", kTRUE);

  TFile* fout = new TFile("BDTToyMC_BayesianLimits.root", "RECREATE");

  TTree* tree = new TTree("LimitTree", "Bayesian limit summary");

  int out_mass = -1;
  double out_NsigTheory = -1.0;
  double out_NbkgTheory = -1.0;

  double out_NexpSig_obs = -1.0;
  double out_NexpSig_m2 = -1.0;
  double out_NexpSig_m1 = -1.0;
  double out_NexpSig_med = -1.0;
  double out_NexpSig_p1 = -1.0;
  double out_NexpSig_p2 = -1.0;

  double out_mu_obs = -1.0;
  double out_mu_m2 = -1.0;
  double out_mu_m1 = -1.0;
  double out_mu_med = -1.0;
  double out_mu_p1 = -1.0;
  double out_mu_p2 = -1.0;

  double out_sigma_obs = -1.0;
  double out_sigma_m2 = -1.0;
  double out_sigma_m1 = -1.0;
  double out_sigma_med = -1.0;
  double out_sigma_p1 = -1.0;
  double out_sigma_p2 = -1.0;

  tree->Branch("mass", &out_mass, "mass/I");
  tree->Branch("NsigTheory", &out_NsigTheory, "NsigTheory/D");
  tree->Branch("NbkgTheory", &out_NbkgTheory, "NbkgTheory/D");

  tree->Branch("NexpSig_obs", &out_NexpSig_obs, "NexpSig_obs/D");
  tree->Branch("NexpSig_exp_m2", &out_NexpSig_m2, "NexpSig_exp_m2/D");
  tree->Branch("NexpSig_exp_m1", &out_NexpSig_m1, "NexpSig_exp_m1/D");
  tree->Branch("NexpSig_exp_med", &out_NexpSig_med, "NexpSig_exp_med/D");
  tree->Branch("NexpSig_exp_p1", &out_NexpSig_p1, "NexpSig_exp_p1/D");
  tree->Branch("NexpSig_exp_p2", &out_NexpSig_p2, "NexpSig_exp_p2/D");

  tree->Branch("mu_obs", &out_mu_obs, "mu_obs/D");
  tree->Branch("mu_exp_m2", &out_mu_m2, "mu_exp_m2/D");
  tree->Branch("mu_exp_m1", &out_mu_m1, "mu_exp_m1/D");
  tree->Branch("mu_exp_med", &out_mu_med, "mu_exp_med/D");
  tree->Branch("mu_exp_p1", &out_mu_p1, "mu_exp_p1/D");
  tree->Branch("mu_exp_p2", &out_mu_p2, "mu_exp_p2/D");

  tree->Branch("sigma_obs", &out_sigma_obs, "sigma_obs/D");
  tree->Branch("sigma_exp_m2", &out_sigma_m2, "sigma_exp_m2/D");
  tree->Branch("sigma_exp_m1", &out_sigma_m1, "sigma_exp_m1/D");
  tree->Branch("sigma_exp_med", &out_sigma_med, "sigma_exp_med/D");
  tree->Branch("sigma_exp_p1", &out_sigma_p1, "sigma_exp_p1/D");
  tree->Branch("sigma_exp_p2", &out_sigma_p2, "sigma_exp_p2/D");

  LimitSummary summary;

  for (int mass : opt.masses) {
    std::cout << "\n\n############################################################\n";
    std::cout << " Running mass m_A = " << mass << " GeV\n";
    std::cout << "############################################################\n";

    Histograms hs;
    AnalysisModel fm;

    std::cout << "[STEP 1] Build histograms\n";

    if (!BuildInputHistograms(mass, hs)) {
      std::cerr << "[ERROR] failed to build input histograms for mass "
                << mass << "\n";
      hs.cleanup();
      continue;
    }

    std::cout << "[INFO] NsigTheory = " << hs.Nsig << "\n";
    std::cout << "[INFO] NbkgTheory = " << hs.Nbkg << "\n";

    std::cout << "[STEP 2] Build RooFit model\n";

    if (!ConstructModel(hs, fm)) {
      std::cerr << "[ERROR] failed to construct model for mass "
                << mass << "\n";
      fm.cleanup();
      hs.cleanup();
      continue;
    }

    TRandom3 rng(opt.seedBase + 100 * mass);

    if (opt.makeOneFitPlots) {
      std::cout << "[STEP 3] One toy fit and diagnostic plots\n";

      std::unique_ptr<RooDataHist> oneToy(
        GenerateToyBinned(fm, rng, hs.Nbkg, Form("oneToy_m%d", mass))
      );

      if (oneToy) {
        fm.Nexp_sig->setVal(0.0);
        fm.Nexp_bkg->setVal(hs.Nbkg);

        std::unique_ptr<RooFitResult> fitOne(
          FitBOnlyModel(fm, *oneToy)
        );

        PrintOneFitSummary(mass, hs, fm, fitOne.get());

        PlotOneToyFit(mass, fm, *oneToy);

        SaveCorrelationPlot(
          fitOne.get(),
          Form("fits_and_limit_plots/step3_correlation_matrix_bonly_ma%d.pdf", mass),
          Form("Correlation matrix, m_{A}=%d GeV", mass)
        );
      }
    }

    if (opt.runPullStudy) {
      std::cout << "[STEP 4] Toy MC pull study\n";
      RunToyPullStudy(mass, hs, fm, rng);
    } else {
      std::cout << "[STEP 4] Toy MC pull study skipped\n";
    }

    std::cout << "[STEP 5] Bayesian limits in Nexp_sig\n";

    const double NexpSigMinPhysical = opt.NexpSigMinPhysicalFactor * hs.Nsig;
    const double NexpSigMaxPhysical = opt.NexpSigMaxPhysicalFactor * hs.Nsig;

    const double NexpSigPlotMin = opt.NexpSigPlotMinFactor * hs.Nsig;
    const double NexpSigPlotMax = opt.NexpSigPlotMaxFactor * hs.Nsig;

    const double NexpBkgMin = opt.NexpBkgMinFactor * hs.Nbkg;
    const double NexpBkgMax = opt.NexpBkgMaxFactor * hs.Nbkg;

    std::unique_ptr<RooDataHist> dataObs(
      GenerateToyBinned(fm, rng, hs.Nbkg, Form("pseudoObs_m%d", mass))
    );

    if (!dataObs) {
      std::cerr << "[ERROR] failed to generate pseudo-observed data\n";
      fm.cleanup();
      hs.cleanup();
      continue;
    }

    TH1D* hPosteriorObs = nullptr;
    double NexpSig95_obs = -1.0;

    bool okObs =
      ComputePosteriorNexpSigMarginalized(fm,
                                           *dataObs,
                                           hs,
                                           NexpSigMinPhysical,
                                           NexpSigMaxPhysical,
                                           NexpSigPlotMin,
                                           NexpSigPlotMax,
                                           NexpBkgMin,
                                           NexpBkgMax,
                                           hPosteriorObs,
                                           NexpSig95_obs);

    if (!okObs || !hPosteriorObs) {
      std::cerr << "[ERROR] failed observed posterior for mass "
                << mass << "\n";
      delete hPosteriorObs;
      fm.cleanup();
      hs.cleanup();
      continue;
    }

    SavePosteriorPlotNexpSig(
      mass,
      hPosteriorObs,
      NexpSig95_obs,
      Form("fits_and_limit_plots/step5_posterior_NexpSig_ma%d.pdf", mass)
    );

    TH1D* hPosteriorMuObs = BuildPosteriorMuFromNexpSig(hPosteriorObs,
                                                        hs.Nsig,
                                                        mass,
                                                        "obs");

    const double mu95_obs_for_plot = NexpSig95_obs / hs.Nsig;

    SavePosteriorPlotMu(
      mass,
      hPosteriorMuObs,
      mu95_obs_for_plot,
      Form("fits_and_limit_plots/step5_posterior_mu_ma%d.pdf", mass)
    );

    fout->cd();
    if (hPosteriorMuObs) {
      hPosteriorMuObs->Write(Form("hPosteriorMuObs_m%d", mass));
    }

    delete hPosteriorMuObs;

    std::vector<double> NexpSig95_toys;
    NexpSig95_toys.reserve(opt.Ntoys_limits);

    for (int itoy = 0; itoy < opt.Ntoys_limits; ++itoy) {
      if ((itoy + 1) % 50 == 0) {
        std::cout << "Toy " << (itoy + 1)
                  << " / " << opt.Ntoys_limits << "\n";
      }

      std::unique_ptr<RooDataHist> toyData(
        GenerateToyBinned(fm, rng, hs.Nbkg, Form("toy_m%d_%d", mass, itoy))
      );

      if (!toyData) continue;

      TH1D* hPosteriorToy = nullptr;
      double NexpSig95_toy = -1.0;

      bool okToy =
        ComputePosteriorNexpSigMarginalized(fm,
                                             *toyData,
                                             hs,
                                             NexpSigMinPhysical,
                                             NexpSigMaxPhysical,
                                             NexpSigPlotMin,
                                             NexpSigPlotMax,
                                             NexpBkgMin,
                                             NexpBkgMax,
                                             hPosteriorToy,
                                             NexpSig95_toy);

      delete hPosteriorToy;

      if (!okToy) continue;
      if (!std::isfinite(NexpSig95_toy) || NexpSig95_toy < 0.0) continue;

      NexpSig95_toys.push_back(NexpSig95_toy);
    }

    if (NexpSig95_toys.empty()) {
      std::cerr << "[ERROR] no valid toys for mass " << mass << "\n";
      delete hPosteriorObs;
      fm.cleanup();
      hs.cleanup();
      continue;
    }

    const double NexpSig_m2  = Quantile(NexpSig95_toys, 0.025);
    const double NexpSig_m1  = Quantile(NexpSig95_toys, 0.160);
    const double NexpSig_med = Quantile(NexpSig95_toys, 0.500);
    const double NexpSig_p1  = Quantile(NexpSig95_toys, 0.840);
    const double NexpSig_p2  = Quantile(NexpSig95_toys, 0.975);

    const double mu_obs = NexpSig95_obs / hs.Nsig;
    const double mu_m2  = NexpSig_m2 / hs.Nsig;
    const double mu_m1  = NexpSig_m1 / hs.Nsig;
    const double mu_med = NexpSig_med / hs.Nsig;
    const double mu_p1  = NexpSig_p1 / hs.Nsig;
    const double mu_p2  = NexpSig_p2 / hs.Nsig;

    const double sigmaNom = sigmaNom_pb[mass];

    const double sigma_obs = mu_obs * sigmaNom;
    const double sigma_m2  = mu_m2 * sigmaNom;
    const double sigma_m1  = mu_m1 * sigmaNom;
    const double sigma_med = mu_med * sigmaNom;
    const double sigma_p1  = mu_p1 * sigmaNom;
    const double sigma_p2  = mu_p2 * sigmaNom;

    summary.mass.push_back(mass);

    summary.NexpSig_obs.push_back(NexpSig95_obs);
    summary.NexpSig_exp_m2.push_back(NexpSig_m2);
    summary.NexpSig_exp_m1.push_back(NexpSig_m1);
    summary.NexpSig_exp_med.push_back(NexpSig_med);
    summary.NexpSig_exp_p1.push_back(NexpSig_p1);
    summary.NexpSig_exp_p2.push_back(NexpSig_p2);

    summary.mu_obs.push_back(mu_obs);
    summary.mu_exp_m2.push_back(mu_m2);
    summary.mu_exp_m1.push_back(mu_m1);
    summary.mu_exp_med.push_back(mu_med);
    summary.mu_exp_p1.push_back(mu_p1);
    summary.mu_exp_p2.push_back(mu_p2);

    summary.sigma_obs.push_back(sigma_obs);
    summary.sigma_exp_m2.push_back(sigma_m2);
    summary.sigma_exp_m1.push_back(sigma_m1);
    summary.sigma_exp_med.push_back(sigma_med);
    summary.sigma_exp_p1.push_back(sigma_p1);
    summary.sigma_exp_p2.push_back(sigma_p2);

    summary.NsigTheory.push_back(hs.Nsig);
    summary.NbkgTheory.push_back(hs.Nbkg);

    out_mass = mass;
    out_NsigTheory = hs.Nsig;
    out_NbkgTheory = hs.Nbkg;

    out_NexpSig_obs = NexpSig95_obs;
    out_NexpSig_m2 = NexpSig_m2;
    out_NexpSig_m1 = NexpSig_m1;
    out_NexpSig_med = NexpSig_med;
    out_NexpSig_p1 = NexpSig_p1;
    out_NexpSig_p2 = NexpSig_p2;

    out_mu_obs = mu_obs;
    out_mu_m2 = mu_m2;
    out_mu_m1 = mu_m1;
    out_mu_med = mu_med;
    out_mu_p1 = mu_p1;
    out_mu_p2 = mu_p2;

    out_sigma_obs = sigma_obs;
    out_sigma_m2 = sigma_m2;
    out_sigma_m1 = sigma_m1;
    out_sigma_med = sigma_med;
    out_sigma_p1 = sigma_p1;
    out_sigma_p2 = sigma_p2;

    tree->Fill();

    fout->cd();
    hPosteriorObs->Write(Form("hPosteriorNexpSigObs_m%d", mass));
    hs.h_sig_r->Write(Form("hSig_m%d", mass));
    hs.h_bkg_r->Write(Form("hBkg_m%d", mass));

    std::cout << "[RESULT] mass=" << mass
              << " NexpSig_obs=" << NexpSig95_obs
              << " NexpSig_exp_med=" << NexpSig_med
              << " mu_exp_med=" << mu_med
              << " sigma_exp_med=" << sigma_med
              << "\n";

    delete hPosteriorObs;

    fm.cleanup();
    hs.cleanup();
  }

  fout->cd();
  tree->Write();

  if (!summary.mass.empty()) {
    MakeFinalLimitPlot(summary.mass,
                       summary.mu_obs,
                       summary.mu_exp_med,
                       summary.mu_exp_m1,
                       summary.mu_exp_p1,
                       summary.mu_exp_m2,
                       summary.mu_exp_p2,
                       "Bayesian 95% upper limit on #mu",
                       "step5_expected_limits_mu");

    MakeFinalLimitPlot(summary.mass,
                       summary.sigma_obs,
                       summary.sigma_exp_med,
                       summary.sigma_exp_m1,
                       summary.sigma_exp_p1,
                       summary.sigma_exp_m2,
                       summary.sigma_exp_p2,
                       "Bayesian 95% upper limit on #sigma [pb]",
                       "step5_expected_limits_sigma");
  }

  fout->Close();

  // Do not delete tree explicitly after closing the ROOT file.
  // This avoids the end-of-job segmentation violation seen previously.
  delete fout;

  std::cout << "\nSaved:\n";
  std::cout << "  BDTToyMC_BayesianLimits.root\n";
  std::cout << "  fits_and_limit_plots/step3_onefit_bonly_bdt_ma*.pdf\n";
  std::cout << "  fits_and_limit_plots/step3_onefit_bonly_bdt_ma*_log.pdf\n";
  std::cout << "  fits_and_limit_plots/step3_correlation_matrix_bonly_ma*.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_posterior_NexpSig_ma*.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_posterior_mu_ma*.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_expected_limits_mu.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_expected_limits_mu_log.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_expected_limits_sigma.pdf\n";
  std::cout << "  fits_and_limit_plots/step5_expected_limits_sigma_log.pdf\n";
}

#endif
