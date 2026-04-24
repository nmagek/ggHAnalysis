// DrawggHOverlay.C
// Overlay signals (M12, M30) + stack ALL backgrounds with THStack and show all in legend.
//
// Usage:
//   root [0] .L DrawggHOverlay.C+
//   root [1] DrawggHOverlay();

#ifndef DRAWGGHOVERLAY_C
#define DRAWGGHOVERLAY_C

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TPad.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace ggHOverlay {

// ---------- helpers ----------
static TH1* GetClone(TFile* f, const char* name, const char* newname)
{
  if (!f) return nullptr;
  TH1* h = dynamic_cast<TH1*>(f->Get(name));
  if (!h) return nullptr;
  TH1* c = dynamic_cast<TH1*>(h->Clone(newname));
  if (!c) return nullptr;
  c->SetDirectory(nullptr);
  return c;
}

static void PrepareHist(TH1* h)
{
  if (!h) return;
  if (h->GetSumw2N() == 0) h->Sumw2(kTRUE);
}

static void ApplyRebin(TH1* h, int r)
{
  if (!h) return;
  if (r > 1) h->Rebin(r);
}

static bool SameBinning(const TH1* a, const TH1* b)
{
  if (!a || !b) return false;
  if (a->GetNbinsX() != b->GetNbinsX()) return false;
  if (a->GetXaxis()->GetXmin() != b->GetXaxis()->GetXmin()) return false;
  if (a->GetXaxis()->GetXmax() != b->GetXaxis()->GetXmax()) return false;
  return true;
}

static double FindSmallestPositive(const TH1* h)
{
  if (!h) return 1e30;
  double m = 1e30;
  for (int b=1; b<=h->GetNbinsX(); ++b) {
    const double v = h->GetBinContent(b);
    if (v > 0.0 && v < m) m = v;
  }
  return m;
}

static void StyleSignal(TH1* h, int color)
{
  if (!h) return;
  h->SetLineWidth(3);
  h->SetFillStyle(0);
  h->SetLineColor(color);
}

static void StyleBkgComponent(TH1* h, int fillColor)
{
  if (!h) return;
  h->SetFillStyle(1001);
  h->SetFillColor(fillColor);
  h->SetLineColor(kBlack);
  h->SetLineWidth(1);
}

static void EnsureDirAndCleanPDFs(const char* outdir)
{
  if (!outdir) return;
  gSystem->mkdir(outdir, kTRUE);

  void* dir = gSystem->OpenDirectory(outdir);
  if (!dir) return;

  const char* ent = nullptr;
  while ((ent = gSystem->GetDirEntry(dir))) {
    std::string n(ent);
    if (n.size() >= 4 && n.substr(n.size()-4) == ".pdf") {
      gSystem->Unlink((std::string(outdir)+"/"+n).c_str());
    }
  }
  gSystem->FreeDirectory(dir);
}

static void SaveCanvasPDF(TCanvas* c, const char* outdir, const char* base)
{
  if (!c) return;
  c->SaveAs((std::string(outdir)+"/"+base+".pdf").c_str());
}

// Build THStack + total background histogram.
// IMPORTANT: we also build the legend entries *in the same order as bkgFiles*,
// and we only add an entry when that component exists for this histogram name.
static void BuildBkgStackAndTotal(const std::vector<TFile*>& bkgFiles,
                                  const std::vector<std::string>& bkgLabels,
                                  const char* hname,
                                  int rebin,
                                  THStack*& stOut,
                                  TH1*& totOut,
                                  std::vector<TH1*>& compsOut,
                                  std::vector<std::string>& compLabelsOut)
{
  stOut = new THStack(Form("st_%s", hname), "");
  totOut = nullptr;
  compsOut.clear();
  compLabelsOut.clear();

  std::vector<int> palette = {
    kGray+1, kAzure-9, kGreen-7, kOrange-2, kMagenta-7, kTeal-7, kPink-6, kYellow-7
  };

  int added = 0;

  for (size_t i=0; i<bkgFiles.size(); ++i) {
    TFile* f = bkgFiles[i];
    if (!f) continue;

    TH1* h = GetClone(f, hname, Form("%s_bkg_%zu", hname, i));
    if (!h) {
      std::cout << "[MISS] " << hname << " in " << f->GetName() << "\n";
      continue;
    }

    PrepareHist(h);
    ApplyRebin(h, rebin);

    if (!totOut) {
      totOut = dynamic_cast<TH1*>(h->Clone(Form("%s_bkgTot", hname)));
      totOut->SetDirectory(nullptr);
      totOut->Reset("ICES");
    } else if (!SameBinning(totOut, h)) {
      std::cout << "[SKIP] binning mismatch for " << hname
                << " in " << f->GetName() << "\n";
      delete h;
      continue;
    }

    StyleBkgComponent(h, palette[added % (int)palette.size()]);
    stOut->Add(h, "HIST");
    totOut->Add(h);

    compsOut.push_back(h);
    const std::string lbl = (i < bkgLabels.size()) ? bkgLabels[i] : (std::string("Bkg_") + std::to_string(i));
    compLabelsOut.push_back(lbl);

    ++added;
  }

  if (added == 0) {
    std::cout << "[WARN] stack EMPTY for " << hname << "\n";
  }
}

static void DrawOne(const char* outdir,
                    const char* hname,
                    const char* title,
                    bool logy,
                    int rebin,
                    TFile* fS12,
                    TFile* fS30,
                    const std::vector<TFile*>& bkgFiles,
                    const std::vector<std::string>& bkgLabels)
{
  // signals
  TH1* s12 = GetClone(fS12, hname, (std::string(hname)+"_sig12").c_str());
  TH1* s30 = GetClone(fS30, hname, (std::string(hname)+"_sig30").c_str());
  if (s12) { PrepareHist(s12); ApplyRebin(s12, rebin); StyleSignal(s12, kRed); }
  if (s30) { PrepareHist(s30); ApplyRebin(s30, rebin); StyleSignal(s30, kGreen+2); }

  // backgrounds
  THStack* st = nullptr;
  TH1* bkgTot = nullptr;
  std::vector<TH1*> bkgComps;
  std::vector<std::string> bkgCompLabels;
  BuildBkgStackAndTotal(bkgFiles, bkgLabels, hname, rebin, st, bkgTot, bkgComps, bkgCompLabels);

  const bool haveBkg = (bkgTot && !bkgComps.empty());
  const bool haveSig = (s12 || s30);

  if (!haveBkg && !haveSig) {
    delete st;
    return;
  }

  // y-range from total background + signals
  double ymax = 0.0;
  if (bkgTot) ymax = std::max(ymax, (double)bkgTot->GetMaximum());
  if (s12)    ymax = std::max(ymax, (double)s12->GetMaximum());
  if (s30)    ymax = std::max(ymax, (double)s30->GetMaximum());
  ymax *= 1.25;
  if (ymax <= 0.0) ymax = 1.0;

  double ymin = 0.0;
  if (logy) {
    double minpos = 1e30;
    minpos = std::min(minpos, FindSmallestPositive(bkgTot));
    minpos = std::min(minpos, FindSmallestPositive(s12));
    minpos = std::min(minpos, FindSmallestPositive(s30));
    ymin = (minpos == 1e30) ? 0.5 : std::min(0.5, minpos);
    if (ymin <= 0.0) ymin = 0.5;
  }

  // base axis hist
  TH1* base = bkgTot ? bkgTot : (s12 ? s12 : s30);

  TCanvas* c = new TCanvas(Form("c_%s", hname), title, 800, 650);
  gPad->SetTicks(1,1);
  gPad->SetLogy(logy ? 1 : 0);

  base->SetTitle(title);
  base->SetMaximum(ymax);
  if (logy) base->SetMinimum(ymin);

  // Draw a frame (using base histogram)
  base->Draw("HIST");

  // Draw the stack on top (filled)
  if (haveBkg) st->Draw("HIST SAME");

  // Overlay signals
  if (s12) s12->Draw("HIST SAME");
  if (s30) s30->Draw("HIST SAME");

  // Legend: signals + ALL backgrounds (components)
  TLegend* leg = new TLegend(0.55, 0.52, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  if (s12) leg->AddEntry(s12, "Signal ggH(m12)", "l");
  if (s30) leg->AddEntry(s30, "Signal ggH(m30)", "l");

  for (size_t i=0; i<bkgComps.size(); ++i) {
    leg->AddEntry(bkgComps[i], bkgCompLabels[i].c_str(), "f");
  }

  leg->Draw();

  SaveCanvasPDF(c, outdir, (std::string(hname) + "_overlay").c_str());

  // cleanup
  delete leg;
  delete c;
  delete s12;
  delete s30;
  delete bkgTot;
  delete st;
  for (auto* h : bkgComps) delete h;
}

} // namespace ggHOverlay

// ---------- main ----------
void DrawggHOverlay(
  const char* sigFile12 = "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-12_2024.root",
  const char* sigFile30 = "ggHAnalysis_signal_GluGluH-01J_HToAATo4B_M-30_2024.root"
)
{
  gStyle->SetOptStat(0);

  const char* outdir = "overlayHists";
  ggHOverlay::EnsureDirAndCleanPDFs(outdir);

  // signals
  TFile* fS12 = TFile::Open(sigFile12, "READ");
  if (!fS12 || fS12->IsZombie()) { std::cout << "ERROR: cannot open " << sigFile12 << "\n"; return; }
  TFile* fS30 = TFile::Open(sigFile30, "READ");
  if (!fS30 || fS30->IsZombie()) { std::cout << "ERROR: cannot open " << sigFile30 << "\n"; fS12->Close(); delete fS12; return; }

  // BACKGROUNDS: list ALL files + labels (1-to-1)
  std::vector<std::string> bkgFileNames = {
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-600to800_2024.root",
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-800to1000_2024.root",
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1000-1200_2024.root",
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1200-1500_2024.root",
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-1500-2000_2024.root",
    "ggHAnalysis_qcd_QCD-4Jets_Bin-HT-2000toinf_2024.root",
  };

  std::vector<std::string> bkgLabels = {
    "QCD HT 600-800",
    "QCD HT 800-1000",
    "QCD HT 1000-1200",
    "QCD HT 1200-1500",
    "QCD HT 1500-2000",
    "QCD HT 2000toInf",
  };

  if (bkgLabels.size() != bkgFileNames.size()) {
    std::cout << "ERROR: bkgLabels.size() must equal bkgFileNames.size()\n";
    fS12->Close(); delete fS12;
    fS30->Close(); delete fS30;
    return;
  }

  std::vector<TFile*> bkgFiles;
  bkgFiles.reserve(bkgFileNames.size());
  for (const auto& fn : bkgFileNames) {
    TFile* fb = TFile::Open(fn.c_str(), "READ");
    if (!fb || fb->IsZombie()) {
      std::cout << "ERROR: cannot open background file: " << fn << "\n";
      if (fb) delete fb;
      bkgFiles.push_back(nullptr);
      continue;
    }
    bkgFiles.push_back(fb);
  }

  // hist list (use your full list; add more if needed)
  struct HDef { const char* name; const char* title; bool logy; int rebin; };
  std::vector<HDef> hists = {
    {"h_nMu",        "N_{#mu}", false, 1},
    {"h_nEle",       "N_{e}", false, 1},
    {"h_nJet",       "N_{jets}", false, 1},
    {"h_nJet_clean", "N_{jets}^{clean}", false, 1},

    {"h_pt_e1",  "e_{1} p_{T}", true, 2},
    {"h_eta_e1", "e_{1} #eta", false, 2},
    {"h_phi_e1", "e_{1} #phi", false, 2},
    {"h_pt_e2",  "e_{2} p_{T}", true, 2},
    {"h_eta_e2", "e_{2} #eta", false, 2},
    {"h_phi_e2", "e_{2} #phi", false, 2},

    {"h_pt_mu1",  "#mu_{1} p_{T}", true, 2},
    {"h_eta_mu1", "#mu_{1} #eta", false, 2},
    {"h_phi_mu1", "#mu_{1} #phi", false, 2},
    {"h_pt_mu2",  "#mu_{2} p_{T}", true, 2},
    {"h_eta_mu2", "#mu_{2} #eta", false, 2},
    {"h_phi_mu2", "#mu_{2} #phi", false, 2},

    {"h_dphi_J1J2",        "|#Delta#phi(J1,J2)| after all cuts", false, 2},
    {"h_dR_J1J2",          "#DeltaR(J1,J2) after all cuts", false, 2},
    {"h_abs_dm_J1J2",      "|m(J1)-m(J2)| after all cuts", false, 2},
    {"h_min_dphi_MET_Jet", "min |#Delta#phi(MET,Jet)| after all cuts", false, 2},
    {"h_HT",               "HT after all cuts", true, 2},
    {"h_Nbjets_after",     "Nbjets after all cuts", false, 1},

    {"h_m_2b",   "m_{bb} after all cuts", false, 2},
    {"h_pt_2b",  "p_{T}^{bb} after all cuts", true, 2},
    {"h_eta_2b", "#eta^{bb} after all cuts", false, 2},

    {"h_MET_after",     "MET after cuts", true, 2},
    {"h_MET_phi_after", "MET #phi after cuts", false, 2},

    {"h_pt_bjet1",  "bjet1 p_{T}", true, 2},
    {"h_eta_bjet1", "bjet1 #eta", false, 2},
    {"h_m_bjet1",   "bjet1 mass", false, 2},
    {"h_pt_bjet2",  "bjet2 p_{T}", true, 2},
    {"h_eta_bjet2", "bjet2 #eta", false, 2},
    {"h_m_bjet2",   "bjet2 mass", false, 2},

    {"h_dphi_MET_bb", "|#Delta#phi(MET,bb)|", false, 2},
    {"h_dphi_MET_J1", "|#Delta#phi(MET,J1)|", false, 2},
  };

  for (const auto& hd : hists) {
    ggHOverlay::DrawOne(outdir, hd.name, hd.title, hd.logy, hd.rebin,
                        fS12, fS30, bkgFiles, bkgLabels);
  }

  // close files
  fS12->Close(); delete fS12;
  fS30->Close(); delete fS30;
  for (auto* fb : bkgFiles) { if (fb) { fb->Close(); delete fb; } }

  std::cout << "Done. PDFs saved under ./" << outdir << "\n";
}

#endif
