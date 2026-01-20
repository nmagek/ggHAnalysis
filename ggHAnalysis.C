#define ggHAnalysis_cxx
#include "ggHAnalysis.h"

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

// PDG codes for the signal
const int PDG_H     = 25;
const int PDG_A     = 36;
const int PDG_B     = 5;
const int PDG_BBAR  = -5;
const int PDG_GLUON = 21;

// ============================================================================
// Main analysis loop
// ============================================================================
void ggHAnalysis::Loop(double sigma_pb, double lumi_fb, const char* outFileName)
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntries();

   const double L_INT_FB        = lumi_fb;     // [fb^-1]
   const double SIGMA_PB        = sigma_pb;    // [pb]
   const double PB_FB_TO_EVENTS = 1.0e3;       // pb*fb^-1 -> events

   const double Nexp    = SIGMA_PB * L_INT_FB * PB_FB_TO_EVENTS;
   const Long64_t Nstat = nentries;

   double wgt = 1.0;
   if (Nstat > 0) wgt = Nexp / static_cast<double>(Nstat);

   std::cout << "\n=== Normalization info ===\n";
   std::cout << "  L_int [fb^-1]   = " << L_INT_FB << "\n";
   std::cout << "  sigma [pb]      = " << SIGMA_PB << "\n";
   std::cout << "  N_exp           = " << Nexp << "\n";
   std::cout << "  N_stat          = " << Nstat << "\n";
   std::cout << "  w = N_exp/N_stat= " << wgt  << "\n";
   std::cout << "==========================\n\n";


   // Higgs window cut
   const float M_HIGGS      = 125.0f;
   const float DM_HIGGS_MAX = 50.0f;

   // finite check helper (for NaN protection)
   auto isFinite = [](float x) -> bool { return std::isfinite(x); };

   // ===================================================
   // GEN analysis switch (based on file name prefix "GluGluH")
   // ===================================================
   bool doGen = false;
   {
      std::string baseName = "";
      bool isGluGluH = false;

      TFile *curFile = fChain->GetCurrentFile();
      if (curFile) {
         baseName = gSystem->BaseName(curFile->GetName());
         if (baseName.compare(0, 7, "GluGluH") == 0) isGluGluH = true;
      }
      doGen = isGluGluH;

      std::cout << (doGen ? "GEN analysis ENABLED" : "GEN analysis DISABLED")
                << " (file=" << baseName << ")\n";
   }

   // =======================
   // GEN histograms (unweighted)
   // =======================
   TH1F *h_pt_H  = new TH1F("h_pt_H",  "p_{T}(H);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_H = new TH1F("h_eta_H", "#eta(H);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_H = new TH1F("h_phi_H", "#phi(H);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_H   = new TH1F("h_m_H",   "m(H);m_{H} [GeV];Entries",     100, 0, 200);

   TH1F *h_pt_A1  = new TH1F("h_pt_A1",  "p_{T}(A1);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_A1 = new TH1F("h_eta_A1", "#eta(A1);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_A1 = new TH1F("h_phi_A1", "#phi(A1);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_A1   = new TH1F("h_m_A1",   "m(A1);m_{A1} [GeV];Entries",     100, 0, 200);

   TH1F *h_pt_A2  = new TH1F("h_pt_A2",  "p_{T}(A2);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_A2 = new TH1F("h_eta_A2", "#eta(A2);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_A2 = new TH1F("h_phi_A2", "#phi(A2);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_A2   = new TH1F("h_m_A2",   "m(A2);m_{A2} [GeV];Entries",     100, 0, 200);

   TH1F *h_dR_AA     = new TH1F("h_dR_AA",     "#DeltaR(A1,A2);#DeltaR;Entries", 60, 0, 6);
   TH1F *h_dR_bb_A1  = new TH1F("h_dR_bb_A1",  "#DeltaR(b,b) from A1;#DeltaR;Entries", 60, 0, 6);
   TH1F *h_dR_bb_A2  = new TH1F("h_dR_bb_A2",  "#DeltaR(b,b) from A2;#DeltaR;Entries", 60, 0, 6);

   TH1F *h_pt_b1  = new TH1F("h_pt_b1",  "p_{T}(b1);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_b1 = new TH1F("h_eta_b1", "#eta(b1);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_b1 = new TH1F("h_phi_b1", "#phi(b1);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_b1   = new TH1F("h_m_b1",   "m(b1);m [GeV];Entries",         100, 0, 50);

   TH1F *h_pt_b2  = new TH1F("h_pt_b2",  "p_{T}(b2);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_b2 = new TH1F("h_eta_b2", "#eta(b2);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_b2 = new TH1F("h_phi_b2", "#phi(b2);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_b2   = new TH1F("h_m_b2",   "m(b2);m [GeV];Entries",         100, 0, 50);

   TH1F *h_pt_b3  = new TH1F("h_pt_b3",  "p_{T}(b3);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_b3 = new TH1F("h_eta_b3", "#eta(b3);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_b3 = new TH1F("h_phi_b3", "#phi(b3);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_b3   = new TH1F("h_m_b3",   "m(b3);m [GeV];Entries",         100, 0, 50);

   TH1F *h_pt_b4  = new TH1F("h_pt_b4",  "p_{T}(b4);p_{T} [GeV];Entries", 100, 0, 500);
   TH1F *h_eta_b4 = new TH1F("h_eta_b4", "#eta(b4);#eta;Entries",         60, -5, 5);
   TH1F *h_phi_b4 = new TH1F("h_phi_b4", "#phi(b4);#phi;Entries",         64, -3.2, 3.2);
   TH1F *h_m_b4   = new TH1F("h_m_b4",   "m(b4);m [GeV];Entries",         100, 0, 50);

   // =======================
   // RECO BEFORE selection (weighted)
   // =======================
   TH1I *h_nEle = new TH1I("h_nEle", "Electron multiplicity;N_{e};Events", 10, 0, 10);
   TH1I *h_nMu  = new TH1I("h_nMu",  "Muon multiplicity;N_{#mu};Events",   10, 0, 10);
   TH1I *h_nJet = new TH1I("h_nJet", "Jet multiplicity;N_{jets};Events",   20, 0, 20);
   TH1I *h_nJet_clean = new TH1I("h_nJet_clean",
                                 "Jet multiplicity after jet–lepton cleaning;N_{jets}^{clean};Events",
                                 20, 0, 20);

   TH1F *h_pt_e1   = new TH1F("h_pt_e1",   "p_{T} of e_{1};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_e1  = new TH1F("h_eta_e1",  "#eta of e_{1};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_e1  = new TH1F("h_phi_e1",  "#phi of e_{1};#phi;Entries",         64, -3.2, 3.2);

   TH1F *h_pt_e2   = new TH1F("h_pt_e2",   "p_{T} of e_{2};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_e2  = new TH1F("h_eta_e2",  "#eta of e_{2};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_e2  = new TH1F("h_phi_e2",  "#phi of e_{2};#phi;Entries",         64, -3.2, 3.2);

   TH1F *h_pt_mu1  = new TH1F("h_pt_mu1",  "p_{T} of #mu_{1};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_mu1 = new TH1F("h_eta_mu1", "#eta of #mu_{1};#eta;Entries",        60, -5.0, 5.0);
   TH1F *h_phi_mu1 = new TH1F("h_phi_mu1", "#phi of #mu_{1};#phi;Entries",        64, -3.2, 3.2);

   TH1F *h_pt_mu2  = new TH1F("h_pt_mu2",  "p_{T} of #mu_{2};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_mu2 = new TH1F("h_eta_mu2", "#eta of #mu_{2};#eta;Entries",        60, -5.0, 5.0);
   TH1F *h_phi_mu2 = new TH1F("h_phi_mu2", "#phi of #mu_{2};#phi;Entries",        64, -3.2, 3.2);

   TH1F *h_pt_J1   = new TH1F("h_pt_J1",  "p_{T}(J1);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_J1  = new TH1F("h_eta_J1", "#eta(J1);#eta;Events",         60, -5, 5);
   TH1F *h_phi_J1  = new TH1F("h_phi_J1", "#phi(J1);#phi;Events",         64, -3.2, 3.2);

   TH1F *h_pt_J2   = new TH1F("h_pt_J2",  "p_{T}(J2);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_J2  = new TH1F("h_eta_J2", "#eta(J2);#eta;Events",         60, -5, 5);
   TH1F *h_phi_J2  = new TH1F("h_phi_J2", "#phi(J2);#phi;Events",         64, -3.2, 3.2);

   TH1F *h_pt_J3   = new TH1F("h_pt_J3",  "p_{T}(J3);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_J3  = new TH1F("h_eta_J3", "#eta(J3);#eta;Events",         60, -5, 5);
   TH1F *h_phi_J3  = new TH1F("h_phi_J3", "#phi(J3);#phi;Events",         64, -3.2, 3.2);

   TH1F *h_pt_J4   = new TH1F("h_pt_J4",  "p_{T}(J4);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_J4  = new TH1F("h_eta_J4", "#eta(J4);#eta;Events",         60, -5, 5);
   TH1F *h_phi_J4  = new TH1F("h_phi_J4", "#phi(J4);#phi;Events",         64, -3.2, 3.2);

   TH1F *h_MET     = new TH1F("h_MET",     "PuppiMET;MET [GeV];Events",   100, 0, 500);
   TH1F *h_MET_phi = new TH1F("h_MET_phi", "PuppiMET #phi;#phi;Events",    64, -3.2, 3.2);

   // =======================
   // RECO AFTER all cuts (weighted)
   // =======================
   TH1F *h_dphi_J1J2 = new TH1F("h_dphi_J1J2", "|#Delta#phi(J1,J2)|;|#Delta#phi|;Events", 64, 0, TMath::Pi());
   TH1F *h_dR_J1J2   = new TH1F("h_dR_J1J2",   "#DeltaR(J1,J2);#DeltaR;Events",           60, 0, 6.0);

   TH1F *h_abs_dm_J1J2 = new TH1F("h_abs_dm_J1J2",
                                  "|m(J1)-m(J2)| (after all cuts);|#Delta m_{jj}| [GeV];Events",
                                  100, 0.0, 200.0);

   TH1F *h_m_2b   = new TH1F("h_m_2b",  "m_{bb} (after all cuts);m_{bb} [GeV];Events", 100, 0, 500);
   TH1F *h_pt_2b  = new TH1F("h_pt_2b", "p_{T}^{bb} (after all cuts);p_{T}^{bb} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_2b = new TH1F("h_eta_2b","#eta^{bb} (after all cuts);#eta^{bb};Events", 60, -5, 5);

   TH1F *h_MET_after     = new TH1F("h_MET_after",     "PuppiMET (after all cuts);MET [GeV];Events", 100, 0, 500);
   TH1F *h_MET_phi_after = new TH1F("h_MET_phi_after", "PuppiMET #phi (after all cuts);#phi;Events",  64, -3.2, 3.2);

   TH1F *h_pt_bjet1  = new TH1F("h_pt_bjet1",  "p_{T}(bjet1) (after all cuts);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_bjet1 = new TH1F("h_eta_bjet1", "#eta(bjet1) (after all cuts);#eta;Events",         60, -5, 5);
   TH1F *h_m_bjet1   = new TH1F("h_m_bjet1",   "m(bjet1) (after all cuts);m [GeV];Events",         100, 0, 200);

   TH1F *h_pt_bjet2  = new TH1F("h_pt_bjet2",  "p_{T}(bjet2) (after all cuts);p_{T} [GeV];Events", 100, 0, 1000);
   TH1F *h_eta_bjet2 = new TH1F("h_eta_bjet2", "#eta(bjet2) (after all cuts);#eta;Events",         60, -5, 5);
   TH1F *h_m_bjet2   = new TH1F("h_m_bjet2",   "m(bjet2) (after all cuts);m [GeV];Events",         100, 0, 200);

   TH1F *h_Nbjets_after = new TH1F("h_Nbjets_after",
                                   "b-tagged jet multiplicity (after all cuts);N_{b\\text{-jets}};Events",
                                   10, 0, 10);

   TH1F *h_HT = new TH1F("h_HT",
                         "H_{T} = #Sigma p_{T}^{jets} (after all cuts);H_{T} [GeV];Events",
                         100, 0.0, 2000.0);

   TH1F *h_dphi_MET_bb = new TH1F("h_dphi_MET_bb",
                                  "|#Delta#phi(MET,bb)| (after all cuts);|#Delta#phi(MET,bb)|;Events",
                                  64, 0.0, TMath::Pi());
   TH1F *h_dphi_MET_J1 = new TH1F("h_dphi_MET_J1",
                                  "|#Delta#phi(MET,J_{1})| (after all cuts);|#Delta#phi(MET,J_{1})|;Events",
                                  64, 0.0, TMath::Pi());

   TH1F *h_dR_J1_e1  = new TH1F("h_dR_J1_e1",
                                "#DeltaR(J_{1}, e_{1});#DeltaR(J_{1},e_{1});Events",
                                60, 0.0, 6.0);
   TH1F *h_dR_J1_mu1 = new TH1F("h_dR_J1_mu1",
                                "#DeltaR(J_{1}, #mu_{1});#DeltaR(J_{1},#mu_{1});Events",
                                60, 0.0, 6.0);

   TH1F *h_dR_J1_e1_clean  = new TH1F("h_dR_J1_e1_clean",
                                      "#DeltaR(J_{1}^{clean}, e_{1});#DeltaR(J_{1}^{clean},e_{1});Events",
                                      60, 0.0, 6.0);
   TH1F *h_dR_J1_mu1_clean = new TH1F("h_dR_J1_mu1_clean",
                                      "#DeltaR(J_{1}^{clean}, #mu_{1});#DeltaR(J_{1}^{clean},#mu_{1});Events",
                                      60, 0.0, 6.0);

   TH1F *h_min_dphi_MET_Jet = new TH1F("h_min_dphi_MET_Jet",
                                       "min |#Delta#phi(MET,Jet)| (after all cuts);min |#Delta#phi|;Events",
                                       64, 0.0, TMath::Pi());

   const float PT_CUT  = 20.0;
   const float ETA_CUT = 2.5;

   struct PtComparator {
      const float* pt;
      PtComparator(const float* p) : pt(p) {}
      bool operator()(int a, int b) const { return pt[a] > pt[b]; }
   };

   // Cutflow counters (unweighted)
   Long64_t nRawEvents = nentries;
   Long64_t nPassCut1  = 0; // veto loose leptons
   Long64_t nPassCut2  = 0; // N(clean jets) >= 2
   Long64_t nPassCut3  = 0; // pT(J1) > 100
   Long64_t nPassCut4  = 0; // 2 b-tag jets
   Long64_t nPassCut5  = 0; // MET < 140
   Long64_t nPassCut6  = 0; // |m_bb - 125| < 50

   // b-tag helper
   auto isBTagged = [&](int idx) -> bool {
      return (Jet_btagUParTAK4probbb[idx] > 0.38);
   };

   // =======================
   // Output file (must exist BEFORE creating the output TTree)
   // =======================
   TFile *f = new TFile(outFileName, "RECREATE");
   if (!f || f->IsZombie()) {
     std::cerr << "ERROR: could not create output file ggHAnalysis_plots.root\n";
     return;
   }
   f->cd();

   
   // =======================
   // Output tree: RECO variables after all cuts (one entry per selected event)
   // =======================
   TTree *tRecoAfterCuts = new TTree("RecoAfterCuts", "Reco variables after all cuts");
   tRecoAfterCuts->SetDirectory(f);

   // Variables to write
   float out_wgt = 0.0f;

   float out_m_2b = 0.0f, out_pt_2b = 0.0f, out_eta_2b = 0.0f;
   float out_met = 0.0f, out_met_phi = 0.0f;

   float out_dphi_j1j2 = 0.0f, out_dR_j1j2 = 0.0f, out_abs_dm_j1j2 = 0.0f;

   float out_pt_bjet1 = 0.0f, out_eta_bjet1 = 0.0f, out_m_bjet1 = 0.0f;
   float out_pt_bjet2 = 0.0f, out_eta_bjet2 = 0.0f, out_m_bjet2 = 0.0f;

   float out_nbjets = 0.0f;
   float out_ht = 0.0f;

   float out_dphi_met_bb = 0.0f, out_dphi_met_j1 = 0.0f;
   float out_min_dphi_met_jet = 0.0f;

   // Optional: store leading cleaned jets too (often useful)
   float out_pt_j1 = 0.0f, out_eta_j1 = 0.0f, out_phi_j1 = 0.0f, out_m_j1 = 0.0f;
   float out_pt_j2 = 0.0f, out_eta_j2 = 0.0f, out_phi_j2 = 0.0f, out_m_j2 = 0.0f;

   // Branches
   tRecoAfterCuts->Branch("wgt", &out_wgt, "wgt/F");

   tRecoAfterCuts->Branch("m_2b",  &out_m_2b,  "m_2b/F");
   tRecoAfterCuts->Branch("pt_2b", &out_pt_2b, "pt_2b/F");
   tRecoAfterCuts->Branch("eta_2b",&out_eta_2b,"eta_2b/F");

   tRecoAfterCuts->Branch("met",     &out_met,     "met/F");
   tRecoAfterCuts->Branch("met_phi", &out_met_phi, "met_phi/F");

   tRecoAfterCuts->Branch("dphi_j1j2",   &out_dphi_j1j2,   "dphi_j1j2/F");
   tRecoAfterCuts->Branch("dR_j1j2",     &out_dR_j1j2,     "dR_j1j2/F");
   tRecoAfterCuts->Branch("abs_dm_j1j2", &out_abs_dm_j1j2, "abs_dm_j1j2/F");

   tRecoAfterCuts->Branch("pt_bjet1",  &out_pt_bjet1,  "pt_bjet1/F");
   tRecoAfterCuts->Branch("eta_bjet1", &out_eta_bjet1, "eta_bjet1/F");
   tRecoAfterCuts->Branch("m_bjet1",   &out_m_bjet1,   "m_bjet1/F");

   tRecoAfterCuts->Branch("pt_bjet2",  &out_pt_bjet2,  "pt_bjet2/F");
   tRecoAfterCuts->Branch("eta_bjet2", &out_eta_bjet2, "eta_bjet2/F");
   tRecoAfterCuts->Branch("m_bjet2",   &out_m_bjet2,   "m_bjet2/F");

   tRecoAfterCuts->Branch("nbjets", &out_nbjets, "nbjets/F");
   tRecoAfterCuts->Branch("ht",     &out_ht,     "ht/F");

   tRecoAfterCuts->Branch("dphi_met_bb", &out_dphi_met_bb, "dphi_met_bb/F");
   tRecoAfterCuts->Branch("dphi_met_j1", &out_dphi_met_j1, "dphi_met_j1/F");
   tRecoAfterCuts->Branch("min_dphi_met_jet", &out_min_dphi_met_jet, "min_dphi_met_jet/F");

   // Optional jet branches
   tRecoAfterCuts->Branch("pt_j1",  &out_pt_j1,  "pt_j1/F");
   tRecoAfterCuts->Branch("eta_j1", &out_eta_j1, "eta_j1/F");
   tRecoAfterCuts->Branch("phi_j1", &out_phi_j1, "phi_j1/F");
   tRecoAfterCuts->Branch("m_j1",   &out_m_j1,   "m_j1/F");

   tRecoAfterCuts->Branch("pt_j2",  &out_pt_j2,  "pt_j2/F");
   tRecoAfterCuts->Branch("eta_j2", &out_eta_j2, "eta_j2/F");
   tRecoAfterCuts->Branch("phi_j2", &out_phi_j2, "phi_j2/F");
   tRecoAfterCuts->Branch("m_j2",   &out_m_j2,   "m_j2/F");

   
   // -----------------------------
   // Event loop
   // -----------------------------
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);

      // Guard against pathological events with NaN MET quantities
      if (!isFinite(PuppiMET_pt) || !isFinite(PuppiMET_phi)) continue;

      auto deltaPhiAbs = [&](float phi1, float phi2) -> float {
    if (!std::isfinite(phi1) || !std::isfinite(phi2)) return 999.f;
    float dphi = phi1 - phi2;
    while (dphi >  (float)TMath::Pi())  dphi -= 2.0f*(float)TMath::Pi();
    while (dphi <= -(float)TMath::Pi()) dphi += 2.0f*(float)TMath::Pi();
    return std::fabs(dphi);
      };


      // GEN block (if enabled)
      if (doGen) {
                 std::map<int, std::vector<int> > higgsToAs;

         for (int i = 0; i < nGenPart; ++i) {
            if (GenPart_pdgId[i] != PDG_A) continue;

            int momIdx = GenPart_genPartIdxMother[i];
            if (momIdx < 0 || momIdx >= nGenPart) continue;
            if (GenPart_pdgId[momIdx] != PDG_H) continue;

            higgsToAs[momIdx].push_back(i);
         }

         int idxH = -1;
         std::vector<int> idxA;

         for (std::map<int, std::vector<int> >::iterator it = higgsToAs.begin();
              it != higgsToAs.end(); ++it) {
            if (it->second.size() >= 2) {
               idxH = it->first;
               idxA = it->second;
               break;
            }
         }

         if (idxH >= 0 && idxA.size() >= 2) {

            int idxA1 = idxA[0];
            int idxA2 = idxA[1];

            h_pt_H->Fill(GenPart_pt[idxH]);
            h_eta_H->Fill(GenPart_eta[idxH]);
            h_phi_H->Fill(GenPart_phi[idxH]);
            h_m_H->Fill(GenPart_mass[idxH]);

            h_pt_A1->Fill(GenPart_pt[idxA1]);
            h_eta_A1->Fill(GenPart_eta[idxA1]);
            h_phi_A1->Fill(GenPart_phi[idxA1]);
            h_m_A1->Fill(GenPart_mass[idxA1]);

            h_pt_A2->Fill(GenPart_pt[idxA2]);
            h_eta_A2->Fill(GenPart_eta[idxA2]);
            h_phi_A2->Fill(GenPart_phi[idxA2]);
            h_m_A2->Fill(GenPart_mass[idxA2]);

            {
               float deta_A = GenPart_eta[idxA1] - GenPart_eta[idxA2];
               float dphi_A = GenPart_phi[idxA1] - GenPart_phi[idxA2];

               while (dphi_A >  TMath::Pi())  dphi_A -= 2.0f*TMath::Pi();
	       while (dphi_A <= -TMath::Pi()) dphi_A += 2.0f*TMath::Pi();


               float dR_AA = std::sqrt(deta_A*deta_A + dphi_A*dphi_A);
               h_dR_AA->Fill(dR_AA);
            }

            // ΔR(b,b) from each A
            for (int which = 0; which < 2; ++which) {
               int idxA_X = (which == 0 ? idxA1 : idxA2);
               TH1F* hist = (which == 0 ? h_dR_bb_A1 : h_dR_bb_A2);

               std::vector<int> b_from_A;

               for (int i = 0; i < nGenPart; ++i) {
                   int pdg = GenPart_pdgId[i];
                  if (std::abs(pdg) != PDG_B) continue;

                  int mom = GenPart_genPartIdxMother[i];
                  if (mom != idxA_X) continue;

                  b_from_A.push_back(i);
               }

               if (b_from_A.size() < 2) continue;

               PtComparator cmpB_local(GenPart_pt);
               std::sort(b_from_A.begin(), b_from_A.end(), cmpB_local);

               int ib1 = b_from_A[0];
               int ib2 = b_from_A[1];

               float deta_b = GenPart_eta[ib1] - GenPart_eta[ib2];
               float dphi_b = GenPart_phi[ib1] - GenPart_phi[ib2];

	      while (dphi_b >  TMath::Pi())  dphi_b -= 2.0f*TMath::Pi();
	      while (dphi_b <= -TMath::Pi()) dphi_b += 2.0f*TMath::Pi();

               float dR_bb = std::sqrt(deta_b*deta_b + dphi_b*dphi_b);
               hist->Fill(dR_bb);
            }

            // b-quarks from all A's sorted by pT
            std::vector<int> b_all_idx;

            for (int i = 0; i < nGenPart; ++i) {
               int pdg = GenPart_pdgId[i];
               if (std::abs(pdg) != PDG_B) continue;

               int mom = GenPart_genPartIdxMother[i];

               bool fromThisHAs = false;
               for (size_t k = 0; k < idxA.size(); ++k) {
                  if (mom == idxA[k]) {
                     fromThisHAs = true;
                     break;
                  }
               }
               if (!fromThisHAs) continue;

               b_all_idx.push_back(i);
            }

            PtComparator cmpB(GenPart_pt);
            std::sort(b_all_idx.begin(), b_all_idx.end(), cmpB);

            if (b_all_idx.size() >= 1) {
               int i1 = b_all_idx[0];
               h_pt_b1->Fill(GenPart_pt[i1]);
               h_eta_b1->Fill(GenPart_eta[i1]);
               h_phi_b1->Fill(GenPart_phi[i1]);
               h_m_b1->Fill(GenPart_mass[i1]);
            }
            if (b_all_idx.size() >= 2) {
               int i2 = b_all_idx[1];
               h_pt_b2->Fill(GenPart_pt[i2]);
               h_eta_b2->Fill(GenPart_eta[i2]);
               h_phi_b2->Fill(GenPart_phi[i2]);
               h_m_b2->Fill(GenPart_mass[i2]);
            }
            if (b_all_idx.size() >= 3) {
               int i3 = b_all_idx[2];
               h_pt_b3->Fill(GenPart_pt[i3]);
               h_eta_b3->Fill(GenPart_eta[i3]);
               h_phi_b3->Fill(GenPart_phi[i3]);
               h_m_b3->Fill(GenPart_mass[i3]);
            }
            if (b_all_idx.size() >= 4) {
               int i4 = b_all_idx[3];
               h_pt_b4->Fill(GenPart_pt[i4]);
               h_eta_b4->Fill(GenPart_eta[i4]);
               h_phi_b4->Fill(GenPart_phi[i4]);
               h_m_b4->Fill(GenPart_mass[i4]);
            }
	    
      }
 }
      // Build loose lepton counts (Cut1) + tight collections for "before selection" plots
      int nLooseMu  = 0;
      int nLooseEle = 0;

      std::vector<int> mu_idx_pass;
      std::vector<int> ele_idx_pass;
      std::vector<int> jet_idx_pass;

      // Muons
      for (int i = 0; i < nMuon; ++i) {
         float pt  = Muon_pt[i];
         float eta = Muon_eta[i];

         bool passLooseMuKin = (pt > 10.0 && std::fabs(eta) < 2.4);
         bool passLooseMuId  = (Muon_looseId[i] != 0);
         if (passLooseMuKin && passLooseMuId) ++nLooseMu;

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Muon_tightId[i] != 0);
         bool passIso = (Muon_pfRelIso04_all[i] < 0.15);

         if (!passKin) continue;
         if (!passId)  continue;
         if (!passIso) continue;

         mu_idx_pass.push_back(i);
      }

      // Electrons
      for (int i = 0; i < nElectron; ++i) {
         float pt  = Electron_pt[i];
         float eta = Electron_eta[i];

         bool passLooseEleKin = (pt > 10.0 && std::fabs(eta) < 2.5);
         bool passLooseEleId  = (Electron_cutBased[i] >= 1);
         if (passLooseEleKin && passLooseEleId) ++nLooseEle;

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Electron_cutBased[i] >= 3);
         bool passIso = (Electron_pfRelIso03_all[i] < 0.15);

         if (!passKin) continue;
         if (!passId)  continue;
         if (!passIso) continue;

         ele_idx_pass.push_back(i);
      }

      // Jets
      for (int i = 0; i < nJet; ++i) {
         float pt  = Jet_pt[i];
         float eta = Jet_eta[i];

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Jet_passJetIdTightLepVeto[i] != 0);

         if (!passKin) continue;
         if (!passId)  continue;

         // Optional: guard obviously bad jet phi values at source
         if (!isFinite(Jet_phi[i]) || !isFinite(Jet_mass[i]) || !isFinite(Jet_pt[i]) || !isFinite(Jet_eta[i])) continue;

         jet_idx_pass.push_back(i);
      }

      // BEFORE-selection histograms (weighted)
      h_nMu->Fill((int)mu_idx_pass.size(), wgt);
      h_nEle->Fill((int)ele_idx_pass.size(), wgt);
      h_nJet->Fill((int)jet_idx_pass.size(), wgt);

      // Sort by pT
      if (!mu_idx_pass.empty()) {
         PtComparator cmpMu(Muon_pt);
         std::sort(mu_idx_pass.begin(), mu_idx_pass.end(), cmpMu);
      }
      if (!ele_idx_pass.empty()) {
         PtComparator cmpEle(Electron_pt);
         std::sort(ele_idx_pass.begin(), ele_idx_pass.end(), cmpEle);
      }
      if (!jet_idx_pass.empty()) {
         PtComparator cmpJet(Jet_pt);
         std::sort(jet_idx_pass.begin(), jet_idx_pass.end(), cmpJet);
      }

      // ----------------------------------------------------
      // Δphi(J1,J2) and ΔR(J1,J2) BEFORE cleaning (raw jets)
      // ----------------------------------------------------
      if (jet_idx_pass.size() >= 2) {
	int j1_raw = jet_idx_pass[0];
	int j2_raw = jet_idx_pass[1];

	// Guard against NaNs
	if (isFinite(Jet_pt[j1_raw]) && isFinite(Jet_eta[j1_raw]) &&
	    isFinite(Jet_phi[j1_raw]) && isFinite(Jet_mass[j1_raw]) &&
	    isFinite(Jet_pt[j2_raw]) && isFinite(Jet_eta[j2_raw]) &&
	    isFinite(Jet_phi[j2_raw]) && isFinite(Jet_mass[j2_raw])) {

	  TLorentzVector J1raw, J2raw;
	  J1raw.SetPtEtaPhiM(Jet_pt[j1_raw], Jet_eta[j1_raw], Jet_phi[j1_raw], Jet_mass[j1_raw]);
	  J2raw.SetPtEtaPhiM(Jet_pt[j2_raw], Jet_eta[j2_raw], Jet_phi[j2_raw], Jet_mass[j2_raw]);

	  float dphi_raw = deltaPhiAbs((float)J1raw.Phi(), (float)J2raw.Phi());
	  float dR_raw   = J1raw.DeltaR(J2raw);

	  h_dphi_J1J2->Fill(dphi_raw, wgt);
	  h_dR_J1J2->Fill(dR_raw, wgt);
	}
      }

      
      // Jet–lepton cleaning (ΔR < 0.4 veto w.r.t. tight leptons)
      std::vector<int> jet_idx_clean;
      jet_idx_clean.reserve(jet_idx_pass.size());

      for (int j : jet_idx_pass) {
         TLorentzVector J;
         J.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);

         bool overlaps = false;

         for (int e : ele_idx_pass) {
            TLorentzVector E;
            E.SetPtEtaPhiM(Electron_pt[e], Electron_eta[e], Electron_phi[e], 0.0);
            if (J.DeltaR(E) < 0.4) { overlaps = true; break; }
         }

         if (!overlaps) {
            for (int m : mu_idx_pass) {
               TLorentzVector M;
               M.SetPtEtaPhiM(Muon_pt[m], Muon_eta[m], Muon_phi[m], Muon_mass[m]);
               if (J.DeltaR(M) < 0.4) { overlaps = true; break; }
            }
         }

         if (!overlaps) jet_idx_clean.push_back(j);
      }
          
      // Fill cleaned jet multiplicity
      h_nJet_clean->Fill((int)jet_idx_clean.size(), wgt);

      // Fill "before selection" kinematics
      if (mu_idx_pass.size() >= 1) {
         int i1 = mu_idx_pass[0];
         h_pt_mu1->Fill(Muon_pt[i1], wgt);
         h_eta_mu1->Fill(Muon_eta[i1], wgt);
         h_phi_mu1->Fill(Muon_phi[i1], wgt);
      }
      if (mu_idx_pass.size() >= 2) {
         int i2 = mu_idx_pass[1];
         h_pt_mu2->Fill(Muon_pt[i2], wgt);
         h_eta_mu2->Fill(Muon_eta[i2], wgt);
         h_phi_mu2->Fill(Muon_phi[i2], wgt);
      }

      if (ele_idx_pass.size() >= 1) {
         int i1 = ele_idx_pass[0];
         h_pt_e1->Fill(Electron_pt[i1], wgt);
         h_eta_e1->Fill(Electron_eta[i1], wgt);
         h_phi_e1->Fill(Electron_phi[i1], wgt);
      }
      if (ele_idx_pass.size() >= 2) {
         int i2 = ele_idx_pass[1];
         h_pt_e2->Fill(Electron_pt[i2], wgt);
         h_eta_e2->Fill(Electron_eta[i2], wgt);
         h_phi_e2->Fill(Electron_phi[i2], wgt);
      }

      if (jet_idx_pass.size() >= 1) {
         int j1b = jet_idx_pass[0];
         h_pt_J1->Fill(Jet_pt[j1b], wgt);
         h_eta_J1->Fill(Jet_eta[j1b], wgt);
         h_phi_J1->Fill(Jet_phi[j1b], wgt);
      }
      if (jet_idx_pass.size() >= 2) {
         int j2b = jet_idx_pass[1];
         h_pt_J2->Fill(Jet_pt[j2b], wgt);
         h_eta_J2->Fill(Jet_eta[j2b], wgt);
         h_phi_J2->Fill(Jet_phi[j2b], wgt);
      }
      if (jet_idx_pass.size() >= 3) {
         int j3b = jet_idx_pass[2];
         h_pt_J3->Fill(Jet_pt[j3b], wgt);
         h_eta_J3->Fill(Jet_eta[j3b], wgt);
         h_phi_J3->Fill(Jet_phi[j3b], wgt);
      }
      if (jet_idx_pass.size() >= 4) {
         int j4b = jet_idx_pass[3];
         h_pt_J4->Fill(Jet_pt[j4b], wgt);
         h_eta_J4->Fill(Jet_eta[j4b], wgt);
         h_phi_J4->Fill(Jet_phi[j4b], wgt);
      }

      h_MET->Fill(PuppiMET_pt, wgt);
      h_MET_phi->Fill(PuppiMET_phi, wgt);
     
      
      
      // Cut1: veto loose leptons
      if (nLooseMu > 0 || nLooseEle > 0) continue;
      ++nPassCut1;

      // Cut2: at least 2 CLEANED jets
      if (jet_idx_clean.size() < 2) continue;
      ++nPassCut2;

      // Cut3: pT(J1) > 100 (cleaned leading jet)
      int j1 = jet_idx_clean[0];
      int j2 = jet_idx_clean[1];
      if (Jet_pt[j1] <= 100.0) continue;
      ++nPassCut3;

      // Cut4: require >=2 b-tagged jets among CLEANED jets
      std::vector<int> bjet_idx;
      bjet_idx.reserve(jet_idx_clean.size());
      for (int idx : jet_idx_clean) {
         if (isBTagged(idx)) bjet_idx.push_back(idx);
      }
      if (bjet_idx.size() < 2) continue;
      ++nPassCut4;

      // Cut5: MET < 140
      if (PuppiMET_pt >= 140.0) continue;
      ++nPassCut5;

      // Sort b-jets by pT
      PtComparator cmpB(Jet_pt);
      std::sort(bjet_idx.begin(), bjet_idx.end(), cmpB);

      int b1 = bjet_idx[0];
      int b2 = bjet_idx[1];

      // Guard b-jets used in Higgs candidate
      if (!isFinite(Jet_pt[b1]) || !isFinite(Jet_eta[b1]) || !isFinite(Jet_phi[b1]) || !isFinite(Jet_mass[b1])) continue;
      if (!isFinite(Jet_pt[b2]) || !isFinite(Jet_eta[b2]) || !isFinite(Jet_phi[b2]) || !isFinite(Jet_mass[b2])) continue;

      TLorentzVector B1, B2, Hcand;
      B1.SetPtEtaPhiM(Jet_pt[b1], Jet_eta[b1], Jet_phi[b1], Jet_mass[b1]);
      B2.SetPtEtaPhiM(Jet_pt[b2], Jet_eta[b2], Jet_phi[b2], Jet_mass[b2]);
      Hcand = B1 + B2;

      const float m_2b = (float)Hcand.M();

      // Cut6: |m_2b - 125| < 50 GeV
      if (std::fabs(m_2b - M_HIGGS) >= DM_HIGGS_MAX) continue;
      ++nPassCut6;

      // Guard leading cleaned jets against NaN inputs before building TLorentzVectors
      if (!isFinite(Jet_pt[j1]) || !isFinite(Jet_eta[j1]) || !isFinite(Jet_phi[j1]) || !isFinite(Jet_mass[j1])) continue;
      if (!isFinite(Jet_pt[j2]) || !isFinite(Jet_eta[j2]) || !isFinite(Jet_phi[j2]) || !isFinite(Jet_mass[j2])) continue;

      // AFTER-CUTS histograms
      TLorentzVector J1, J2;
      J1.SetPtEtaPhiM(Jet_pt[j1], Jet_eta[j1], Jet_phi[j1], Jet_mass[j1]);
      J2.SetPtEtaPhiM(Jet_pt[j2], Jet_eta[j2], Jet_phi[j2], Jet_mass[j2]);

      float dphi = deltaPhiAbs((float)J1.Phi(), (float)J2.Phi());
      float dR   = J1.DeltaR(J2);

      float dmjj = std::fabs((float)J1.M() - (float)J2.M());
      h_abs_dm_J1J2->Fill(dmjj, wgt);

      h_m_2b->Fill(m_2b, wgt);
      h_pt_2b->Fill((float)Hcand.Pt(), wgt);
      h_eta_2b->Fill((float)Hcand.Eta(), wgt);

      h_MET_after->Fill(PuppiMET_pt, wgt);
      h_MET_phi_after->Fill(PuppiMET_phi, wgt);

      h_pt_bjet1->Fill(Jet_pt[b1], wgt);
      h_eta_bjet1->Fill(Jet_eta[b1], wgt);
      h_m_bjet1->Fill(Jet_mass[b1], wgt);

      h_pt_bjet2->Fill(Jet_pt[b2], wgt);
      h_eta_bjet2->Fill(Jet_eta[b2], wgt);
      h_m_bjet2->Fill(Jet_mass[b2], wgt);

      h_Nbjets_after->Fill((int)bjet_idx.size(), wgt);

      double HT = 0.0;
      for (int idx : jet_idx_clean) HT += Jet_pt[idx];
      h_HT->Fill(HT, wgt);

      float dphi_met_bb = deltaPhiAbs((float)PuppiMET_phi, (float)Hcand.Phi());
      float dphi_met_j1 = deltaPhiAbs((float)PuppiMET_phi, (float)J1.Phi());
      h_dphi_MET_bb->Fill(dphi_met_bb, wgt);
      h_dphi_MET_J1->Fill(dphi_met_j1, wgt);

      float minDphi = 999.f;
      for (int idx : jet_idx_clean) {
         if (!isFinite(Jet_phi[idx])) continue;
         float dph = deltaPhiAbs((float)PuppiMET_phi, (float)Jet_phi[idx]);
         if (dph < minDphi) minDphi = dph;
      }
      if (minDphi < 900.f) h_min_dphi_MET_Jet->Fill(minDphi, wgt);

      // -----------------------------
      // Fill output tree (after all cuts)
      // -----------------------------
      out_wgt = (float)wgt;

      out_m_2b  = (float)m_2b;
      out_pt_2b = (float)Hcand.Pt();
      out_eta_2b = (float)Hcand.Eta();

      out_met = (float)PuppiMET_pt;
      out_met_phi = (float)PuppiMET_phi;

      out_dphi_j1j2 = (float)dphi;
      out_dR_j1j2   = (float)dR;
      out_abs_dm_j1j2 = (float)dmjj;

      out_pt_bjet1  = (float)Jet_pt[b1];
      out_eta_bjet1 = (float)Jet_eta[b1];
      out_m_bjet1   = (float)Jet_mass[b1];

      out_pt_bjet2  = (float)Jet_pt[b2];
      out_eta_bjet2 = (float)Jet_eta[b2];
      out_m_bjet2   = (float)Jet_mass[b2];

      out_nbjets = (float)bjet_idx.size();
      out_ht = (float)HT;

      out_dphi_met_bb = (float)dphi_met_bb;
      out_dphi_met_j1 = (float)dphi_met_j1;
      out_min_dphi_met_jet = (minDphi < 900.f ? minDphi : -1.f);

      // Optional leading cleaned jets
      out_pt_j1  = (float)Jet_pt[j1];
      out_eta_j1 = (float)Jet_eta[j1];
      out_phi_j1 = (float)Jet_phi[j1];
      out_m_j1   = (float)Jet_mass[j1];

      out_pt_j2  = (float)Jet_pt[j2];
      out_eta_j2 = (float)Jet_eta[j2];
      out_phi_j2 = (float)Jet_phi[j2];
      out_m_j2   = (float)Jet_mass[j2];

      tRecoAfterCuts->Fill();
      
      if (!ele_idx_pass.empty()) {
         int e1 = ele_idx_pass[0];
         TLorentzVector E1;
         E1.SetPtEtaPhiM(Electron_pt[e1], Electron_eta[e1], Electron_phi[e1], 0.0);
         h_dR_J1_e1->Fill(J1.DeltaR(E1), wgt);
         h_dR_J1_e1_clean->Fill(J1.DeltaR(E1), wgt);
      }
      if (!mu_idx_pass.empty()) {
         int m1 = mu_idx_pass[0];
         TLorentzVector M1;
         M1.SetPtEtaPhiM(Muon_pt[m1], Muon_eta[m1], Muon_phi[m1], Muon_mass[m1]);
         h_dR_J1_mu1->Fill(J1.DeltaR(M1), wgt);
         h_dR_J1_mu1_clean->Fill(J1.DeltaR(M1), wgt);
      }
   }

   // Print cutflow
   auto eff = [&](Long64_t n) -> double {
      if (nRawEvents <= 0) return 0.0;
      return (double)n / (double)nRawEvents;
   };
   auto wNev = [&](Long64_t n) -> double { return (double)n * wgt; };

   std::cout << "---------------------------------------------------------------\n";
   std::cout << "Cutflow (unweighted counts, efficiencies, weighted yields)\n";
   std::cout << "---------------------------------------------------------------\n";

   auto printRow = [&](const char* label, Long64_t n, double eps) {
      std::cout << std::left << std::setw(22) << label
                << std::right << std::setprecision(1) << std::setw(15) << (double)n
                << std::setprecision(3) << std::setw(15) << eps
                << std::setprecision(1) << std::setw(15) << wNev(n)
                << "\n";
   };

   printRow("Raw",                 nRawEvents, 1.000);
   printRow("Cut1: veto loose ℓ",  nPassCut1,  eff(nPassCut1));
   printRow("Cut2: N_{jet}>=2",    nPassCut2,  eff(nPassCut2));
   printRow("Cut3: p_{T}(J1)>100", nPassCut3,  eff(nPassCut3));
   printRow("Cut4: 2 b-tag jets",  nPassCut4,  eff(nPassCut4));
   printRow("Cut5: MET<140 GeV",   nPassCut5,  eff(nPassCut5));
   printRow("Cut6: |m_bb-125|<50", nPassCut6,  eff(nPassCut6));

   std::cout << "---------------------------------------------------------------\n";

   //save histograms
   f->cd();

   // GEN (unweighted) -- will be empty unless you implement GEN filling
   h_pt_H->Write();   h_eta_H->Write();   h_phi_H->Write();   h_m_H->Write();
   h_pt_A1->Write();  h_eta_A1->Write();  h_phi_A1->Write();  h_m_A1->Write();
   h_pt_A2->Write();  h_eta_A2->Write();  h_phi_A2->Write();  h_m_A2->Write();
   h_dR_AA->Write();
   h_dR_bb_A1->Write();
   h_dR_bb_A2->Write();
   h_pt_b1->Write();  h_eta_b1->Write();  h_phi_b1->Write();  h_m_b1->Write();
   h_pt_b2->Write();  h_eta_b2->Write();  h_phi_b2->Write();  h_m_b2->Write();
   h_pt_b3->Write();  h_eta_b3->Write();  h_phi_b3->Write();  h_m_b3->Write();
   h_pt_b4->Write();  h_eta_b4->Write();  h_phi_b4->Write();  h_m_b4->Write();

   // RECO before selection (weighted)
   h_nEle->Write();
   h_nMu->Write();
   h_nJet->Write();
   h_nJet_clean->Write();

   h_pt_e1->Write();  h_eta_e1->Write();  h_phi_e1->Write();
   h_pt_e2->Write();  h_eta_e2->Write();  h_phi_e2->Write();

   h_pt_mu1->Write(); h_eta_mu1->Write(); h_phi_mu1->Write();
   h_pt_mu2->Write(); h_eta_mu2->Write(); h_phi_mu2->Write();

   h_pt_J1->Write();  h_eta_J1->Write();  h_phi_J1->Write();
   h_pt_J2->Write();  h_eta_J2->Write();  h_phi_J2->Write();
   h_pt_J3->Write();  h_eta_J3->Write();  h_phi_J3->Write();
   h_pt_J4->Write();  h_eta_J4->Write();  h_phi_J4->Write();

   h_MET->Write();
   h_MET_phi->Write();

   h_dR_J1_e1->Write();
   h_dR_J1_mu1->Write();
   h_dR_J1_e1_clean->Write();
   h_dR_J1_mu1_clean->Write();

   // RECO after all cuts (weighted)
   h_dphi_J1J2->Write();
   h_dR_J1J2->Write();
   h_abs_dm_J1J2->Write();
   h_min_dphi_MET_Jet->Write();

   h_m_2b->Write();
   h_pt_2b->Write();
   h_eta_2b->Write();

   h_MET_after->Write();
   h_MET_phi_after->Write();

   h_pt_bjet1->Write();
   h_eta_bjet1->Write();
   h_m_bjet1->Write();

   h_pt_bjet2->Write();
   h_eta_bjet2->Write();
   h_m_bjet2->Write();

   h_Nbjets_after->Write();
   h_HT->Write();

   h_dphi_MET_bb->Write();
   h_dphi_MET_J1->Write();

   

   tRecoAfterCuts->Write();

   f->Close();

   std::cout << "Wrote histograms to ggHAnalysis_plots.root\n";
}
