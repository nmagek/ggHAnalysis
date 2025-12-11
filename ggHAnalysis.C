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
void ggHAnalysis::Loop()
{
   if (fChain == 0) return;

   // Use GetEntries() (not GetEntriesFast) to avoid the 9e18 sentinel for TChain
   Long64_t nentries = fChain->GetEntries();

   // -----------------------------
   // Physics normalization inputs
   // -----------------------------
   const double L_INT_FB        = 108.96;  // Integrated luminosity [fb^-1]
   const double SIGMA_GGH_PB    = 48.58;   // ggH cross section [pb]
   const double PB_FB_TO_EVENTS = 1.0e3;   // pb * fb^-1 → events

   const double Nexp    = SIGMA_GGH_PB * L_INT_FB * PB_FB_TO_EVENTS;
   const Long64_t Nstat = nentries;
   
   // Normalization weight for cutflow + RECO histograms
   double wgt(1.0);
   if (Nstat > 0) {
      wgt = Nexp / static_cast<double>(Nstat);   // w = Nexp / Nstat
   }

   std::cout << "\n=== Normalization info ===\n";
   std::cout << "  L_int        = " << L_INT_FB      << " fb^-1\n";
   std::cout << "  sigma_ggH    = " << SIGMA_GGH_PB  << " pb\n";
   std::cout << "  N_exp        = " << Nexp          << " events (theory)\n";
   std::cout << "  N_stat       = " << Nstat         << " events (in ntuple)\n";
   std::cout << "  w = N_exp/N_stat = " << wgt       << "\n";
   std::cout << "==========================\n\n";

   // ===================================================
   // GEN analysis switch (based on file name prefix "GluGluH")
   // ===================================================
   bool doGen = false;  // will be decided from filename

   // =======================
   // GEN histograms (unweighted)
   // =======================
   TH1F *h_pt_H   = new TH1F("h_pt_H",   "p_{T} of Higgs (H);p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_H  = new TH1F("h_eta_H",  "#eta of Higgs (H);#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_H  = new TH1F("h_phi_H",  "#phi of Higgs (H);#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_H    = new TH1F("h_m_H",    "Mass of Higgs (H);m_{H} [GeV];Entries",  100, 0.0, 200.0);

   TH1F *h_pt_A1   = new TH1F("h_pt_A1",   "p_{T} of A_{1} (from H);p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_A1  = new TH1F("h_eta_A1",  "#eta of A_{1} (from H);#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_A1  = new TH1F("h_phi_A1",  "#phi of A_{1} (from H);#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_A1    = new TH1F("h_m_A1",    "Mass of A_{1} (from H);m_{A} [GeV];Entries",  100, 0.0, 50.0);

   TH1F *h_pt_A2   = new TH1F("h_pt_A2",   "p_{T} of A_{2} (from H);p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_A2  = new TH1F("h_eta_A2",  "#eta of A_{2} (from H);#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_A2  = new TH1F("h_phi_A2",  "#phi of A_{2} (from H);#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_A2    = new TH1F("h_m_A2",    "Mass of A_{2} (from H);m_{A} [GeV];Entries",  100, 0.0, 50.0);

   TH1F *h_dR_AA  = new TH1F("h_dR_AA",
                             "#DeltaR(A_{1},A_{2});#DeltaR(A,A);Entries",
                             60, 0.0, 6.0);

   TH1F *h_dR_bb_A1 = new TH1F("h_dR_bb_A1",
                               "#DeltaR(b,b) from A_{1};#DeltaR(b,b);Entries",
                               60, 0.0, 6.0);

   TH1F *h_dR_bb_A2 = new TH1F("h_dR_bb_A2",
                               "#DeltaR(b,b) from A_{2};#DeltaR(b,b);Entries",
                               60, 0.0, 6.0);

   TH1F *h_pt_b1   = new TH1F("h_pt_b1",   "p_{T} of b_{1};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_b1  = new TH1F("h_eta_b1",  "#eta of b_{1};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_b1  = new TH1F("h_phi_b1",  "#phi of b_{1};#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_b1    = new TH1F("h_m_b1",    "Mass of b_{1};m_{b} [GeV];Entries",  100, 0.0, 10.0);

   TH1F *h_pt_b2   = new TH1F("h_pt_b2",   "p_{T} of b_{2};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_b2  = new TH1F("h_eta_b2",  "#eta of b_{2};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_b2  = new TH1F("h_phi_b2",  "#phi of b_{2};#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_b2    = new TH1F("h_m_b2",    "Mass of b_{2};m_{b} [GeV];Entries",  100, 0.0, 10.0);

   TH1F *h_pt_b3   = new TH1F("h_pt_b3",   "p_{T} of b_{3};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_b3  = new TH1F("h_eta_b3",  "#eta of b_{3};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_b3  = new TH1F("h_phi_b3",  "#phi of b_{3};#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_b3    = new TH1F("h_m_b3",    "Mass of b_{3};m_{b} [GeV];Entries",  100, 0.0, 10.0);

   TH1F *h_pt_b4   = new TH1F("h_pt_b4",   "p_{T} of b_{4};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_b4  = new TH1F("h_eta_b4",  "#eta of b_{4};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_b4  = new TH1F("h_phi_b4",  "#phi of b_{4};#phi;Entries",         16, -3.2, 3.2);
   TH1F *h_m_b4    = new TH1F("h_m_b4",    "Mass of b_{4};m_{b} [GeV];Entries",  100, 0.0, 10.0);

   // =======================
   // RECO histograms (BEFORE selection) – WEIGHTED with wgt
   // =======================
   TH1I *h_nEle  = new TH1I("h_nEle",  "Electron multiplicity (p_{T}>20, |#eta|<2.5, ID+ISO);N_{e};Events",    10, 0, 10);
   TH1I *h_nMu   = new TH1I("h_nMu",   "Muon multiplicity (p_{T}>20, |#eta|<2.5, ID+ISO);N_{#mu};Events",      10, 0, 10);
   TH1I *h_nJet  = new TH1I("h_nJet",  "Jet multiplicity (p_{T}>20, |#eta|<2.5, tightLepVeto);N_{jets};Events",20, 0, 20);

   TH1F *h_pt_e1   = new TH1F("h_pt_e1",   "p_{T} of e_{1};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_e1  = new TH1F("h_eta_e1",  "#eta of e_{1};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_e1  = new TH1F("h_phi_e1",  "#phi of e_{1};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_pt_e2   = new TH1F("h_pt_e2",   "p_{T} of e_{2};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_e2  = new TH1F("h_eta_e2",  "#eta of e_{2};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_e2  = new TH1F("h_phi_e2",  "#phi of e_{2};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_pt_mu1  = new TH1F("h_pt_mu1",  "p_{T} of #mu_{1};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_mu1 = new TH1F("h_eta_mu1", "#eta of #mu_{1};#eta;Entries",        60, -5.0, 5.0);
   TH1F *h_phi_mu1 = new TH1F("h_phi_mu1", "#phi of #mu_{1};#phi;Entries",        16, -3.2, 3.2);

   TH1F *h_pt_mu2  = new TH1F("h_pt_mu2",  "p_{T} of #mu_{2};p_{T} [GeV];Entries", 100, 0.0, 500.0);
   TH1F *h_eta_mu2 = new TH1F("h_eta_mu2", "#eta of #mu_{2};#eta;Entries",        60, -5.0, 5.0);
   TH1F *h_phi_mu2 = new TH1F("h_phi_mu2", "#phi of #mu_{2};#phi;Entries",        16, -3.2, 3.2);

   TH1F *h_pt_J1   = new TH1F("h_pt_J1",   "p_{T} of J_{1};p_{T} [GeV];Entries", 100, 0.0, 1000.0);
   TH1F *h_eta_J1  = new TH1F("h_eta_J1",  "#eta of J_{1};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_J1  = new TH1F("h_phi_J1",  "#phi of J_{1};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_pt_J2   = new TH1F("h_pt_J2",   "p_{T} of J_{2};p_{T} [GeV];Entries", 100, 0.0, 1000.0);
   TH1F *h_eta_J2  = new TH1F("h_eta_J2",  "#eta of J_{2};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_J2  = new TH1F("h_phi_J2",  "#phi of J_{2};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_pt_J3   = new TH1F("h_pt_J3",   "p_{T} of J_{3};p_{T} [GeV];Entries", 100, 0.0, 1000.0);
   TH1F *h_eta_J3  = new TH1F("h_eta_J3",  "#eta of J_{3};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_J3  = new TH1F("h_phi_J3",  "#phi of J_{3};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_pt_J4   = new TH1F("h_pt_J4",   "p_{T} of J_{4};p_{T} [GeV];Entries", 100, 0.0, 1000.0);
   TH1F *h_eta_J4  = new TH1F("h_eta_J4",  "#eta of J_{4};#eta;Entries",         60, -5.0, 5.0);
   TH1F *h_phi_J4  = new TH1F("h_phi_J4",  "#phi of J_{4};#phi;Entries",         16, -3.2, 3.2);

   TH1F *h_MET     = new TH1F("h_MET",     "Puppi MET (before sel);p_{T}^{miss} [GeV];Events", 100, 0.0, 500.0);
   TH1F *h_MET_phi = new TH1F("h_MET_phi", "Puppi MET #phi (before sel);#phi^{miss};Events",   64, -3.2, 3.2);

   // =======================
   // RECO histograms (AFTER ALL CUTS) – WEIGHTED
   // =======================

   // |Δφ(J1,J2)| AFTER ALL CUTS (weighted)
   TH1F *h_dphi_J1J2 = new TH1F("h_dphi_J1J2",
                                "|#Delta#phi(J_{1},J_{2})| (after all cuts);|#Delta#phi|;Events",
                                64, 0.0, TMath::Pi());

   // Reconstructed 2b system (Higgs candidate) AFTER ALL CUTS (weighted)
   TH1F *h_m_2b   = new TH1F("h_m_2b",
                             "Invariant mass of two b-tagged jets (after all cuts);m_{bb} [GeV];Events",
                             100, 0.0, 500.0);

   TH1F *h_pt_2b  = new TH1F("h_pt_2b",
                             "p_{T} of two b-tagged jets (after all cuts);p_{T}^{bb} [GeV];Events",
                             100, 0.0, 1000.0);

   TH1F *h_eta_2b = new TH1F("h_eta_2b",
                             "#eta of two b-tagged jets (after all cuts);#eta^{bb};Events",
                             60, -5.0, 5.0);

   // MET after all cuts
   TH1F *h_MET_after = new TH1F("h_MET_after",
                                "Puppi MET (after all cuts);p_{T}^{miss} [GeV];Events",
                                100, 0.0, 500.0);

   // b-tagged jets (leading & subleading) after all cuts
   TH1F *h_pt_bjet1  = new TH1F("h_pt_bjet1",
                                "Leading b-tagged jet p_{T} (after all cuts);p_{T} [GeV];Events",
                                100, 0.0, 1000.0);
   TH1F *h_eta_bjet1 = new TH1F("h_eta_bjet1",
                                "Leading b-tagged jet #eta (after all cuts);#eta;Events",
                                60, -5.0, 5.0);
   TH1F *h_m_bjet1   = new TH1F("h_m_bjet1",
                                "Leading b-tagged jet mass (after all cuts);m_{j} [GeV];Events",
                                100, 0.0, 100.0);

   TH1F *h_pt_bjet2  = new TH1F("h_pt_bjet2",
                                "Subleading b-tagged jet p_{T} (after all cuts);p_{T} [GeV];Events",
                                100, 0.0, 1000.0);
   TH1F *h_eta_bjet2 = new TH1F("h_eta_bjet2",
                                "Subleading b-tagged jet #eta (after all cuts);#eta;Events",
                                60, -5.0, 5.0);
   TH1F *h_m_bjet2   = new TH1F("h_m_bjet2",
                                "Subleading b-tagged jet mass (after all cuts);m_{j} [GeV];Events",
                                100, 0.0, 100.0);

   // b-tag multiplicity after all cuts
   TH1F *h_Nbjets_after = new TH1F("h_Nbjets_after",
                                   "b-tagged jet multiplicity (after all cuts);N_{b\\text{-jets}};Events",
                                   10, 0, 10);

   // event hardness variables after all cuts
   TH1F *h_HT    = new TH1F("h_HT",
                            "H_{T} = #Sigma p_{T}^{jets} (after all cuts);H_{T} [GeV];Events",
                            100, 0.0, 2000.0);
   TH1F *h_HT_2b = new TH1F("h_HT_2b",
                            "H_{T}^{2b} = p_{T}(b_{1})+p_{T}(b_{2}) (after all cuts);H_{T}^{2b} [GeV];Events",
                            100, 0.0, 2000.0);
   TH1F *h_ST    = new TH1F("h_ST",
                            "S_{T} = H_{T} + p_{T}^{miss} (after all cuts);S_{T} [GeV];Events",
                            100, 0.0, 2500.0);

   // MET-related angles after all cuts
   TH1F *h_dphi_MET_bb = new TH1F("h_dphi_MET_bb",
                                  "|#Delta#phi(MET,bb)| (after all cuts);|#Delta#phi(MET,bb)|;Events",
                                  64, 0.0, TMath::Pi());
   TH1F *h_dphi_MET_J1 = new TH1F("h_dphi_MET_J1",
                                  "|#Delta#phi(MET,J_{1})| (after all cuts);|#Delta#phi(MET,J_{1})|;Events",
                                  64, 0.0, TMath::Pi());

   // ΔR(J1, e1) and ΔR(J1, μ1) (diagnostics)
   TH1F *h_dR_J1_e1  = new TH1F("h_dR_J1_e1",
                                "#DeltaR(J_{1}, e_{1});#DeltaR(J_{1},e_{1});Events",
                                60, 0.0, 6.0);
   TH1F *h_dR_J1_mu1 = new TH1F("h_dR_J1_mu1",
                                "#DeltaR(J_{1}, #mu_{1});#DeltaR(J_{1},#mu_{1});Events",
                                60, 0.0, 6.0);

   // ΔR(J1, e1) and ΔR(J1, μ1) AFTER cleaning (using cleaned leading jet)
   TH1F *h_dR_J1_e1_clean  = new TH1F("h_dR_J1_e1_clean",
                                   "#DeltaR(J_{1}^{clean}, e_{1});#DeltaR(J_{1}^{clean},e_{1});Events",
                                   60, 0.0, 6.0);
   TH1F *h_dR_J1_mu1_clean = new TH1F("h_dR_J1_mu1_clean",
                                   "#DeltaR(J_{1}^{clean}, #mu_{1});#DeltaR(J_{1}^{clean},#mu_{1});Events",
                                   60, 0.0, 6.0);

   const float PI      = 3.14159265;
   const float PT_CUT  = 20.0;
   const float ETA_CUT = 2.5;

   Long64_t nbytes = 0, nb = 0;

   struct PtComparator {
      const float* pt;
      PtComparator(const float* p) : pt(p) {}
      bool operator()(int a, int b) const { return pt[a] > pt[b]; }
   };

   // Cutflow counters (unweighted counts)
   Long64_t nRawEvents = nentries;
   Long64_t nPassCut1  = 0; // veto loose leptons
   Long64_t nPassCut2  = 0; // Njets >= 2
   Long64_t nPassCut3  = 0; // pT(J1) > 100 GeV
   Long64_t nPassCut4  = 0; // 2 b-tagged jets
   Long64_t nPassCut5  = 0; // MET < 140 GeV

   // Helper for b-tag WP 
   auto isBTagged = [&](int idx) -> bool {
      // Use existing branch: Jet_btagUParTAK4probbb
      return (Jet_btagUParTAK4probbb[idx] > 0.38);  // adjust threshold to your WP
   };

   // Decide once if we do GEN based on file name (starts with "GluGluH")
   std::string baseName = "";
   bool isGluGluH = false;
   
   TFile *curFile = fChain->GetCurrentFile();
   if (curFile) {
     const char *fullName = curFile->GetName();   // e.g. /path/to/GluGluH-signal.root
     baseName = gSystem->BaseName(fullName);      // e.g. GluGluH-signal.root
     if (baseName.compare(0, 7, "GluGluH") == 0) {
       isGluGluH = true;
     }
   }
   
   doGen = isGluGluH;
   
   if (doGen) {
     std::cout << "GEN analysis ENABLED (file \"" << baseName
	       << "\" starts with \"GluGluH\")" << std::endl;
   } else {
     std::cout << "GEN analysis DISABLED (file does not start with \"GluGluH\")" << std::endl;
   }
   
   // =========================
   // Event loop
   // =========================
   for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // ==========================================================
      // GEN-LEVEL analysis (only if doGen == true) – UNWEIGHTED
      // ==========================================================
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

               while (dphi_A >  PI)  dphi_A -= 2.0f*PI;
               while (dphi_A <= -PI) dphi_A += 2.0f*PI;

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

               while (dphi_b >  PI)  dphi_b -= 2.0f*PI;
               while (dphi_b <= -PI) dphi_b += 2.0f*PI;

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
         } // if valid Higgs
      } // end if(doGen)
      
      // ===========================================
      // RECO objects
      // ===========================================
      std::vector<int> ele_idx_pass;
      std::vector<int> mu_idx_pass;
      std::vector<int> jet_idx_pass;

      int nLooseMu  = 0;
      int nLooseEle = 0;
      int nLooseTau = 0; // placeholder (no tau branches yet)

      // ---- Muons ----
      for (int i = 0; i < nMuon; ++i) {
         float pt  = Muon_pt[i];
         float eta = Muon_eta[i];

         bool passLooseMuKin = (pt > 10.0 && std::fabs(eta) < 2.5);
         bool passLooseMuId  = (Muon_looseId[i] != 0);
         if (passLooseMuKin && passLooseMuId) {
            ++nLooseMu;
         }

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Muon_tightId[i] != 0);
         bool passIso = (Muon_pfRelIso04_all[i] < 0.15);

         if (!passKin) continue;
         if (!passId)  continue;
         if (!passIso) continue;

         mu_idx_pass.push_back(i);
      }

      // ---- Electrons ----
      for (int i = 0; i < nElectron; ++i) {
         float pt  = Electron_pt[i];
         float eta = Electron_eta[i];

         bool passLooseEleKin = (pt > 10.0 && std::fabs(eta) < 2.5);
         bool passLooseEleId  = (Electron_cutBased[i] >= 1); // 'veto' WP
         if (passLooseEleKin && passLooseEleId) {
            ++nLooseEle;
         }

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Electron_cutBased[i] >= 3);         // medium/tight
         bool passIso = (Electron_pfRelIso03_all[i] < 0.15);

         if (!passKin) continue;
         if (!passId)  continue;
         if (!passIso) continue;

         ele_idx_pass.push_back(i);
      }

      // ---- Taus (placeholder) ----
      // If you add tau branches, fill nLooseTau here.

      // ---- Jets ----
      for (int i = 0; i < nJet; ++i) {
         float pt  = Jet_pt[i];
         float eta = Jet_eta[i];

         bool passKin = (pt >= PT_CUT && std::fabs(eta) <= ETA_CUT);
         bool passId  = (Jet_passJetIdTightLepVeto[i] != 0);

         if (!passKin) continue;
         if (!passId)  continue;

         jet_idx_pass.push_back(i);
      }

      // ---------- BEFORE-SELECTION HISTOGRAMS (WEIGHTED) ----------
      h_nMu->Fill( (int)mu_idx_pass.size(),  wgt );
      h_nEle->Fill( (int)ele_idx_pass.size(), wgt );
      h_nJet->Fill( (int)jet_idx_pass.size(), wgt );

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

      // Muons
      if (mu_idx_pass.size() >= 1) {
         int i1 = mu_idx_pass[0];
         h_pt_mu1->Fill(Muon_pt[i1],  wgt);
         h_eta_mu1->Fill(Muon_eta[i1], wgt);
         h_phi_mu1->Fill(Muon_phi[i1], wgt);
      }
      if (mu_idx_pass.size() >= 2) {
         int i2 = mu_idx_pass[1];
         h_pt_mu2->Fill(Muon_pt[i2],  wgt);
         h_eta_mu2->Fill(Muon_eta[i2], wgt);
         h_phi_mu2->Fill(Muon_phi[i2], wgt);
      }

      // Electrons
      if (ele_idx_pass.size() >= 1) {
         int i1 = ele_idx_pass[0];
         h_pt_e1->Fill(Electron_pt[i1],  wgt);
         h_eta_e1->Fill(Electron_eta[i1], wgt);
         h_phi_e1->Fill(Electron_phi[i1], wgt);
      }
      if (ele_idx_pass.size() >= 2) {
         int i2 = ele_idx_pass[1];
         h_pt_e2->Fill(Electron_pt[i2],  wgt);
         h_eta_e2->Fill(Electron_eta[i2], wgt);
         h_phi_e2->Fill(Electron_phi[i2], wgt);
      }

      // Jets
      if (jet_idx_pass.size() >= 1) {
         int j1 = jet_idx_pass[0];
         h_pt_J1->Fill(Jet_pt[j1],  wgt);
         h_eta_J1->Fill(Jet_eta[j1], wgt);
         h_phi_J1->Fill(Jet_phi[j1], wgt);
      }
      if (jet_idx_pass.size() >= 2) {
         int j2 = jet_idx_pass[1];
         h_pt_J2->Fill(Jet_pt[j2],  wgt);
         h_eta_J2->Fill(Jet_eta[j2], wgt);
         h_phi_J2->Fill(Jet_phi[j2], wgt);
      }
      if (jet_idx_pass.size() >= 3) {
         int j3 = jet_idx_pass[2];
         h_pt_J3->Fill(Jet_pt[j3],  wgt);
         h_eta_J3->Fill(Jet_eta[j3], wgt);
         h_phi_J3->Fill(Jet_phi[j3], wgt);
      }
      if (jet_idx_pass.size() >= 4) {
         int j4 = jet_idx_pass[3];
         h_pt_J4->Fill(Jet_pt[j4],  wgt);
         h_eta_J4->Fill(Jet_eta[j4], wgt);
         h_phi_J4->Fill(Jet_phi[j4], wgt);
      }

      // MET (before selection)
      h_MET->Fill(PuppiMET_pt,  wgt);
      h_MET_phi->Fill(PuppiMET_phi, wgt);

      // -----------------------------------------------
      // ΔR(J1,e1) and ΔR(J1,μ1) BEFORE cleaning
      // -----------------------------------------------
      if (!jet_idx_pass.empty()) {
         int j1_idx_raw = jet_idx_pass[0];

         if (!ele_idx_pass.empty()) {
            int e1_idx = ele_idx_pass[0];
            float dEta = Jet_eta[j1_idx_raw] - Electron_eta[e1_idx];
            float dPhi = Jet_phi[j1_idx_raw] - Electron_phi[e1_idx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            h_dR_J1_e1->Fill(dR, wgt);
         }

         if (!mu_idx_pass.empty()) {
            int m1_idx = mu_idx_pass[0];
            float dEta = Jet_eta[j1_idx_raw] - Muon_eta[m1_idx];
            float dPhi = Jet_phi[j1_idx_raw] - Muon_phi[m1_idx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            h_dR_J1_mu1->Fill(dR, wgt);
         }
      }

      // -------------------------------------------------------------
      // jet–lepton cleaning: remove jets with ΔR<0.4 to ANY lepton
      // -------------------------------------------------------------
      std::vector<int> jet_idx_clean;
      for (size_t ijet = 0; ijet < jet_idx_pass.size(); ++ijet) {
         int jidx = jet_idx_pass[ijet];
         bool overlap = false;

         // check vs selected muons
         for (size_t imu = 0; imu < mu_idx_pass.size(); ++imu) {
            int midx = mu_idx_pass[imu];
            float dEta = Jet_eta[jidx] - Muon_eta[midx];
            float dPhi = Jet_phi[jidx] - Muon_phi[midx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            if (dR < 0.4) { overlap = true; break; }
         }
         if (overlap) continue;

         // check vs selected electrons
         for (size_t ie = 0; ie < ele_idx_pass.size(); ++ie) {
            int eidx = ele_idx_pass[ie];
            float dEta = Jet_eta[jidx] - Electron_eta[eidx];
            float dPhi = Jet_phi[jidx] - Electron_phi[eidx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            if (dR < 0.4) { overlap = true; break; }
         }
         if (overlap) continue;

         // keep jet
         jet_idx_clean.push_back(jidx);
      }

      // -----------------------------------------------
      // ΔR(J1_clean, e1) and ΔR(J1_clean, μ1) AFTER cleaning
      // -----------------------------------------------
      if (!jet_idx_clean.empty()) {
         int j1_idx_clean = jet_idx_clean[0];

         // with leading electron (if any)
         if (!ele_idx_pass.empty()) {
            int e1_idx = ele_idx_pass[0];
            float dEta = Jet_eta[j1_idx_clean] - Electron_eta[e1_idx];
            float dPhi = Jet_phi[j1_idx_clean] - Electron_phi[e1_idx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            h_dR_J1_e1_clean->Fill(dR, wgt);
         }

         // with leading muon (if any)
         if (!mu_idx_pass.empty()) {
            int m1_idx = mu_idx_pass[0];
            float dEta = Jet_eta[j1_idx_clean] - Muon_eta[m1_idx];
            float dPhi = Jet_phi[j1_idx_clean] - Muon_phi[m1_idx];
            while (dPhi >  PI)  dPhi -= 2.0f*PI;
            while (dPhi <= -PI) dPhi += 2.0f*PI;
            float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            h_dR_J1_mu1_clean->Fill(dR, wgt);
         }
      }

      // =========================
      // EVENT SELECTION (cutflow, unweighted)
      // =========================

      int nLooseLep = nLooseMu + nLooseEle; // tau placeholder ignored for now

      // Cut1: veto loose leptons
      if (nLooseLep != 0) continue;
      ++nPassCut1;

      // Cut2: at least 2 jets (after cleaning)
      if (jet_idx_clean.size() < 2) continue;
      ++nPassCut2;

      // Cut3: leading jet pT > 100 GeV
      int j1_idx = jet_idx_clean[0];
      float ptJ1 = Jet_pt[j1_idx];
      if (ptJ1 < 100.0) continue;
      ++nPassCut3;

      // Cut4: two b-tagged jets (double b-tag)
      int j2_idx = jet_idx_clean[1];

      bool j1_b = isBTagged(j1_idx);
      bool j2_b = isBTagged(j2_idx);
      
      if (!j1_b || !j2_b) continue;
      ++nPassCut4;

      // Cut5: MET < 140 GeV
      if (PuppiMET_pt > 140.0) continue;
      ++nPassCut5;

      // =========================
      // FINAL EVENT-LEVEL KINEMATICS (after all cuts)
      // =========================

      // |Δφ(J1,J2)| AFTER all cuts (weighted)
      float dphi = Jet_phi[j1_idx] - Jet_phi[j2_idx];
      while (dphi >  PI)  dphi -= 2.0f*PI;
      while (dphi <= -PI) dphi += 2.0f*PI;
      float abs_dphi = std::fabs(dphi);
      h_dphi_J1J2->Fill(abs_dphi, wgt);

      // Reconstruct 2b system (Higgs candidate) from the two b-tagged jets (weighted)
      TLorentzVector j1p4, j2p4, Hcand;
      j1p4.SetPtEtaPhiM(Jet_pt[j1_idx],
                        Jet_eta[j1_idx],
                        Jet_phi[j1_idx],
                        Jet_mass[j1_idx]);

      j2p4.SetPtEtaPhiM(Jet_pt[j2_idx],
                        Jet_eta[j2_idx],
                        Jet_phi[j2_idx],
                        Jet_mass[j2_idx]);

      Hcand = j1p4 + j2p4;

      h_m_2b->Fill(Hcand.M(),    wgt);
      h_pt_2b->Fill(Hcand.Pt(),  wgt);
      h_eta_2b->Fill(Hcand.Eta(), wgt);

      // MET AFTER ALL CUTS
      h_MET_after->Fill(PuppiMET_pt, wgt);

      // b-tag multiplicity after all cuts & b-jet ordering
      std::vector<int> bjet_indices;
      for (size_t ii = 0; ii < jet_idx_clean.size(); ++ii) {
         int jidx = jet_idx_clean[ii];
         if (isBTagged(jidx)) {
            bjet_indices.push_back(jidx);
         }
      }

      h_Nbjets_after->Fill((int)bjet_indices.size(), wgt);

      if (!bjet_indices.empty()) {
         PtComparator cmpJetAll(Jet_pt);
         std::sort(bjet_indices.begin(), bjet_indices.end(), cmpJetAll);
      }

      if (bjet_indices.size() >= 1) {
         int b1 = bjet_indices[0];
         h_pt_bjet1->Fill(Jet_pt[b1],   wgt);
         h_eta_bjet1->Fill(Jet_eta[b1], wgt);
         h_m_bjet1->Fill(Jet_mass[b1],  wgt);
      }
      if (bjet_indices.size() >= 2) {
         int b2 = bjet_indices[1];
         h_pt_bjet2->Fill(Jet_pt[b2],   wgt);
         h_eta_bjet2->Fill(Jet_eta[b2], wgt);
         h_m_bjet2->Fill(Jet_mass[b2],  wgt);
      }

      // Event hardness: HT = sum pT(jets), HT_2b = pT(b1)+pT(b2), ST = HT + MET
      double HT = 0.0;
      for (size_t ii = 0; ii < jet_idx_clean.size(); ++ii) {
         int jidx = jet_idx_clean[ii];
         HT += Jet_pt[jidx];
      }

      double HT_2b = 0.0;
      if (bjet_indices.size() >= 2) {
         int b1 = bjet_indices[0];
         int b2 = bjet_indices[1];
         HT_2b = Jet_pt[b1] + Jet_pt[b2];
      } else {
         // fallback: use the two b-tag jets that passed cut4
         HT_2b = Jet_pt[j1_idx] + Jet_pt[j2_idx];
      }

      double ST = HT + PuppiMET_pt;

      h_HT->Fill(HT, wgt);
      h_HT_2b->Fill(HT_2b, wgt);
      h_ST->Fill(ST, wgt);

      // MET-related angles
      double dphi_MET_bb = PuppiMET_phi - Hcand.Phi();
      while (dphi_MET_bb >  PI)  dphi_MET_bb -= 2.0*PI;
      while (dphi_MET_bb <= -PI) dphi_MET_bb += 2.0*PI;
      h_dphi_MET_bb->Fill(std::fabs(dphi_MET_bb), wgt);

      double dphi_MET_J1 = PuppiMET_phi - Jet_phi[j1_idx];
      while (dphi_MET_J1 >  PI)  dphi_MET_J1 -= 2.0*PI;
      while (dphi_MET_J1 <= -PI) dphi_MET_J1 += 2.0*PI;
      h_dphi_MET_J1->Fill(std::fabs(dphi_MET_J1), wgt);

   } // end event loop

   // =======================
   // Cutflow table (unweighted Nevents + WNevents)
   // =======================
   std::cout << "\nEvent flow (cutflow):\n";
   std::cout << "---------------------------------------------------------------\n";
   std::cout << std::left
             << std::setw(22) << "Cut"
             << std::setw(15) << "Nevents"
             << std::setw(15) << "ε_cut"
             << std::setw(15) << "WNevents"
             << "\n";
   std::cout << "---------------------------------------------------------------\n";

   std::cout << std::fixed;

   auto eff = [&](Long64_t n) -> double {
      return (nRawEvents > 0) ? double(n) / double(nRawEvents) : 0.0;
   };
   auto wNev = [&](Long64_t n) -> double {
      return double(n) * wgt;
   };

   auto printRow = [&](const char* label, Long64_t n, double eps) {
      std::cout << std::left << std::setw(22) << label;

      // Nevents with 1 decimal
      std::cout << std::right << std::setprecision(1)
                << std::setw(15) << static_cast<double>(n);

      // ε_cut with 3 decimals
      std::cout << std::setprecision(3)
                << std::setw(15) << eps;

      // WNevents with 1 decimal
      std::cout << std::setprecision(1)
                << std::setw(15) << wNev(n)
                << "\n";
   };

   printRow("Raw",                nRawEvents, 1.000);
   printRow("Cut1: veto loose ℓ", nPassCut1, eff(nPassCut1));
   printRow("Cut2: N_{jet}>=2",   nPassCut2, eff(nPassCut2));
   printRow("Cut3: p_{T}(J1)>100",nPassCut3, eff(nPassCut3));
   printRow("Cut4: 2 b-tag jets", nPassCut4, eff(nPassCut4));
   printRow("Cut5: MET<140 GeV",  nPassCut5, eff(nPassCut5));

   std::cout << "---------------------------------------------------------------\n";

   // =======================
   // Save histograms
   // =======================
   TFile *f = new TFile("ggHAnalysis_plots.root", "RECREATE");

   // GEN (unweighted)
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

   // RECO (before selection, weighted)
   h_nEle->Write();
   h_nMu->Write();
   h_nJet->Write();

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

   // ΔR(J1,ℓ1) histos (before & after cleaning)
   h_dR_J1_e1->Write();
   h_dR_J1_mu1->Write();
   h_dR_J1_e1_clean->Write();
   h_dR_J1_mu1_clean->Write();

   // RECO after all cuts (weighted)
   h_dphi_J1J2->Write();

   h_m_2b->Write();
   h_pt_2b->Write();
   h_eta_2b->Write();

   h_MET_after->Write();

   h_pt_bjet1->Write();
   h_eta_bjet1->Write();
   h_m_bjet1->Write();

   h_pt_bjet2->Write();
   h_eta_bjet2->Write();
   h_m_bjet2->Write();

   h_Nbjets_after->Write();

   h_HT->Write();
   h_HT_2b->Write();
   h_ST->Write();

   h_dphi_MET_bb->Write();
   h_dphi_MET_J1->Write();

   f->Close();

   std::cout << "Wrote GEN (unweighted) + RECO (weighted) histograms to ggHAnalysis_plots.root" << std::endl;
}
