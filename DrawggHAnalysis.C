

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

// Detach histogram from file so it remains valid after file close
static void Detach(TH1 *h) {
   if (h) h->SetDirectory(nullptr);
}

void DrawggHAnalysis(int mode = 0)
{
   gStyle->SetOptStat(1111);

   TFile *f = TFile::Open("ggHAnalysis_plots.root", "READ");
   if (!f || f->IsZombie()) {
      std::cout << "ERROR: cannot open ggHAnalysis_plots.root\n";
      return;
   }

   std::cout << "Opened file: " << f->GetName() << "\n";
   f->ls();

   // =======================
   // Get GEN histograms
   // =======================
   TH1F *h_pt_H   = (TH1F*)f->Get("h_pt_H");
   TH1F *h_eta_H  = (TH1F*)f->Get("h_eta_H");
   TH1F *h_phi_H  = (TH1F*)f->Get("h_phi_H");
   TH1F *h_m_H    = (TH1F*)f->Get("h_m_H");

   TH1F *h_pt_A1  = (TH1F*)f->Get("h_pt_A1");
   TH1F *h_eta_A1 = (TH1F*)f->Get("h_eta_A1");
   TH1F *h_phi_A1 = (TH1F*)f->Get("h_phi_A1");
   TH1F *h_m_A1   = (TH1F*)f->Get("h_m_A1");

   TH1F *h_pt_A2  = (TH1F*)f->Get("h_pt_A2");
   TH1F *h_eta_A2 = (TH1F*)f->Get("h_eta_A2");
   TH1F *h_phi_A2 = (TH1F*)f->Get("h_phi_A2");
   TH1F *h_m_A2   = (TH1F*)f->Get("h_m_A2");

   TH1F *h_dR_AA    = (TH1F*)f->Get("h_dR_AA");
   TH1F *h_dR_bb_A1 = (TH1F*)f->Get("h_dR_bb_A1");
   TH1F *h_dR_bb_A2 = (TH1F*)f->Get("h_dR_bb_A2");

   TH1F *h_pt_b1  = (TH1F*)f->Get("h_pt_b1");
   TH1F *h_eta_b1 = (TH1F*)f->Get("h_eta_b1");
   TH1F *h_phi_b1 = (TH1F*)f->Get("h_phi_b1");
   TH1F *h_m_b1   = (TH1F*)f->Get("h_m_b1");

   TH1F *h_pt_b2  = (TH1F*)f->Get("h_pt_b2");
   TH1F *h_eta_b2 = (TH1F*)f->Get("h_eta_b2");
   TH1F *h_phi_b2 = (TH1F*)f->Get("h_phi_b2");
   TH1F *h_m_b2   = (TH1F*)f->Get("h_m_b2");

   TH1F *h_pt_b3  = (TH1F*)f->Get("h_pt_b3");
   TH1F *h_eta_b3 = (TH1F*)f->Get("h_eta_b3");
   TH1F *h_phi_b3 = (TH1F*)f->Get("h_phi_b3");
   TH1F *h_m_b3   = (TH1F*)f->Get("h_m_b3");

   TH1F *h_pt_b4  = (TH1F*)f->Get("h_pt_b4");
   TH1F *h_eta_b4 = (TH1F*)f->Get("h_eta_b4");
   TH1F *h_phi_b4 = (TH1F*)f->Get("h_phi_b4");
   TH1F *h_m_b4   = (TH1F*)f->Get("h_m_b4");

   // Detach GEN
   Detach(h_pt_H); Detach(h_eta_H); Detach(h_phi_H); Detach(h_m_H);
   Detach(h_pt_A1); Detach(h_eta_A1); Detach(h_phi_A1); Detach(h_m_A1);
   Detach(h_pt_A2); Detach(h_eta_A2); Detach(h_phi_A2); Detach(h_m_A2);
   Detach(h_dR_AA); Detach(h_dR_bb_A1); Detach(h_dR_bb_A2);
   Detach(h_pt_b1); Detach(h_eta_b1); Detach(h_phi_b1); Detach(h_m_b1);
   Detach(h_pt_b2); Detach(h_eta_b2); Detach(h_phi_b2); Detach(h_m_b2);
   Detach(h_pt_b3); Detach(h_eta_b3); Detach(h_phi_b3); Detach(h_m_b3);
   Detach(h_pt_b4); Detach(h_eta_b4); Detach(h_phi_b4); Detach(h_m_b4);

   // =======================
   // Get RECO-before histograms
   // IMPORTANT: these are TH1I in the file, so retrieve as TH1*
   // =======================
   TH1  *h_nEle = (TH1*)f->Get("h_nEle");
   TH1  *h_nMu  = (TH1*)f->Get("h_nMu");
   TH1  *h_nJet = (TH1*)f->Get("h_nJet");
   TH1  *h_nJet_clean = (TH1*)f->Get("h_nJet_clean");

   TH1F *h_pt_e1  = (TH1F*)f->Get("h_pt_e1");
   TH1F *h_eta_e1 = (TH1F*)f->Get("h_eta_e1");
   TH1F *h_phi_e1 = (TH1F*)f->Get("h_phi_e1");

   TH1F *h_pt_e2  = (TH1F*)f->Get("h_pt_e2");
   TH1F *h_eta_e2 = (TH1F*)f->Get("h_eta_e2");
   TH1F *h_phi_e2 = (TH1F*)f->Get("h_phi_e2");

   TH1F *h_pt_mu1  = (TH1F*)f->Get("h_pt_mu1");
   TH1F *h_eta_mu1 = (TH1F*)f->Get("h_eta_mu1");
   TH1F *h_phi_mu1 = (TH1F*)f->Get("h_phi_mu1");

   TH1F *h_pt_mu2  = (TH1F*)f->Get("h_pt_mu2");
   TH1F *h_eta_mu2 = (TH1F*)f->Get("h_eta_mu2");
   TH1F *h_phi_mu2 = (TH1F*)f->Get("h_phi_mu2");

   TH1F *h_pt_J1  = (TH1F*)f->Get("h_pt_J1");
   TH1F *h_eta_J1 = (TH1F*)f->Get("h_eta_J1");
   TH1F *h_phi_J1 = (TH1F*)f->Get("h_phi_J1");

   TH1F *h_pt_J2  = (TH1F*)f->Get("h_pt_J2");
   TH1F *h_eta_J2 = (TH1F*)f->Get("h_eta_J2");
   TH1F *h_phi_J2 = (TH1F*)f->Get("h_phi_J2");

   TH1F *h_pt_J3  = (TH1F*)f->Get("h_pt_J3");
   TH1F *h_eta_J3 = (TH1F*)f->Get("h_eta_J3");
   TH1F *h_phi_J3 = (TH1F*)f->Get("h_phi_J3");

   TH1F *h_pt_J4  = (TH1F*)f->Get("h_pt_J4");
   TH1F *h_eta_J4 = (TH1F*)f->Get("h_eta_J4");
   TH1F *h_phi_J4 = (TH1F*)f->Get("h_phi_J4");

   TH1F *h_MET     = (TH1F*)f->Get("h_MET");
   TH1F *h_MET_phi = (TH1F*)f->Get("h_MET_phi");

   TH1F *h_dR_J1_e1        = (TH1F*)f->Get("h_dR_J1_e1");
   TH1F *h_dR_J1_mu1       = (TH1F*)f->Get("h_dR_J1_mu1");
   TH1F *h_dR_J1_e1_clean  = (TH1F*)f->Get("h_dR_J1_e1_clean");
   TH1F *h_dR_J1_mu1_clean = (TH1F*)f->Get("h_dR_J1_mu1_clean");

   // Detach RECO-before
   Detach(h_nEle); Detach(h_nMu); Detach(h_nJet); Detach(h_nJet_clean);
   Detach(h_pt_e1); Detach(h_eta_e1); Detach(h_phi_e1);
   Detach(h_pt_e2); Detach(h_eta_e2); Detach(h_phi_e2);
   Detach(h_pt_mu1); Detach(h_eta_mu1); Detach(h_phi_mu1);
   Detach(h_pt_mu2); Detach(h_eta_mu2); Detach(h_phi_mu2);
   Detach(h_pt_J1); Detach(h_eta_J1); Detach(h_phi_J1);
   Detach(h_pt_J2); Detach(h_eta_J2); Detach(h_phi_J2);
   Detach(h_pt_J3); Detach(h_eta_J3); Detach(h_phi_J3);
   Detach(h_pt_J4); Detach(h_eta_J4); Detach(h_phi_J4);
   Detach(h_MET); Detach(h_MET_phi);
   Detach(h_dR_J1_e1); Detach(h_dR_J1_mu1);
   Detach(h_dR_J1_e1_clean); Detach(h_dR_J1_mu1_clean);

   // =======================
   // Get RECO-after histograms
   // =======================
   TH1F *h_dphi_J1J2 = (TH1F*)f->Get("h_dphi_J1J2");
   TH1F *h_dR_J1J2   = (TH1F*)f->Get("h_dR_J1J2");
   TH1F *h_abs_dm_J1J2 = (TH1F*)f->Get("h_abs_dm_J1J2"); // NEW

   TH1F *h_min_dphi_MET_Jet = (TH1F*)f->Get("h_min_dphi_MET_Jet");

   TH1F *h_m_2b   = (TH1F*)f->Get("h_m_2b");
   TH1F *h_pt_2b  = (TH1F*)f->Get("h_pt_2b");
   TH1F *h_eta_2b = (TH1F*)f->Get("h_eta_2b");

   TH1F *h_MET_after     = (TH1F*)f->Get("h_MET_after");
   TH1F *h_MET_phi_after = (TH1F*)f->Get("h_MET_phi_after");

   TH1F *h_pt_bjet1  = (TH1F*)f->Get("h_pt_bjet1");
   TH1F *h_eta_bjet1 = (TH1F*)f->Get("h_eta_bjet1");
   TH1F *h_m_bjet1   = (TH1F*)f->Get("h_m_bjet1");

   TH1F *h_pt_bjet2  = (TH1F*)f->Get("h_pt_bjet2");
   TH1F *h_eta_bjet2 = (TH1F*)f->Get("h_eta_bjet2");
   TH1F *h_m_bjet2   = (TH1F*)f->Get("h_m_bjet2");

   TH1F *h_Nbjets_after = (TH1F*)f->Get("h_Nbjets_after");
   TH1F *h_HT           = (TH1F*)f->Get("h_HT");

   TH1F *h_dphi_MET_bb = (TH1F*)f->Get("h_dphi_MET_bb");
   TH1F *h_dphi_MET_J1 = (TH1F*)f->Get("h_dphi_MET_J1");

   // Detach RECO-after
   Detach(h_dphi_J1J2); Detach(h_dR_J1J2); Detach(h_abs_dm_J1J2);
   Detach(h_min_dphi_MET_Jet);
   Detach(h_m_2b); Detach(h_pt_2b); Detach(h_eta_2b);
   Detach(h_MET_after); Detach(h_MET_phi_after);
   Detach(h_pt_bjet1); Detach(h_eta_bjet1); Detach(h_m_bjet1);
   Detach(h_pt_bjet2); Detach(h_eta_bjet2); Detach(h_m_bjet2);
   Detach(h_Nbjets_after); Detach(h_HT);
   Detach(h_dphi_MET_bb); Detach(h_dphi_MET_J1);

   // Close file
   f->Close();
   delete f;

   // =======================
   // MODE 0: GEN
   // =======================
   if (mode == 0) {
      TCanvas *c_H = new TCanvas("c_H", "Higgs (GEN)", 900, 800);
      c_H->Divide(2,2);
      c_H->cd(1); if (h_pt_H)  h_pt_H->Draw("hist");
      c_H->cd(2); if (h_eta_H) h_eta_H->Draw("hist");
      c_H->cd(3); if (h_phi_H) h_phi_H->Draw("hist");
      c_H->cd(4); if (h_m_H)   h_m_H->Draw("hist");

      TCanvas *c_A = new TCanvas("c_A", "A_{1} and A_{2} (GEN)", 1200, 900);
      c_A->Divide(3,3);
      c_A->cd(1); if (h_pt_A1)  h_pt_A1->Draw("hist");
      c_A->cd(2); if (h_eta_A1) h_eta_A1->Draw("hist");
      c_A->cd(3); if (h_phi_A1) h_phi_A1->Draw("hist");
      c_A->cd(4); if (h_pt_A2)  h_pt_A2->Draw("hist");
      c_A->cd(5); if (h_eta_A2) h_eta_A2->Draw("hist");
      c_A->cd(6); if (h_phi_A2) h_phi_A2->Draw("hist");
      c_A->cd(7); if (h_m_A1)   h_m_A1->Draw("hist");
      c_A->cd(8); if (h_m_A2)   h_m_A2->Draw("hist");
      c_A->cd(9); if (h_dR_AA)  h_dR_AA->Draw("hist");

      TCanvas *c_dR_bb = new TCanvas("c_dR_bb", "#DeltaR(b,b) from A_{1} and A_{2}", 900, 600);
      c_dR_bb->Divide(2,1);
      c_dR_bb->cd(1); if (h_dR_bb_A1) h_dR_bb_A1->Draw("hist");
      c_dR_bb->cd(2); if (h_dR_bb_A2) h_dR_bb_A2->Draw("hist");

      TCanvas *c_4b = new TCanvas("c_4b", "b_{1..4} (GEN)", 1600, 1200);
      c_4b->Divide(4,4);
      c_4b->cd(1);  if (h_pt_b1)  h_pt_b1->Draw("hist");
      c_4b->cd(2);  if (h_eta_b1) h_eta_b1->Draw("hist");
      c_4b->cd(3);  if (h_phi_b1) h_phi_b1->Draw("hist");
      c_4b->cd(4);  if (h_m_b1)   h_m_b1->Draw("hist");
      c_4b->cd(5);  if (h_pt_b2)  h_pt_b2->Draw("hist");
      c_4b->cd(6);  if (h_eta_b2) h_eta_b2->Draw("hist");
      c_4b->cd(7);  if (h_phi_b2) h_phi_b2->Draw("hist");
      c_4b->cd(8);  if (h_m_b2)   h_m_b2->Draw("hist");
      c_4b->cd(9);  if (h_pt_b3)  h_pt_b3->Draw("hist");
      c_4b->cd(10); if (h_eta_b3) h_eta_b3->Draw("hist");
      c_4b->cd(11); if (h_phi_b3) h_phi_b3->Draw("hist");
      c_4b->cd(12); if (h_m_b3)   h_m_b3->Draw("hist");
      c_4b->cd(13); if (h_pt_b4)  h_pt_b4->Draw("hist");
      c_4b->cd(14); if (h_eta_b4) h_eta_b4->Draw("hist");
      c_4b->cd(15); if (h_phi_b4) h_phi_b4->Draw("hist");
      c_4b->cd(16); if (h_m_b4)   h_m_b4->Draw("hist");

      return;
   }

   // =======================
   // MODE 1: RECO BEFORE selection  (FIXED: full canvases, no placeholder return)
   // =======================
   if (mode == 1) {
      TCanvas *c_mult_before = new TCanvas("c_mult_before",
         "Multiplicity: muons, electrons, jets (before sel)", 1200, 400);
      c_mult_before->Divide(4,1);
      c_mult_before->cd(1); if (h_nMu)  { h_nMu->SetLineWidth(2);  h_nMu->Draw("hist"); }
      c_mult_before->cd(2); if (h_nEle) { h_nEle->SetLineWidth(2); h_nEle->Draw("hist"); }
      c_mult_before->cd(3); if (h_nJet) { h_nJet->SetLineWidth(2); h_nJet->Draw("hist"); }
      c_mult_before->cd(4); if (h_nJet_clean) { h_nJet_clean->SetLineWidth(2); h_nJet_clean->Draw("hist"); }

      TCanvas *c_ele_before = new TCanvas("c_ele_before", "Electrons (before sel)", 1200, 800);
      c_ele_before->Divide(3,2);
      c_ele_before->cd(1); if (h_pt_e1)  { h_pt_e1->SetLineWidth(2);  h_pt_e1->Draw("hist"); }
      c_ele_before->cd(2); if (h_eta_e1) { h_eta_e1->SetLineWidth(2); h_eta_e1->Draw("hist"); }
      c_ele_before->cd(3); if (h_phi_e1) { h_phi_e1->SetLineWidth(2); h_phi_e1->Draw("hist"); }
      c_ele_before->cd(4); if (h_pt_e2)  { h_pt_e2->SetLineWidth(2);  h_pt_e2->Draw("hist"); }
      c_ele_before->cd(5); if (h_eta_e2) { h_eta_e2->SetLineWidth(2); h_eta_e2->Draw("hist"); }
      c_ele_before->cd(6); if (h_phi_e2) { h_phi_e2->SetLineWidth(2); h_phi_e2->Draw("hist"); }

      TCanvas *c_mu_before = new TCanvas("c_mu_before", "Muons (before sel)", 1200, 800);
      c_mu_before->Divide(3,2);
      c_mu_before->cd(1); if (h_pt_mu1)  { h_pt_mu1->SetLineWidth(2);  h_pt_mu1->Draw("hist"); }
      c_mu_before->cd(2); if (h_eta_mu1) { h_eta_mu1->SetLineWidth(2); h_eta_mu1->Draw("hist"); }
      c_mu_before->cd(3); if (h_phi_mu1) { h_phi_mu1->SetLineWidth(2); h_phi_mu1->Draw("hist"); }
      c_mu_before->cd(4); if (h_pt_mu2)  { h_pt_mu2->SetLineWidth(2);  h_pt_mu2->Draw("hist"); }
      c_mu_before->cd(5); if (h_eta_mu2) { h_eta_mu2->SetLineWidth(2); h_eta_mu2->Draw("hist"); }
      c_mu_before->cd(6); if (h_phi_mu2) { h_phi_mu2->SetLineWidth(2); h_phi_mu2->Draw("hist"); }

      TCanvas *c_jets_before = new TCanvas("c_jets_before", "Jets J1-J4 (before sel)", 1400, 900);
      c_jets_before->Divide(4,3);
      c_jets_before->cd(1);  if (h_pt_J1)  { h_pt_J1->SetLineWidth(2);  h_pt_J1->Draw("hist"); }
      c_jets_before->cd(2);  if (h_eta_J1) { h_eta_J1->SetLineWidth(2); h_eta_J1->Draw("hist"); }
      c_jets_before->cd(3);  if (h_phi_J1) { h_phi_J1->SetLineWidth(2); h_phi_J1->Draw("hist"); }
      c_jets_before->cd(4);  if (h_pt_J2)  { h_pt_J2->SetLineWidth(2);  h_pt_J2->Draw("hist"); }
      c_jets_before->cd(5);  if (h_eta_J2) { h_eta_J2->SetLineWidth(2); h_eta_J2->Draw("hist"); }
      c_jets_before->cd(6);  if (h_phi_J2) { h_phi_J2->SetLineWidth(2); h_phi_J2->Draw("hist"); }
      c_jets_before->cd(7);  if (h_pt_J3)  { h_pt_J3->SetLineWidth(2);  h_pt_J3->Draw("hist"); }
      c_jets_before->cd(8);  if (h_eta_J3) { h_eta_J3->SetLineWidth(2); h_eta_J3->Draw("hist"); }
      c_jets_before->cd(9);  if (h_phi_J3) { h_phi_J3->SetLineWidth(2); h_phi_J3->Draw("hist"); }
      c_jets_before->cd(10); if (h_pt_J4)  { h_pt_J4->SetLineWidth(2);  h_pt_J4->Draw("hist"); }
      c_jets_before->cd(11); if (h_eta_J4) { h_eta_J4->SetLineWidth(2); h_eta_J4->Draw("hist"); }
      c_jets_before->cd(12); if (h_phi_J4) { h_phi_J4->SetLineWidth(2); h_phi_J4->Draw("hist"); }

      TCanvas *c_met_before = new TCanvas("c_met_before", "MET (before sel)", 900, 450);
      c_met_before->Divide(2,1);
      c_met_before->cd(1); if (h_MET)     { h_MET->SetLineWidth(2);     h_MET->Draw("hist"); }
      c_met_before->cd(2); if (h_MET_phi) { h_MET_phi->SetLineWidth(2); h_MET_phi->Draw("hist"); }

      TCanvas *c_dR_clean = new TCanvas("c_dR_clean",
         "#DeltaR(J1, lepton) before/after jet-lepton cleaning", 1200, 800);
      c_dR_clean->Divide(2,2);
      c_dR_clean->cd(1); if (h_dR_J1_e1)        { h_dR_J1_e1->SetLineWidth(2);        h_dR_J1_e1->Draw("hist"); }
      c_dR_clean->cd(2); if (h_dR_J1_mu1)       { h_dR_J1_mu1->SetLineWidth(2);       h_dR_J1_mu1->Draw("hist"); }
      c_dR_clean->cd(3); if (h_dR_J1_e1_clean)  { h_dR_J1_e1_clean->SetLineWidth(2);  h_dR_J1_e1_clean->Draw("hist"); }
      c_dR_clean->cd(4); if (h_dR_J1_mu1_clean) { h_dR_J1_mu1_clean->SetLineWidth(2); h_dR_J1_mu1_clean->Draw("hist"); }

      return;
   }

   // =======================
   // MODE 2: RECO AFTER all cuts
   // =======================
   if (mode == 2) {
      TCanvas *c_dphi = new TCanvas("c_dphi_J1J2", "|#Delta#phi(J1,J2)| after all cuts", 600, 500);
      if (h_dphi_J1J2) { h_dphi_J1J2->SetLineWidth(2); h_dphi_J1J2->Draw("hist"); }

      TCanvas *c_dR = new TCanvas("c_dR_J1J2", "#DeltaR(J1,J2) after all cuts", 600, 500);
      if (h_dR_J1J2) { h_dR_J1J2->SetLineWidth(2); h_dR_J1J2->Draw("hist"); }

      // NEW: |m(J1)-m(J2)|
      TCanvas *c_absdm = new TCanvas("c_abs_dm_J1J2", "|m(J1)-m(J2)| after all cuts", 600, 500);
      if (h_abs_dm_J1J2) { h_abs_dm_J1J2->SetLineWidth(2); h_abs_dm_J1J2->Draw("hist"); }
      else std::cout << "WARNING: h_abs_dm_J1J2 not found in file.\n";

      TCanvas *c_minDphi = new TCanvas("c_min_dphi_MET_Jet", "min |#Delta#phi(MET,Jet)| after all cuts", 600, 500);
      if (h_min_dphi_MET_Jet) { h_min_dphi_MET_Jet->SetLineWidth(2); h_min_dphi_MET_Jet->Draw("hist"); }

      TCanvas *c_Hreco = new TCanvas("c_Hreco", "Higgs (bb) after all cuts", 1200, 400);
      c_Hreco->Divide(3,1);
      c_Hreco->cd(1); if (h_m_2b)  { h_m_2b->SetLineWidth(2);  h_m_2b->Draw("hist"); }
      c_Hreco->cd(2); if (h_pt_2b) { h_pt_2b->SetLineWidth(2); h_pt_2b->Draw("hist"); }
      c_Hreco->cd(3); if (h_eta_2b){ h_eta_2b->SetLineWidth(2);h_eta_2b->Draw("hist"); }

      TCanvas *c_MET_after = new TCanvas("c_MET_after", "MET after all cuts", 900, 450);
      c_MET_after->Divide(2,1);
      c_MET_after->cd(1); if (h_MET_after)     { h_MET_after->SetLineWidth(2);     h_MET_after->Draw("hist"); }
      c_MET_after->cd(2); if (h_MET_phi_after) { h_MET_phi_after->SetLineWidth(2); h_MET_phi_after->Draw("hist"); }

      TCanvas *c_bjets = new TCanvas("c_bjets", "b-jets after all cuts", 1200, 800);
      c_bjets->Divide(3,2);
      c_bjets->cd(1); if (h_pt_bjet1)  { h_pt_bjet1->SetLineWidth(2);  h_pt_bjet1->Draw("hist"); }
      c_bjets->cd(2); if (h_eta_bjet1) { h_eta_bjet1->SetLineWidth(2); h_eta_bjet1->Draw("hist"); }
      c_bjets->cd(3); if (h_m_bjet1)   { h_m_bjet1->SetLineWidth(2);   h_m_bjet1->Draw("hist"); }
      c_bjets->cd(4); if (h_pt_bjet2)  { h_pt_bjet2->SetLineWidth(2);  h_pt_bjet2->Draw("hist"); }
      c_bjets->cd(5); if (h_eta_bjet2) { h_eta_bjet2->SetLineWidth(2); h_eta_bjet2->Draw("hist"); }
      c_bjets->cd(6); if (h_m_bjet2)   { h_m_bjet2->SetLineWidth(2);   h_m_bjet2->Draw("hist"); }

      TCanvas *c_HT = new TCanvas("c_HT", "HT after all cuts", 600, 500);
      if (h_HT) { h_HT->SetLineWidth(2); h_HT->Draw("hist"); }

      TCanvas *c_Nb = new TCanvas("c_Nbjets_after", "Nbjets after all cuts", 600, 500);
      if (h_Nbjets_after) { h_Nbjets_after->SetLineWidth(2); h_Nbjets_after->Draw("hist"); }

      TCanvas *c_METangles = new TCanvas("c_MET_angles", "MET angles after all cuts", 800, 400);
      c_METangles->Divide(2,1);
      c_METangles->cd(1); if (h_dphi_MET_bb) { h_dphi_MET_bb->SetLineWidth(2); h_dphi_MET_bb->Draw("hist"); }
      c_METangles->cd(2); if (h_dphi_MET_J1) { h_dphi_MET_J1->SetLineWidth(2); h_dphi_MET_J1->Draw("hist"); }

      return;
   }

   std::cout << "Unknown mode " << mode << " (use 0 GEN, 1 RECO-before, 2 RECO-after)\n";
}
