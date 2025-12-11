#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

// Small helper: detach histogram from file so it's not deleted when file closes
void Detach(TH1 *h) {
   if (h) h->SetDirectory(0);
}

// mode = 0 → draw only GEN canvases
// mode = 1 → draw only RECO BEFORE selection canvases
// mode = 2 → draw only RECO AFTER ALL CUTS canvases
void DrawggHAnalysis(int mode = 0)
{
   gStyle->SetOptStat(1111); // keep stats box

   // Open the file produced by ggHAnalysis::Loop()
   TFile *f = TFile::Open("ggHAnalysis_plots.root");
   if (!f || f->IsZombie()) {
      std::cout << "DrawggHAnalysis: could not open ggHAnalysis_plots.root" << std::endl;
      return;
   }

   std::cout << "Opened file: " << f->GetName() << std::endl;
   f->ls();  // debug: list keys

   // =======================
   // Get GEN histograms
   // =======================
   TH1F *h_pt_H   = (TH1F*)f->Get("h_pt_H");
   TH1F *h_eta_H  = (TH1F*)f->Get("h_eta_H");
   TH1F *h_phi_H  = (TH1F*)f->Get("h_phi_H");
   TH1F *h_m_H    = (TH1F*)f->Get("h_m_H");

   TH1F *h_pt_A1   = (TH1F*)f->Get("h_pt_A1");
   TH1F *h_eta_A1  = (TH1F*)f->Get("h_eta_A1");
   TH1F *h_phi_A1  = (TH1F*)f->Get("h_phi_A1");
   TH1F *h_m_A1    = (TH1F*)f->Get("h_m_A1");

   TH1F *h_pt_A2   = (TH1F*)f->Get("h_pt_A2");
   TH1F *h_eta_A2  = (TH1F*)f->Get("h_eta_A2");
   TH1F *h_phi_A2  = (TH1F*)f->Get("h_phi_A2");
   TH1F *h_m_A2    = (TH1F*)f->Get("h_m_A2");

   TH1F *h_dR_AA   = (TH1F*)f->Get("h_dR_AA");

   TH1F *h_dR_bb_A1 = (TH1F*)f->Get("h_dR_bb_A1");
   TH1F *h_dR_bb_A2 = (TH1F*)f->Get("h_dR_bb_A2");

   TH1F *h_pt_b1   = (TH1F*)f->Get("h_pt_b1");
   TH1F *h_eta_b1  = (TH1F*)f->Get("h_eta_b1");
   TH1F *h_phi_b1  = (TH1F*)f->Get("h_phi_b1");
   TH1F *h_m_b1    = (TH1F*)f->Get("h_m_b1");

   TH1F *h_pt_b2   = (TH1F*)f->Get("h_pt_b2");
   TH1F *h_eta_b2  = (TH1F*)f->Get("h_eta_b2");
   TH1F *h_phi_b2  = (TH1F*)f->Get("h_phi_b2");
   TH1F *h_m_b2    = (TH1F*)f->Get("h_m_b2");

   TH1F *h_pt_b3   = (TH1F*)f->Get("h_pt_b3");
   TH1F *h_eta_b3  = (TH1F*)f->Get("h_eta_b3");
   TH1F *h_phi_b3  = (TH1F*)f->Get("h_phi_b3");
   TH1F *h_m_b3    = (TH1F*)f->Get("h_m_b3");

   TH1F *h_pt_b4   = (TH1F*)f->Get("h_pt_b4");
   TH1F *h_eta_b4  = (TH1F*)f->Get("h_eta_b4");
   TH1F *h_phi_b4  = (TH1F*)f->Get("h_phi_b4");
   TH1F *h_m_b4    = (TH1F*)f->Get("h_m_b4");

   // Detach GEN histos
   Detach(h_pt_H);   Detach(h_eta_H);   Detach(h_phi_H);   Detach(h_m_H);
   Detach(h_pt_A1);  Detach(h_eta_A1);  Detach(h_phi_A1);  Detach(h_m_A1);
   Detach(h_pt_A2);  Detach(h_eta_A2);  Detach(h_phi_A2);  Detach(h_m_A2);
   Detach(h_dR_AA);
   Detach(h_dR_bb_A1); Detach(h_dR_bb_A2);
   Detach(h_pt_b1);  Detach(h_eta_b1);  Detach(h_phi_b1);  Detach(h_m_b1);
   Detach(h_pt_b2);  Detach(h_eta_b2);  Detach(h_phi_b2);  Detach(h_m_b2);
   Detach(h_pt_b3);  Detach(h_eta_b3);  Detach(h_phi_b3);  Detach(h_m_b3);
   Detach(h_pt_b4);  Detach(h_eta_b4);  Detach(h_phi_b4);  Detach(h_m_b4);

   // =======================
   // Get RECO histograms (before selection)
   // =======================
   TH1 *h_nEle  = (TH1*)f->Get("h_nEle");
   TH1 *h_nMu   = (TH1*)f->Get("h_nMu");
   TH1 *h_nJet  = (TH1*)f->Get("h_nJet");

   TH1F *h_pt_e1   = (TH1F*)f->Get("h_pt_e1");
   TH1F *h_eta_e1  = (TH1F*)f->Get("h_eta_e1");
   TH1F *h_phi_e1  = (TH1F*)f->Get("h_phi_e1");

   TH1F *h_pt_e2   = (TH1F*)f->Get("h_pt_e2");
   TH1F *h_eta_e2  = (TH1F*)f->Get("h_eta_e2");
   TH1F *h_phi_e2  = (TH1F*)f->Get("h_phi_e2");

   TH1F *h_pt_mu1  = (TH1F*)f->Get("h_pt_mu1");
   TH1F *h_eta_mu1 = (TH1F*)f->Get("h_eta_mu1");
   TH1F *h_phi_mu1 = (TH1F*)f->Get("h_phi_mu1");

   TH1F *h_pt_mu2  = (TH1F*)f->Get("h_pt_mu2");
   TH1F *h_eta_mu2 = (TH1F*)f->Get("h_eta_mu2");
   TH1F *h_phi_mu2 = (TH1F*)f->Get("h_phi_mu2");

   TH1F *h_pt_J1   = (TH1F*)f->Get("h_pt_J1");
   TH1F *h_eta_J1  = (TH1F*)f->Get("h_eta_J1");
   TH1F *h_phi_J1  = (TH1F*)f->Get("h_phi_J1");

   TH1F *h_pt_J2   = (TH1F*)f->Get("h_pt_J2");
   TH1F *h_eta_J2  = (TH1F*)f->Get("h_eta_J2");
   TH1F *h_phi_J2  = (TH1F*)f->Get("h_phi_J2");

   TH1F *h_pt_J3   = (TH1F*)f->Get("h_pt_J3");
   TH1F *h_eta_J3  = (TH1F*)f->Get("h_eta_J3");
   TH1F *h_phi_J3  = (TH1F*)f->Get("h_phi_J3");

   TH1F *h_pt_J4   = (TH1F*)f->Get("h_pt_J4");
   TH1F *h_eta_J4  = (TH1F*)f->Get("h_eta_J4");
   TH1F *h_phi_J4  = (TH1F*)f->Get("h_phi_J4");

   TH1F *h_MET     = (TH1F*)f->Get("h_MET");
   TH1F *h_MET_phi = (TH1F*)f->Get("h_MET_phi");

   // ΔR(J1,ℓ1) histos (before and after jet–lepton cleaning)
   TH1F *h_dR_J1_e1        = (TH1F*)f->Get("h_dR_J1_e1");
   TH1F *h_dR_J1_mu1       = (TH1F*)f->Get("h_dR_J1_mu1");
   TH1F *h_dR_J1_e1_clean  = (TH1F*)f->Get("h_dR_J1_e1_clean");
   TH1F *h_dR_J1_mu1_clean = (TH1F*)f->Get("h_dR_J1_mu1_clean");

   // Detach RECO-before histos
   Detach(h_nEle); Detach(h_nMu); Detach(h_nJet);

   Detach(h_pt_e1);  Detach(h_eta_e1);  Detach(h_phi_e1);
   Detach(h_pt_e2);  Detach(h_eta_e2);  Detach(h_phi_e2);

   Detach(h_pt_mu1); Detach(h_eta_mu1); Detach(h_phi_mu1);
   Detach(h_pt_mu2); Detach(h_eta_mu2); Detach(h_phi_mu2);

   Detach(h_pt_J1);  Detach(h_eta_J1);  Detach(h_phi_J1);
   Detach(h_pt_J2);  Detach(h_eta_J2);  Detach(h_phi_J2);
   Detach(h_pt_J3);  Detach(h_eta_J3);  Detach(h_phi_J3);
   Detach(h_pt_J4);  Detach(h_eta_J4);  Detach(h_phi_J4);

   Detach(h_MET);    Detach(h_MET_phi);

   Detach(h_dR_J1_e1);
   Detach(h_dR_J1_mu1);
   Detach(h_dR_J1_e1_clean);
   Detach(h_dR_J1_mu1_clean);

   // =======================
   // Get RECO histograms (after all cuts)
   // =======================
   TH1F *h_dphi_J1J2 = (TH1F*)f->Get("h_dphi_J1J2");

   TH1F *h_m_2b   = (TH1F*)f->Get("h_m_2b");
   TH1F *h_pt_2b  = (TH1F*)f->Get("h_pt_2b");
   TH1F *h_eta_2b = (TH1F*)f->Get("h_eta_2b");

   TH1F *h_MET_after = (TH1F*)f->Get("h_MET_after");

   TH1F *h_pt_bjet1  = (TH1F*)f->Get("h_pt_bjet1");
   TH1F *h_eta_bjet1 = (TH1F*)f->Get("h_eta_bjet1");
   TH1F *h_m_bjet1   = (TH1F*)f->Get("h_m_bjet1");

   TH1F *h_pt_bjet2  = (TH1F*)f->Get("h_pt_bjet2");
   TH1F *h_eta_bjet2 = (TH1F*)f->Get("h_eta_bjet2");
   TH1F *h_m_bjet2   = (TH1F*)f->Get("h_m_bjet2");

   TH1F *h_Nbjets_after = (TH1F*)f->Get("h_Nbjets_after");

   TH1F *h_HT    = (TH1F*)f->Get("h_HT");
   TH1F *h_HT_2b = (TH1F*)f->Get("h_HT_2b");
   TH1F *h_ST    = (TH1F*)f->Get("h_ST");

   TH1F *h_dphi_MET_bb = (TH1F*)f->Get("h_dphi_MET_bb");
   TH1F *h_dphi_MET_J1 = (TH1F*)f->Get("h_dphi_MET_J1");

   // Detach RECO-after histos
   Detach(h_dphi_J1J2);
   Detach(h_m_2b);   Detach(h_pt_2b);   Detach(h_eta_2b);
   Detach(h_MET_after);
   Detach(h_pt_bjet1); Detach(h_eta_bjet1); Detach(h_m_bjet1);
   Detach(h_pt_bjet2); Detach(h_eta_bjet2); Detach(h_m_bjet2);
   Detach(h_Nbjets_after);
   Detach(h_HT); Detach(h_HT_2b); Detach(h_ST);
   Detach(h_dphi_MET_bb); Detach(h_dphi_MET_J1);

   // Quick debug: check one histogram
   if (h_pt_H) {
      std::cout << "h_pt_H entries = " << h_pt_H->GetEntries() << std::endl;
   } else {
      std::cout << "WARNING: h_pt_H not found in file!" << std::endl;
   }

   // =======================
   // MODE 0: GEN CANVASES
   // =======================
   if (mode == 0) {

      // 1) Higgs canvas (GEN)
      TCanvas *c_H = new TCanvas("c_H", "Higgs (GEN)", 900, 800);
      c_H->Divide(2,2);
      c_H->cd(1); if (h_pt_H)  h_pt_H->Draw("hist");
      c_H->cd(2); if (h_eta_H) h_eta_H->Draw("hist");
      c_H->cd(3); if (h_phi_H) h_phi_H->Draw("hist");
      c_H->cd(4); if (h_m_H)   h_m_H->Draw("hist");

      // 2) A bosons canvas
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

      // 3) Four b-quarks (b1–b4)
      TCanvas *c_4b = new TCanvas("c_4b",
         "b_{1}, b_{2}, b_{3}, b_{4} (ordered by p_{T})", 1600, 1200);
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

      // 4) ΔR(b,b) from each A boson
      TCanvas *c_dR_bb = new TCanvas("c_dR_bb",
         "#DeltaR(b,b) from A_{1} and A_{2}", 900, 600);
      c_dR_bb->Divide(2,1);

      c_dR_bb->cd(1);
      if (h_dR_bb_A1) {
         h_dR_bb_A1->SetLineColor(kRed);
         h_dR_bb_A1->SetLineWidth(2);
         h_dR_bb_A1->SetTitle("#DeltaR(b,b) from A_{1};#DeltaR(b,b);Entries");
         h_dR_bb_A1->Draw("hist");
      }

      c_dR_bb->cd(2);
      if (h_dR_bb_A2) {
         h_dR_bb_A2->SetLineColor(kBlue);
         h_dR_bb_A2->SetLineWidth(2);
         h_dR_bb_A2->SetTitle("#DeltaR(b,b) from A_{2};#DeltaR(b,b);Entries");
         h_dR_bb_A2->Draw("hist");
      }

      // DO NOT close file here; histos are already detached in any case
      return; // GEN only
   }

   // =======================
   // MODE 1: RECO CANVASES BEFORE SELECTION
   // =======================
   if (mode == 1) {

      // 1) Multiplicity BEFORE selection
      TCanvas *c_mult_before = new TCanvas("c_mult_before",
                                    "Multiplicity: muons, electrons, jets (before sel)", 1200, 400);
      c_mult_before->Divide(3,1);

      c_mult_before->cd(1);
      if (h_nMu) {
         h_nMu->SetLineColor(kRed);
         h_nMu->SetLineWidth(2);
         h_nMu->SetTitle("Muon multiplicity (before sel);N_{#mu};Events");
         h_nMu->Draw("hist");
      }

      c_mult_before->cd(2);
      if (h_nEle) {
         h_nEle->SetLineColor(kBlue);
         h_nEle->SetLineWidth(2);
         h_nEle->SetTitle("Electron multiplicity (before sel);N_{e};Events");
         h_nEle->Draw("hist");
      }

      c_mult_before->cd(3);
      if (h_nJet) {
         h_nJet->SetLineColor(kGreen+2);
         h_nJet->SetLineWidth(2);
         h_nJet->SetTitle("Jet multiplicity (before sel);N_{jets};Events");
         h_nJet->Draw("hist");
      }

      // 2) Electrons BEFORE selection
      TCanvas *c_eleReco_before = new TCanvas("c_eleReco_before",
         "Leading / subleading electrons (before sel)", 1200, 800);
      c_eleReco_before->Divide(3,2);

      c_eleReco_before->cd(1);
      if (h_pt_e1) {
         h_pt_e1->SetLineColor(kRed);
         h_pt_e1->SetLineWidth(2);
         h_pt_e1->SetTitle("e_{1} p_{T};p_{T} [GeV];Entries");
         h_pt_e1->Draw("hist");
      }

      c_eleReco_before->cd(2);
      if (h_eta_e1) {
         h_eta_e1->SetLineColor(kRed);
         h_eta_e1->SetLineWidth(2);
         h_eta_e1->SetTitle("e_{1} #eta;#eta;Entries");
         h_eta_e1->Draw("hist");
      }

      c_eleReco_before->cd(3);
      if (h_phi_e1) {
         h_phi_e1->SetLineColor(kRed);
         h_phi_e1->SetLineWidth(2);
         h_phi_e1->SetTitle("e_{1} #phi;#phi;Entries");
         h_phi_e1->Draw("hist");
      }

      c_eleReco_before->cd(4);
      if (h_pt_e2) {
         h_pt_e2->SetLineColor(kBlue);
         h_pt_e2->SetLineWidth(2);
         h_pt_e2->SetTitle("e_{2} p_{T};p_{T} [GeV];Entries");
         h_pt_e2->Draw("hist");
      }

      c_eleReco_before->cd(5);
      if (h_eta_e2) {
         h_eta_e2->SetLineColor(kBlue);
         h_eta_e2->SetLineWidth(2);
         h_eta_e2->SetTitle("e_{2} #eta;#eta;Entries");
         h_eta_e2->Draw("hist");
      }

      c_eleReco_before->cd(6);
      if (h_phi_e2) {
         h_phi_e2->SetLineColor(kBlue);
         h_phi_e2->SetLineWidth(2);
         h_phi_e2->SetTitle("e_{2} #phi;#phi;Entries");
         h_phi_e2->Draw("hist");
      }

      // 3) Muons BEFORE selection
      TCanvas *c_muReco_before = new TCanvas("c_muReco_before",
         "Leading / subleading muons (before sel)", 1200, 800);
      c_muReco_before->Divide(3,2);

      c_muReco_before->cd(1);
      if (h_pt_mu1) {
         h_pt_mu1->SetLineColor(kRed);
         h_pt_mu1->SetLineWidth(2);
         h_pt_mu1->SetTitle("#mu_{1} p_{T};p_{T} [GeV];Entries");
         h_pt_mu1->Draw("hist");
      }

      c_muReco_before->cd(2);
      if (h_eta_mu1) {
         h_eta_mu1->SetLineColor(kRed);
         h_eta_mu1->SetLineWidth(2);
         h_eta_mu1->SetTitle("#mu_{1} #eta;#eta;Entries");
         h_eta_mu1->Draw("hist");
      }

      c_muReco_before->cd(3);
      if (h_phi_mu1) {
         h_phi_mu1->SetLineColor(kRed);
         h_phi_mu1->SetLineWidth(2);
         h_phi_mu1->SetTitle("#mu_{1} #phi;#phi;Entries");
         h_phi_mu1->Draw("hist");
      }

      c_muReco_before->cd(4);
      if (h_pt_mu2) {
         h_pt_mu2->SetLineColor(kBlue);
         h_pt_mu2->SetLineWidth(2);
         h_pt_mu2->SetTitle("#mu_{2} p_{T};p_{T} [GeV];Entries");
         h_pt_mu2->Draw("hist");
      }

      c_muReco_before->cd(5);
      if (h_eta_mu2) {
         h_eta_mu2->SetLineColor(kBlue);
         h_eta_mu2->SetLineWidth(2);
         h_eta_mu2->SetTitle("#mu_{2} #eta;#eta;Entries");
         h_eta_mu2->Draw("hist");
      }

      c_muReco_before->cd(6);
      if (h_phi_mu2) {
         h_phi_mu2->SetLineColor(kBlue);
         h_phi_mu2->SetLineWidth(2);
         h_phi_mu2->SetTitle("#mu_{2} #phi;#phi;Entries");
         h_phi_mu2->Draw("hist");
      }

      // 4) Jets J1..J4 BEFORE selection
      TCanvas *c_jetReco_before = new TCanvas("c_jetReco_before",
         "Jets J_{1}..J_{4} (before sel)", 1200, 1000);
      c_jetReco_before->Divide(3,4);

      c_jetReco_before->cd(1);
      if (h_pt_J1) {
         h_pt_J1->SetLineColor(kRed);
         h_pt_J1->SetLineWidth(2);
         h_pt_J1->SetTitle("J_{1} p_{T};p_{T} [GeV];Entries");
         h_pt_J1->Draw("hist");
      }

      c_jetReco_before->cd(2);
      if (h_eta_J1) {
         h_eta_J1->SetLineColor(kRed);
         h_eta_J1->SetLineWidth(2);
         h_eta_J1->SetTitle("J_{1} #eta;#eta;Entries");
         h_eta_J1->Draw("hist");
      }

      c_jetReco_before->cd(3);
      if (h_phi_J1) {
         h_phi_J1->SetLineColor(kRed);
         h_phi_J1->SetLineWidth(2);
         h_phi_J1->SetTitle("J_{1} #phi;#phi;Entries");
         h_phi_J1->Draw("hist");
      }

      c_jetReco_before->cd(4);
      if (h_pt_J2) {
         h_pt_J2->SetLineColor(kBlue);
         h_pt_J2->SetLineWidth(2);
         h_pt_J2->SetTitle("J_{2} p_{T};p_{T} [GeV];Entries");
         h_pt_J2->Draw("hist");
      }

      c_jetReco_before->cd(5);
      if (h_eta_J2) {
         h_eta_J2->SetLineColor(kBlue);
         h_eta_J2->SetLineWidth(2);
         h_eta_J2->SetTitle("J_{2} #eta;#eta;Entries");
         h_eta_J2->Draw("hist");
      }

      c_jetReco_before->cd(6);
      if (h_phi_J2) {
         h_phi_J2->SetLineColor(kBlue);
         h_phi_J2->SetLineWidth(2);
         h_phi_J2->SetTitle("J_{2} #phi;#phi;Entries");
         h_phi_J2->Draw("hist");
      }

      c_jetReco_before->cd(7);
      if (h_pt_J3) {
         h_pt_J3->SetLineColor(kGreen+2);
         h_pt_J3->SetLineWidth(2);
         h_pt_J3->SetTitle("J_{3} p_{T};p_{T} [GeV];Entries");
         h_pt_J3->Draw("hist");
      }

      c_jetReco_before->cd(8);
      if (h_eta_J3) {
         h_eta_J3->SetLineColor(kGreen+2);
         h_eta_J3->SetLineWidth(2);
         h_eta_J3->SetTitle("J_{3} #eta;#eta;Entries");
         h_eta_J3->Draw("hist");
      }

      c_jetReco_before->cd(9);
      if (h_phi_J3) {
         h_phi_J3->SetLineColor(kGreen+2);
         h_phi_J3->SetLineWidth(2);
         h_phi_J3->SetTitle("J_{3} #phi;#phi;Entries");
         h_phi_J3->Draw("hist");
      }

      c_jetReco_before->cd(10);
      if (h_pt_J4) {
         h_pt_J4->SetLineColor(kMagenta);
         h_pt_J4->SetLineWidth(2);
         h_pt_J4->SetTitle("J_{4} p_{T};p_{T} [GeV];Entries");
         h_pt_J4->Draw("hist");
      }

      c_jetReco_before->cd(11);
      if (h_eta_J4) {
         h_eta_J4->SetLineColor(kMagenta);
         h_eta_J4->SetLineWidth(2);
         h_eta_J4->SetTitle("J_{4} #eta;#eta;Entries");
         h_eta_J4->Draw("hist");
      }

      c_jetReco_before->cd(12);
      if (h_phi_J4) {
         h_phi_J4->SetLineColor(kMagenta);
         h_phi_J4->SetLineWidth(2);
         h_phi_J4->SetTitle("J_{4} #phi;#phi;Entries");
         h_phi_J4->Draw("hist");
      }

      // 5) MET & MET φ BEFORE selection
      TCanvas *c_MET_before = new TCanvas("c_MET_before",
         "Puppi MET (before sel)", 800, 400);
      c_MET_before->Divide(2,1);

      c_MET_before->cd(1);
      if (h_MET) {
         h_MET->SetLineColor(kBlack);
         h_MET->SetLineWidth(2);
         h_MET->SetTitle("Puppi MET (before sel);p_{T}^{miss} [GeV];Events");
         h_MET->Draw("hist");
      }

      c_MET_before->cd(2);
      if (h_MET_phi) {
         h_MET_phi->SetLineColor(kBlack);
         h_MET_phi->SetLineWidth(2);
         h_MET_phi->SetTitle("Puppi MET #phi (before sel);#phi^{miss};Events");
         h_MET_phi->Draw("hist");
      }

      // 6) ΔR(J1,ℓ1) before & after cleaning
      TCanvas *c_dR_J1_lep = new TCanvas("c_dR_J1_lep",
         "#DeltaR(J_{1}, #ell_{1}) before/after cleaning", 1000, 800);
      c_dR_J1_lep->Divide(2,2);

      // (1) ΔR(J1,e1) BEFORE cleaning
      c_dR_J1_lep->cd(1);
      if (h_dR_J1_e1) {
         h_dR_J1_e1->SetLineColor(kRed);
         h_dR_J1_e1->SetLineWidth(2);
         h_dR_J1_e1->SetTitle("#DeltaR(J_{1},e_{1}) BEFORE cleaning;#DeltaR(J_{1},e_{1});Events");
         h_dR_J1_e1->Draw("hist");
      }

      // (2) ΔR(J1,μ1) BEFORE cleaning
      c_dR_J1_lep->cd(2);
      if (h_dR_J1_mu1) {
         h_dR_J1_mu1->SetLineColor(kBlue);
         h_dR_J1_mu1->SetLineWidth(2);
         h_dR_J1_mu1->SetTitle("#DeltaR(J_{1},#mu_{1}) BEFORE cleaning;#DeltaR(J_{1},#mu_{1});Events");
         h_dR_J1_mu1->Draw("hist");
      }

      // (3) ΔR(J1_clean,e1) AFTER cleaning
      c_dR_J1_lep->cd(3);
      if (h_dR_J1_e1_clean) {
         h_dR_J1_e1_clean->SetLineColor(kRed+2);
         h_dR_J1_e1_clean->SetLineWidth(2);
         h_dR_J1_e1_clean->SetTitle("#DeltaR(J_{1}^{clean},e_{1});#DeltaR(J_{1}^{clean},e_{1});Events");
         h_dR_J1_e1_clean->Draw("hist");
      }

      // (4) ΔR(J1_clean,μ1) AFTER cleaning
      c_dR_J1_lep->cd(4);
      if (h_dR_J1_mu1_clean) {
         h_dR_J1_mu1_clean->SetLineColor(kBlue+2);
         h_dR_J1_mu1_clean->SetLineWidth(2);
         h_dR_J1_mu1_clean->SetTitle("#DeltaR(J_{1}^{clean},#mu_{1});#DeltaR(J_{1}^{clean},#mu_{1});Events");
         h_dR_J1_mu1_clean->Draw("hist");
      }

      return; // RECO before selection only
   }

   // =======================
   // MODE 2: RECO CANVASES AFTER ALL CUTS
   // =======================
   if (mode == 2) {

      // 1) |Δφ(J1,J2)| AFTER ALL CUTS
      TCanvas *c_dphi = new TCanvas("c_dphi_J1J2",
         "|#Delta#phi(J_{1},J_{2})| after all cuts", 600, 500);
      if (h_dphi_J1J2) {
         h_dphi_J1J2->SetLineWidth(2);
         h_dphi_J1J2->SetTitle("|#Delta#phi(J_{1},J_{2})| (after all cuts);|#Delta#phi|;Events");
         h_dphi_J1J2->Draw("hist");
      }

      // 2) Reconstructed Higgs (2b system) AFTER ALL CUTS
      TCanvas *c_Hreco = new TCanvas("c_Hreco",
         "Reconstructed Higgs from two b-tagged jets (after all cuts)", 1200, 400);
      c_Hreco->Divide(3,1);

      c_Hreco->cd(1);
      if (h_m_2b) {
         h_m_2b->SetLineWidth(2);
         h_m_2b->SetTitle("m_{bb};m_{bb} [GeV];Events");
         h_m_2b->Draw("hist");
      }

      c_Hreco->cd(2);
      if (h_pt_2b) {
         h_pt_2b->SetLineWidth(2);
         h_pt_2b->SetTitle("p_{T}^{bb};p_{T}^{bb} [GeV];Events");
         h_pt_2b->Draw("hist");
      }

      c_Hreco->cd(3);
      if (h_eta_2b) {
         h_eta_2b->SetLineWidth(2);
         h_eta_2b->SetTitle("#eta^{bb};#eta^{bb};Events");
         h_eta_2b->Draw("hist");
      }

      // 3) MET AFTER ALL CUTS
      TCanvas *c_MET_after = new TCanvas("c_MET_after",
         "Puppi MET (after all cuts)", 600, 500);
      if (h_MET_after) {
         h_MET_after->SetLineColor(kBlack);
         h_MET_after->SetLineWidth(2);
         h_MET_after->SetTitle("Puppi MET (after all cuts);p_{T}^{miss} [GeV];Events");
         h_MET_after->Draw("hist");
      }

      // 4) b-tagged jets (leading & subleading) AFTER ALL CUTS
      TCanvas *c_bjets = new TCanvas("c_bjets_after",
         "Leading / subleading b-tagged jets (after all cuts)", 1200, 800);
      c_bjets->Divide(3,2);

      // Leading b-tagged jet
      c_bjets->cd(1);
      if (h_pt_bjet1) {
         h_pt_bjet1->SetLineColor(kRed);
         h_pt_bjet1->SetLineWidth(2);
         h_pt_bjet1->SetTitle("Leading b-jet p_{T};p_{T} [GeV];Events");
         h_pt_bjet1->Draw("hist");
      }

      c_bjets->cd(2);
      if (h_eta_bjet1) {
         h_eta_bjet1->SetLineColor(kRed);
         h_eta_bjet1->SetLineWidth(2);
         h_eta_bjet1->SetTitle("Leading b-jet #eta;#eta;Events");
         h_eta_bjet1->Draw("hist");
      }

      c_bjets->cd(3);
      if (h_m_bjet1) {
         h_m_bjet1->SetLineColor(kRed);
         h_m_bjet1->SetLineWidth(2);
         h_m_bjet1->SetTitle("Leading b-jet mass;m_{j} [GeV];Events");
         h_m_bjet1->Draw("hist");
      }

      // Subleading b-tagged jet
      c_bjets->cd(4);
      if (h_pt_bjet2) {
         h_pt_bjet2->SetLineColor(kBlue);
         h_pt_bjet2->SetLineWidth(2);
         h_pt_bjet2->SetTitle("Subleading b-jet p_{T};p_{T} [GeV];Events");
         h_pt_bjet2->Draw("hist");
      }

      c_bjets->cd(5);
      if (h_eta_bjet2) {
         h_eta_bjet2->SetLineColor(kBlue);
         h_eta_bjet2->SetLineWidth(2);
         h_eta_bjet2->SetTitle("Subleading b-jet #eta;#eta;Events");
         h_eta_bjet2->Draw("hist");
      }

      c_bjets->cd(6);
      if (h_m_bjet2) {
         h_m_bjet2->SetLineColor(kBlue);
         h_m_bjet2->SetLineWidth(2);
         h_m_bjet2->SetTitle("Subleading b-jet mass;m_{j} [GeV];Events");
         h_m_bjet2->Draw("hist");
      }

      // 5) Event hardness: HT, HT_2b, ST
      TCanvas *c_HT = new TCanvas("c_HT",
         "Event hardness variables (after all cuts)", 1200, 400);
      c_HT->Divide(3,1);

      c_HT->cd(1);
      if (h_HT) {
         h_HT->SetLineWidth(2);
         h_HT->SetTitle("H_{T} = #Sigma p_{T}^{jets};H_{T} [GeV];Events");
         h_HT->Draw("hist");
      }

      c_HT->cd(2);
      if (h_HT_2b) {
         h_HT_2b->SetLineWidth(2);
         h_HT_2b->SetTitle("H_{T}^{2b} = p_{T}(b_{1})+p_{T}(b_{2});H_{T}^{2b} [GeV];Events");
         h_HT_2b->Draw("hist");
      }

      c_HT->cd(3);
      if (h_ST) {
         h_ST->SetLineWidth(2);
         h_ST->SetTitle("S_{T} = H_{T} + p_{T}^{miss};S_{T} [GeV];Events");
         h_ST->Draw("hist");
      }

      // 6) b-tag multiplicity after all cuts
      TCanvas *c_Nb = new TCanvas("c_Nbjets_after",
         "b-tagged jet multiplicity (after all cuts)", 600, 500);
      if (h_Nbjets_after) {
         h_Nbjets_after->SetLineWidth(2);
         h_Nbjets_after->SetTitle("N_{b-jets} (after all cuts);N_{b-jets};Events");
         h_Nbjets_after->Draw("hist");
      }

      // 7) MET-related angles after all cuts
      TCanvas *c_METangles = new TCanvas("c_MET_angles",
         "MET-related angles (after all cuts)", 800, 400);
      c_METangles->Divide(2,1);

      c_METangles->cd(1);
      if (h_dphi_MET_bb) {
         h_dphi_MET_bb->SetLineWidth(2);
         h_dphi_MET_bb->SetTitle("|#Delta#phi(MET,bb)| (after all cuts);|#Delta#phi(MET,bb)|;Events");
         h_dphi_MET_bb->Draw("hist");
      }

      c_METangles->cd(2);
      if (h_dphi_MET_J1) {
         h_dphi_MET_J1->SetLineWidth(2);
         h_dphi_MET_J1->SetTitle("|#Delta#phi(MET,J_{1})| (after all cuts);|#Delta#phi(MET,J_{1})|;Events");
         h_dphi_MET_J1->Draw("hist");
      }

      return; // RECO after all cuts only
   }

   std::cout << "DrawggHAnalysis: unknown mode " << mode
             << " (use 0 for GEN, 1 for RECO-before, 2 for RECO-after)" << std::endl;
}
