#ifndef PDGHELPERS_C
#define PDGHELPERS_C

#include "ggHAnalysis.h"   // so ggHAnalysis is known
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TRandom.h>
#include <TTree.h>

// PDG codes for the signal
const int PDG_H     = 25;
const int PDG_A     = 36;
const int PDG_B     = 5;
const int PDG_BBAR  = -5;
const int PDG_GLUON = 21;

// ----------------------------------------------------------------------
// PDGName: map PDG code to a readable string
// ----------------------------------------------------------------------
std::string PDGName(int pdg)
{
   switch (pdg)
   {
      // Fundamental particles
      case  21: return "g";        // gluon
      case  25: return "H";        // Higgs
      case  36: return "A";        // CP-odd Higgs
      case   1: return "d";
      case  -1: return "dbar";
      case   2: return "u";
      case  -2: return "ubar";
      case   5: return "b";
      case  -5: return "bbar";

      // Leptons
      case  11: return "e-";
      case -11: return "e+";
      case  13: return "mu-";
      case -13: return "mu+";
      case  15: return "tau-";
      case -15: return "tau+";

      // B mesons
      case  511: return "B0";
      case -511: return "B0bar";
      case  521: return "B+";
      case -521: return "B-";
      case  523: return "B*+";
      case -523: return "B*-";

      // D mesons
      case  411: return "D+";
      case -411: return "D-";
      case  421: return "D0";
      case -421: return "D0bar";
      case  413: return "D*+";
      case -413: return "D*-";
      case  423: return "D*0";
      case -423: return "D*0bar";

      default:
         return "pdg" + std::to_string(pdg);
   }
}

// ============================================================================
// Dump GEN info for one event into a text file (optional helper)
// ============================================================================
void PrintGenInfo(ggHAnalysis *obj)
{
   if (!obj || !obj->fChain) {
      std::cout << "PrintGenInfo: invalid ggHAnalysis object or fChain" << std::endl;
      return;
   }

   TTree *chain = obj->fChain;
   Long64_t nentries = chain->GetEntriesFast();
   if (nentries == 0) {
      std::cout << "PrintGenInfo: no events in the file" << std::endl;
      return;
   }

   const int NPRINT = 10;

   for (int k = 0; k < NPRINT; ++k) {

      Long64_t evt = (Long64_t)(gRandom->Rndm() * nentries);
      chain->GetEntry(evt);

      std::string fname = "GenInfo_event" + std::to_string(evt) + ".txt";
      std::ofstream out(fname.c_str());

      if (!out.is_open()) {
         std::cout << "PrintGenInfo: could not open " << fname << std::endl;
         continue;
      }

      out << "GEN particle dump for RANDOM event " << evt << "\n\n";

      out << std::left
          << std::setw(5)  << "Idx"
          << std::setw(10) << "Name"
          << std::setw(8)  << "PDG"
          << std::setw(8)  << "Status"
          << std::setw(8)  << "MomIdx"
          << std::setw(12) << "MomName"
          << std::setw(12) << "pT[GeV]"
          << std::setw(12) << "Mass[GeV]"
          << "\n";

      out << std::string(5+10+8+8+8+12+12+12, '-') << "\n";

      for (int i = 0; i < obj->nGenPart; ++i) {

         int pdg    = obj->GenPart_pdgId[i];
         int status = obj->GenPart_status[i];
         int momIdx = obj->GenPart_genPartIdxMother[i];

         std::string momName = "none";
         if (momIdx >= 0 && momIdx < obj->nGenPart) {
            momName = PDGName(obj->GenPart_pdgId[momIdx]);
         }

         out << std::left
             << std::setw(5)  << i
             << std::setw(10) << PDGName(pdg)
             << std::setw(8)  << pdg
             << std::setw(8)  << status
             << std::setw(8)  << momIdx
             << std::setw(12) << momName
             << std::setw(12) << obj->GenPart_pt[i]
             << std::setw(12) << obj->GenPart_mass[i]
             << "\n";
      }

      out.close();

      std::cout << "Printed random event " << evt
                << " into file " << fname << std::endl;
   }
}

#endif // PDGHELPERS_C
