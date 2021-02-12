#include <iostream>
#include <algorithm>

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>


// This test file does a number of basic checks on the Hits files.
// Clearly, a SPECIFIC structure of the file is to be assumed

int main(int argc, char **argv)
{
   std::cout << "*********************************\n";
   std::cout << "     " << argv[0] << "       \n"; 
   std::cout << "*********************************\n";
   std::cout << "\n\n" << std::endl;
   
   if (argc != 2) {
     std::cout << "USAGE: " << argv[0] << " [INPUTFILE] \n";
     std::cout << "where [INPUTFILE] is the podio file containing simulation hits\n";
     return -1;
   }

   podio::ROOTReader l_reader;
   l_reader.openFile(argv[1]);
   if (! l_reader.isValid()){
     std::cout << "Error: cannot open file " << argv[1] << std::endl;
     return -1;
   }

   podio::EventStore l_store;
   l_store.setReader(&l_reader);
   
   // Find the number of events

   unsigned int nevents  = l_reader.getEntries();
   std::cout << "The file contains " << nevents << " events " << std::endl;

   // Now opening the output file

   TFile f_output("sim_output_validation.hist.root","recreate");

   std::vector<Float_t> s_avg_hits;
   std::vector<Float_t> c_avg_hits;

   unsigned int evt_avg = 0;
   unsigned int s_nhits_avg = 0;
   unsigned int c_nhits_avg = 0;

   float s_energy;
   float c_energy;

   TH1F * h_s_energy = new TH1F("h_s_energy","",100,1,-1);
   TH1F * h_c_energy = new TH1F("h_c_energy","",100,1,-1);
   
   // Now looping on all events

   for (unsigned int i_evt = 0; i_evt < nevents; ++i_evt){
     if (i_evt%100 == 0) std::cout << "Processed " << i_evt << " events" << std::endl;
     
     // Trying to access the scintillation hit collections
     auto & s_hitColl = l_store.get<edm4hep::SimCalorimeterHitCollection>("S_caloHits");
     auto & c_hitColl = l_store.get<edm4hep::SimCalorimeterHitCollection>("C_caloHits");
     auto & aux_hitColl = l_store.get<edm4hep::SimCalorimeterHitCollection>("Auxiliary_infoHits");

     // Computing average number of hits in 10 events

     ++evt_avg;
     s_nhits_avg += s_hitColl.size();
     c_nhits_avg += c_hitColl.size();

     if (i_evt%10 == 0 && i_evt != 0){
       s_avg_hits.push_back(float(s_nhits_avg)/float(evt_avg));
       c_avg_hits.push_back(float(c_nhits_avg)/float(evt_avg));
       evt_avg = 0;
       s_nhits_avg = 0;
       c_nhits_avg = 0;
     }

     // compute the total energy and plot it

     s_energy = 0;
     c_energy = 0;

     for (auto&  hit : s_hitColl){
       s_energy += hit.getEnergy();
     }

     for (auto&  hit : c_hitColl){
       c_energy += hit.getEnergy();
     }

     h_s_energy->Fill(s_energy);
     h_c_energy->Fill(c_energy);

     l_store.clear();
     l_reader.endOfEvent();
     
   }

   // Now prepare the TGraphs. First, understand the range of the x and y axis, then produce the plot.

   std::vector<float> myx_s;
   std::vector<float> myx_c;

   for (unsigned int i = 0; i < s_avg_hits.size(); ++i){
     myx_s.push_back(10. * float(i));
   }

   for (unsigned int i = 0; i < c_avg_hits.size(); ++i){
     myx_c.push_back(10. * float(i));
   }

   TGraph * g_s = new TGraph(s_avg_hits.size(),&myx_s[0],&s_avg_hits[0]);
   TGraph * g_c = new TGraph(c_avg_hits.size(),&myx_c[0],&c_avg_hits[0]);

   g_s->SetName("g_hits_stability_s");
   g_c->SetName("g_hits_stability_c");

   TH1F * h_s_hits = g_s->GetHistogram();
   TH1F * h_c_hits = g_c->GetHistogram();

   
   
   // Now need to make the actual plot of the stability of the number of hits vs event number

   f_output.cd();
   h_s_energy->Write();
   h_c_energy->Write();
   g_s->Write();
   g_c->Write();
   f_output.Close();

   l_reader.closeFile();

  return 0;
}
