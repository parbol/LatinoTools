#define SkimMeBaby_cxx
#include "SkimMeBaby.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include "TLorentzVector.h" 
#define MMUON 105.658369e-3
#define MELECTRON 0.5109989461e-3
#include "mt2_bisect.h"

Int_t maxEntries = -1;

void SkimMeBaby::Loop()
{

   if (fChain == 0) return;
  
   //Opening output file
   TFile* out = TFile::Open(outputSample, "recreate");
   fChain->LoadTree(0); 
   outputtree = (TTree*) fChain->GetTree()->CloneTree(0);

   //Variables only present in MC
   if (!outputSample.Contains("Run2017")) {
     fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
     outputtree->SetBranchStatus("Generator_weight", 1);
     fChain->SetBranchAddress("baseW", &baseW, &b_baseW);
     outputtree->SetBranchStatus("baseW", 1);
     fChain->SetBranchAddress("Xsec", &Xsec, &b_Xsec);
     outputtree->SetBranchStatus("Xsec", 1);
   }  

   disconnect(); 
   outputtree->SetBranchStatus("run", 1);
   outputtree->SetBranchStatus("event", 1);
   outputtree->SetBranchStatus("Electron_dxy", 1);
   outputtree->SetBranchStatus("Electron_dz", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP80", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP90", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17Iso_WPL", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP80", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP90", 1);
   outputtree->SetBranchStatus("Electron_mvaFall17noIso_WPL", 1);
   outputtree->SetBranchStatus("MET_phi", 1);
   outputtree->SetBranchStatus("MET_pt", 1);
   outputtree->SetBranchStatus("Flag_HBHENoiseFilter", 1);
   outputtree->SetBranchStatus("Flag_HBHENoiseIsoFilter", 1);
   outputtree->SetBranchStatus("Flag_CSCTightHaloFilter", 1);
   outputtree->SetBranchStatus("Flag_CSCTightHaloTrkMuUnvetoFilter", 1);
   outputtree->SetBranchStatus("Flag_CSCTightHalo2015Filter", 1);
   outputtree->SetBranchStatus("Flag_globalTightHalo2016Filter", 1);
   outputtree->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1);
   outputtree->SetBranchStatus("Flag_HcalStripHaloFilter", 1);
   outputtree->SetBranchStatus("Flag_hcalLaserEventFilter", 1);
   outputtree->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
   outputtree->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 1);
   outputtree->SetBranchStatus("Flag_ecalBadCalibFilter", 1);
   outputtree->SetBranchStatus("Flag_goodVertices", 1);
   outputtree->SetBranchStatus("Flag_eeBadScFilter", 1);
   outputtree->SetBranchStatus("Flag_ecalLaserCorrFilter", 1);
   outputtree->SetBranchStatus("Flag_trkPOGFilters", 1);
   outputtree->SetBranchStatus("Flag_chargedHadronTrackResolutionFilter", 1);
   outputtree->SetBranchStatus("Flag_muonBadTrackFilter", 1);
   outputtree->SetBranchStatus("Flag_BadChargedCandidateFilter", 1);
   outputtree->SetBranchStatus("Flag_BadPFMuonFilter", 1);
   outputtree->SetBranchStatus("Flag_BadChargedCandidateSummer16Filter", 1);
   outputtree->SetBranchStatus("Flag_BadPFMuonSummer16Filter", 1);
   outputtree->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1);
   outputtree->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1);
   outputtree->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1);
   outputtree->SetBranchStatus("Flag_METFilters", 1);
   outputtree->SetBranchStatus("nLepton", 1);
   outputtree->SetBranchStatus("Lepton_pdgId", 1);
   outputtree->SetBranchStatus("Lepton_electronIdx", 1);
   outputtree->SetBranchStatus("Lepton_muonIdx", 1);
   outputtree->SetBranchStatus("Lepton_pt", 1);
   outputtree->SetBranchStatus("Lepton_eta", 1);
   outputtree->SetBranchStatus("Lepton_phi", 1);
   outputtree->SetBranchStatus("Lepton_isTightMuon_cut_Tight80x_HWWW", 1);
   outputtree->SetBranchStatus("puWeight", 1);
   outputtree->SetBranchStatus("mll", 1);
   outputtree->SetBranchStatus("dphill", 1);
   outputtree->SetBranchStatus("ptll", 1);
   outputtree->SetBranchStatus("pt1", 1);
   outputtree->SetBranchStatus("pt2", 1);
   outputtree->SetBranchStatus("mth", 1);
   outputtree->SetBranchStatus("mT2", 1);
   outputtree->SetBranchStatus("channel", 1);
   outputtree->SetBranchStatus("drll", 1);
   outputtree->SetBranchStatus("njet", 1);
   outputtree->SetBranchStatus("eta1", 1);
   outputtree->SetBranchStatus("eta2", 1);
   outputtree->SetBranchStatus("phi1", 1);
   outputtree->SetBranchStatus("phi2", 1);

   Int_t nTight = 0;
   Float_t dxycut = 0.;
   Float_t dzcut = 0.;

   Float_t dphi = -999.;
   Float_t deta = -999.;
   Float_t mt2 = -999.;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   if(maxEntries != -1) nentries = maxEntries;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%10000 == 0) printf("Entry number %d \n", jentry);
   
      if(nLepton < 2) continue;
      std::vector<int> leptonIndex;

      nTight = 0;

      for(Int_t i=0; i<nLepton; i++) {

	//Electrons
         if(abs(Lepton_pdgId[i]) == 11) {
            if(Lepton_pt[i] > 20.0) {
              dxycut = 0.02;
            } else {
              dxycut = 0.01;
            }   
            if(Electron_mvaFall17Iso_WP80[i] > 0.5 && Electron_dz[Lepton_electronIdx[i]]< 0.1 && Electron_dxy[Lepton_electronIdx[i]]< dxycut) {
                nTight ++;
                leptonIndex.push_back(i);
            } 
         }

         //Muons 
         if(abs(Lepton_pdgId[i]) == 13) {
            if(Lepton_isTightMuon_cut_Tight80x_HWWW[i] > 0.5) {
                nTight ++;
                leptonIndex.push_back(i);
            } 
         }
      }
      
      if(nTight < 2) continue;
      if(leptonIndex.size() < 2) continue;
      TLorentzVector lepton1, lepton2;

      pt1 = Lepton_pt[leptonIndex[0]];
      eta1 = Lepton_eta[leptonIndex[0]];
      phi1 = Lepton_phi[leptonIndex[0]];
      pt2 = Lepton_pt[leptonIndex[1]];
      eta2 = Lepton_eta[leptonIndex[1]];
      phi2 = Lepton_phi[leptonIndex[1]];

      if(abs(Lepton_pdgId[leptonIndex[0]]) == 11) {  
          lepton1.SetPtEtaPhiM(pt1, eta1, phi1, MELECTRON);
      } else {
          lepton1.SetPtEtaPhiM(pt1, eta1, phi1, MMUON);
      }
      if(abs(Lepton_pdgId[leptonIndex[1]]) == 11) {  
          lepton2.SetPtEtaPhiM(pt2, eta2, phi2, MELECTRON);
      } else {
          lepton2.SetPtEtaPhiM(pt2, eta2, phi2, MMUON);
      }
      mll = (lepton1 + lepton2).M();
      ptll = (lepton1 + lepton2).Pt();
      dphi = lepton1.DeltaPhi(lepton2);
      deta = fabs(lepton1.Eta() - lepton2.Eta());

      double pa[3]; 
      double pb[3];
      if(abs(Lepton_pdgId[leptonIndex[0]]) == 11) {  
         pa[0] = MELECTRON;
         pa[1] = lepton1.Px();     
         pa[2] = lepton1.Py();
      } else {     
         pa[0] = MMUON;
         pa[1] = lepton1.Px();     
         pa[2] = lepton1.Py();
      }
      if(abs(Lepton_pdgId[leptonIndex[1]]) == 11) {  
         pb[0] = MELECTRON;
         pb[1] = lepton2.Px();     
         pb[2] = lepton2.Py();
      } else {     
         pb[0] = MMUON;
         pb[1] = lepton2.Px();     
         pb[2] = lepton2.Py();
      }

      double pmiss[3] = {0, MET_pt*cos(MET_phi), MET_pt*sin(MET_phi)};
      double mn = 0.0;
      mt2_bisect::mt2 mt2_event;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2 = (float) mt2_event.get_mt2();


      outputtree->Fill();

   }

   outputtree->Write();
   out->Close();
}



