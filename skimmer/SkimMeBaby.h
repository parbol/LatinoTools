//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  9 19:12:43 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano//Fall2017_nAOD_v1_Full2017/MCl1loose2017__MCformulas/nanoLatino_WLLJJToLNu_M-50_QCD_3Jet__part43.root
//////////////////////////////////////////////////////////

#ifndef SkimMeBaby_h
#define SkimMeBaby_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class SkimMeBaby {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *outputtree;
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         Electron_dxy[100];   //[nElectron]
   Float_t         Electron_dz[100];   //[nElectron]
   Float_t         Generator_weight;
   Bool_t          Electron_mvaFall17Iso_WP80[100];   //[nElectron]
   Bool_t          Electron_mvaFall17Iso_WP90[100];   //[nElectron]
   Bool_t          Electron_mvaFall17Iso_WPL[100];   //[nElectron]
   Bool_t          Electron_mvaFall17noIso_WP80[100];   //[nElectron]
   Bool_t          Electron_mvaFall17noIso_WP90[100];   //[nElectron]
   Bool_t          Electron_mvaFall17noIso_WPL[100];   //[nElectron]
   Float_t         MET_phi;
   Float_t         MET_pt;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   UInt_t          nLepton;
   Int_t           Lepton_pdgId[100];   //[nLepton]
   Int_t           Lepton_electronIdx[100];   //[nLepton]
   Int_t           Lepton_muonIdx[100];   //[nLepton]
   Float_t         Lepton_pt[100];   //[nLepton]
   Float_t         Lepton_eta[100];   //[nLepton]
   Float_t         Lepton_phi[100];   //[nLepton]
   Int_t           Lepton_isTightMuon_cut_Tight80x_HWWW[100];   //[nLepton]
   Float_t         puWeight;
   Float_t         mll;
   Float_t         dphill;
   Float_t         ptll;
   Float_t         pt1;
   Float_t         pt2;
   Float_t         mth;
   Float_t         mT2;
   Float_t         channel;
   Float_t         drll;
   Float_t         njet;
   Float_t         Jet_btagCSVV2[100];
   Float_t         eta1;
   Float_t         eta2;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         baseW;
   Float_t         Xsec;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Mu18_Mu9_SameSign;
   Bool_t          HLT_Mu18_Mu9_SameSign_DZ;
   Bool_t          HLT_Mu18_Mu9;
   Bool_t          HLT_Mu18_Mu9_DZ;
   Bool_t          HLT_Mu20_Mu10_SameSign;
   Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
   Bool_t          HLT_Mu20_Mu10;
   Bool_t          HLT_Mu20_Mu10_DZ;
   Bool_t          HLT_Mu23_Mu12_SameSign;
   Bool_t          HLT_Mu23_Mu12_SameSign_DZ;
   Bool_t          HLT_Mu23_Mu12;
   Bool_t          HLT_Mu23_Mu12_DZ;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Electron_mvaFall17Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17noIso_WPL;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_nLepton;   //!
   TBranch        *b_Lepton_pdgId;   //!
   TBranch        *b_Lepton_electronIdx;   //!
   TBranch        *b_Lepton_muonIdx;   //!
   TBranch        *b_Lepton_pt;   //!
   TBranch        *b_Lepton_eta;   //!
   TBranch        *b_Lepton_phi;   //!
   TBranch        *b_Lepton_isTightMuon_cut_Tight80x_HWWW;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_mll;   //!
   TBranch        *b_dphill;   //!
   TBranch        *b_ptll;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_mth;   //!
   TBranch        *b_mT2;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_drll;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_baseW;   //!
   TBranch        *b_Xsec;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu18_Mu9;   //!
   TBranch        *b_HLT_Mu18_Mu9_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10;   //!
   TBranch        *b_HLT_Mu20_Mu10_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12;   //!
   TBranch        *b_HLT_Mu23_Mu12_DZ;   //!

   SkimMeBaby(TString, TString);
   virtual ~SkimMeBaby();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   TString inputSample, outputSample;
   virtual void     Show(Long64_t entry = -1);
   void disconnect();
};

#endif

#ifdef SkimMeBaby_cxx
SkimMeBaby::SkimMeBaby(TString inputSample_, TString outputSample_) : fChain(0)
{
  inputSample = inputSample_;
  outputSample  = outputSample_;
  TFile* f = TFile::Open(inputSample, "READ");
  TTree* tree = (TTree*) f->Get("Events");
  if(tree == NULL) return;
  Init(tree);
  Loop();
}


SkimMeBaby::~SkimMeBaby()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SkimMeBaby::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SkimMeBaby::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SkimMeBaby::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_mvaFall17Iso_WP80", Electron_mvaFall17Iso_WP80, &b_Electron_mvaFall17Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17Iso_WP90", Electron_mvaFall17Iso_WP90, &b_Electron_mvaFall17Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17Iso_WPL", Electron_mvaFall17Iso_WPL, &b_Electron_mvaFall17Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17noIso_WP80", Electron_mvaFall17noIso_WP80, &b_Electron_mvaFall17noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17noIso_WP90", Electron_mvaFall17noIso_WP90, &b_Electron_mvaFall17noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17noIso_WPL", Electron_mvaFall17noIso_WPL, &b_Electron_mvaFall17noIso_WPL);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("nLepton", &nLepton, &b_nLepton);
   fChain->SetBranchAddress("Lepton_pdgId", Lepton_pdgId, &b_Lepton_pdgId);
   fChain->SetBranchAddress("Lepton_electronIdx", Lepton_electronIdx, &b_Lepton_electronIdx);
   fChain->SetBranchAddress("Lepton_muonIdx", Lepton_muonIdx, &b_Lepton_muonIdx);
   fChain->SetBranchAddress("Lepton_pt", Lepton_pt, &b_Lepton_pt);
   fChain->SetBranchAddress("Lepton_eta", Lepton_eta, &b_Lepton_eta);
   fChain->SetBranchAddress("Lepton_phi", Lepton_phi, &b_Lepton_phi);
   fChain->SetBranchAddress("Lepton_isTightMuon_cut_Tight80x_HWWW", Lepton_isTightMuon_cut_Tight80x_HWWW, &b_Lepton_isTightMuon_cut_Tight80x_HWWW);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
   fChain->SetBranchAddress("dphill", &dphill, &b_dphill);
   fChain->SetBranchAddress("ptll", &ptll, &b_ptll);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("mT2", &mT2, &b_mT2);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("drll", &drll, &b_drll);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ, &b_HLT_Mu18_Mu9_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9, &b_HLT_Mu18_Mu9);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ, &b_HLT_Mu18_Mu9_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign, &b_HLT_Mu23_Mu12_SameSign);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ, &b_HLT_Mu23_Mu12_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12, &b_HLT_Mu23_Mu12);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ, &b_HLT_Mu23_Mu12_DZ);
   
   Notify();
}

Bool_t SkimMeBaby::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SkimMeBaby::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SkimMeBaby::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void SkimMeBaby::disconnect() {

  outputtree->SetBranchStatus("run", 0);
  outputtree->SetBranchStatus("luminosityBlock", 0);
  outputtree->SetBranchStatus("event", 0);
  outputtree->SetBranchStatus("CaloMET_phi", 0);
  outputtree->SetBranchStatus("CaloMET_pt", 0);
  outputtree->SetBranchStatus("CaloMET_sumEt", 0);
  outputtree->SetBranchStatus("nElectron", 0);
  outputtree->SetBranchStatus("Electron_deltaEtaSC", 0);
  outputtree->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 0);
  outputtree->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 0);
  outputtree->SetBranchStatus("Electron_dr03TkSumPt", 0);
  outputtree->SetBranchStatus("Electron_dxy", 0);
  outputtree->SetBranchStatus("Electron_dxyErr", 0);
  outputtree->SetBranchStatus("Electron_dz", 0);
  outputtree->SetBranchStatus("Electron_dzErr", 0);
  outputtree->SetBranchStatus("Electron_eCorr", 0);
  outputtree->SetBranchStatus("Electron_eInvMinusPInv", 0);
  outputtree->SetBranchStatus("Electron_energyErr", 0);
  outputtree->SetBranchStatus("Electron_eta", 0);
  outputtree->SetBranchStatus("Electron_hoe", 0);
  outputtree->SetBranchStatus("Electron_ip3d", 0);
  outputtree->SetBranchStatus("Electron_mass", 0);
  outputtree->SetBranchStatus("Electron_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("Electron_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso", 0);
  outputtree->SetBranchStatus("Electron_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("Electron_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("Electron_phi", 0);
  outputtree->SetBranchStatus("Electron_pt", 0);
  outputtree->SetBranchStatus("Electron_r9", 0);
  outputtree->SetBranchStatus("Electron_sieie", 0);
  outputtree->SetBranchStatus("Electron_sip3d", 0);
  outputtree->SetBranchStatus("Electron_mvaTTH", 0);
  outputtree->SetBranchStatus("Electron_charge", 0);
  outputtree->SetBranchStatus("Electron_cutBased", 0);
  outputtree->SetBranchStatus("Electron_jetIdx", 0);
  outputtree->SetBranchStatus("Electron_pdgId", 0);
  outputtree->SetBranchStatus("Electron_photonIdx", 0);
  outputtree->SetBranchStatus("Electron_tightCharge", 0);
  outputtree->SetBranchStatus("Electron_vidNestedWPBitmap", 0);
  outputtree->SetBranchStatus("Electron_convVeto", 0);
  outputtree->SetBranchStatus("Electron_cutBased_HEEP", 0);
  outputtree->SetBranchStatus("Electron_isPFcand", 0);
  outputtree->SetBranchStatus("Electron_lostHits", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP80", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP90", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WPL", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP80", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP90", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WPL", 0);
  outputtree->SetBranchStatus("nFatJet", 0);
  outputtree->SetBranchStatus("FatJet_area", 0);
  outputtree->SetBranchStatus("FatJet_btagCMVA", 0);
  outputtree->SetBranchStatus("FatJet_btagCSVV2", 0);
  outputtree->SetBranchStatus("FatJet_btagDeepB", 0);
  outputtree->SetBranchStatus("FatJet_btagHbb", 0);
  outputtree->SetBranchStatus("FatJet_eta", 0);
  outputtree->SetBranchStatus("FatJet_mass", 0);
  outputtree->SetBranchStatus("FatJet_msoftdrop", 0);
  outputtree->SetBranchStatus("FatJet_n2b1", 0);
  outputtree->SetBranchStatus("FatJet_n3b1", 0);
  outputtree->SetBranchStatus("FatJet_phi", 0);
  outputtree->SetBranchStatus("FatJet_pt", 0);
  outputtree->SetBranchStatus("FatJet_tau1", 0);
  outputtree->SetBranchStatus("FatJet_tau2", 0);
  outputtree->SetBranchStatus("FatJet_tau3", 0);
  outputtree->SetBranchStatus("FatJet_tau4", 0);
  outputtree->SetBranchStatus("FatJet_jetId", 0);
  outputtree->SetBranchStatus("FatJet_subJetIdx1", 0);
  outputtree->SetBranchStatus("FatJet_subJetIdx2", 0);
  outputtree->SetBranchStatus("nGenJetAK8", 0);
  outputtree->SetBranchStatus("GenJetAK8_eta", 0);
  outputtree->SetBranchStatus("GenJetAK8_mass", 0);
  outputtree->SetBranchStatus("GenJetAK8_phi", 0);
  outputtree->SetBranchStatus("GenJetAK8_pt", 0);
  outputtree->SetBranchStatus("nGenJet", 0);
  outputtree->SetBranchStatus("GenJet_eta", 0);
  outputtree->SetBranchStatus("GenJet_mass", 0);
  outputtree->SetBranchStatus("GenJet_phi", 0);
  outputtree->SetBranchStatus("GenJet_pt", 0);
  outputtree->SetBranchStatus("nGenPart", 0);
  outputtree->SetBranchStatus("GenPart_eta", 0);
  outputtree->SetBranchStatus("GenPart_mass", 0);
  outputtree->SetBranchStatus("GenPart_phi", 0);
  outputtree->SetBranchStatus("GenPart_pt", 0);
  outputtree->SetBranchStatus("GenPart_genPartIdxMother", 0);
  outputtree->SetBranchStatus("GenPart_pdgId", 0);
  outputtree->SetBranchStatus("GenPart_status", 0);
  outputtree->SetBranchStatus("GenPart_statusFlags", 0);
  outputtree->SetBranchStatus("nSubGenJetAK8", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_eta", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_mass", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_phi", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_pt", 0);
  outputtree->SetBranchStatus("Generator_binvar", 0);
  outputtree->SetBranchStatus("Generator_scalePDF", 0);
  outputtree->SetBranchStatus("Generator_weight", 0);
  outputtree->SetBranchStatus("Generator_x1", 0);
  outputtree->SetBranchStatus("Generator_x2", 0);
  outputtree->SetBranchStatus("Generator_xpdf1", 0);
  outputtree->SetBranchStatus("Generator_xpdf2", 0);
  outputtree->SetBranchStatus("Generator_id1", 0);
  outputtree->SetBranchStatus("Generator_id2", 0);
  outputtree->SetBranchStatus("nGenVisTau", 0);
  outputtree->SetBranchStatus("GenVisTau_eta", 0);
  outputtree->SetBranchStatus("GenVisTau_mass", 0);
  outputtree->SetBranchStatus("GenVisTau_phi", 0);
  outputtree->SetBranchStatus("GenVisTau_pt", 0);
  outputtree->SetBranchStatus("GenVisTau_charge", 0);
  outputtree->SetBranchStatus("GenVisTau_genPartIdxMother", 0);
  outputtree->SetBranchStatus("GenVisTau_status", 0);
  outputtree->SetBranchStatus("genWeight", 0);
  outputtree->SetBranchStatus("LHEWeight_originalXWGTUP", 0);
  outputtree->SetBranchStatus("nLHEPdfWeight", 0);
  outputtree->SetBranchStatus("LHEPdfWeight", 0);
  outputtree->SetBranchStatus("nLHEScaleWeight", 0);
  outputtree->SetBranchStatus("LHEScaleWeight", 0);
  outputtree->SetBranchStatus("nIsoTrack", 0);
  outputtree->SetBranchStatus("IsoTrack_dxy", 0);
  outputtree->SetBranchStatus("IsoTrack_dz", 0);
  outputtree->SetBranchStatus("IsoTrack_eta", 0);
  outputtree->SetBranchStatus("IsoTrack_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("IsoTrack_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("IsoTrack_phi", 0);
  outputtree->SetBranchStatus("IsoTrack_pt", 0);
  outputtree->SetBranchStatus("IsoTrack_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("IsoTrack_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("IsoTrack_pdgId", 0);
  outputtree->SetBranchStatus("IsoTrack_isHighPurityTrack", 0);
  outputtree->SetBranchStatus("IsoTrack_isPFcand", 0);
  outputtree->SetBranchStatus("nJet", 0);
  outputtree->SetBranchStatus("Jet_area", 0);
  outputtree->SetBranchStatus("Jet_btagCMVA", 0);
  outputtree->SetBranchStatus("Jet_btagCSVV2", 0);
  outputtree->SetBranchStatus("Jet_btagDeepB", 0);
  outputtree->SetBranchStatus("Jet_btagDeepC", 0);
  outputtree->SetBranchStatus("Jet_btagDeepFlavB", 0);
  outputtree->SetBranchStatus("Jet_chEmEF", 0);
  outputtree->SetBranchStatus("Jet_chHEF", 0);
  outputtree->SetBranchStatus("Jet_eta", 0);
  outputtree->SetBranchStatus("Jet_mass", 0);
  outputtree->SetBranchStatus("Jet_neEmEF", 0);
  outputtree->SetBranchStatus("Jet_neHEF", 0);
  outputtree->SetBranchStatus("Jet_phi", 0);
  outputtree->SetBranchStatus("Jet_pt", 0);
  outputtree->SetBranchStatus("Jet_qgl", 0);
  outputtree->SetBranchStatus("Jet_rawFactor", 0);
  outputtree->SetBranchStatus("Jet_bReg", 0);
  outputtree->SetBranchStatus("Jet_electronIdx1", 0);
  outputtree->SetBranchStatus("Jet_electronIdx2", 0);
  outputtree->SetBranchStatus("Jet_jetId", 0);
  outputtree->SetBranchStatus("Jet_muonIdx1", 0);
  outputtree->SetBranchStatus("Jet_muonIdx2", 0);
  outputtree->SetBranchStatus("Jet_nConstituents", 0);
  outputtree->SetBranchStatus("Jet_nElectrons", 0);
  outputtree->SetBranchStatus("Jet_nMuons", 0);
  outputtree->SetBranchStatus("Jet_puId", 0);
  outputtree->SetBranchStatus("LHE_HT", 0);
  outputtree->SetBranchStatus("LHE_HTIncoming", 0);
  outputtree->SetBranchStatus("LHE_Vpt", 0);
  outputtree->SetBranchStatus("LHE_Njets", 0);
  outputtree->SetBranchStatus("LHE_Nb", 0);
  outputtree->SetBranchStatus("LHE_Nc", 0);
  outputtree->SetBranchStatus("LHE_Nuds", 0);
  outputtree->SetBranchStatus("LHE_Nglu", 0);
  outputtree->SetBranchStatus("LHE_NpNLO", 0);
  outputtree->SetBranchStatus("LHE_NpLO", 0);
  outputtree->SetBranchStatus("nLHEPart", 0);
  outputtree->SetBranchStatus("LHEPart_pt", 0);
  outputtree->SetBranchStatus("LHEPart_eta", 0);
  outputtree->SetBranchStatus("LHEPart_phi", 0);
  outputtree->SetBranchStatus("LHEPart_mass", 0);
  outputtree->SetBranchStatus("LHEPart_pdgId", 0);
  outputtree->SetBranchStatus("GenMET_phi", 0);
  outputtree->SetBranchStatus("GenMET_pt", 0);
  outputtree->SetBranchStatus("MET_MetUnclustEnUpDeltaX", 0);
  outputtree->SetBranchStatus("MET_MetUnclustEnUpDeltaY", 0);
  outputtree->SetBranchStatus("MET_covXX", 0);
  outputtree->SetBranchStatus("MET_covXY", 0);
  outputtree->SetBranchStatus("MET_covYY", 0);
  outputtree->SetBranchStatus("MET_phi", 0);
  outputtree->SetBranchStatus("MET_pt", 0);
  outputtree->SetBranchStatus("MET_significance", 0);
  outputtree->SetBranchStatus("MET_sumEt", 0);
  outputtree->SetBranchStatus("nMuon", 0);
  outputtree->SetBranchStatus("Muon_dxy", 0);
  outputtree->SetBranchStatus("Muon_dxyErr", 0);
  outputtree->SetBranchStatus("Muon_dz", 0);
  outputtree->SetBranchStatus("Muon_dzErr", 0);
  outputtree->SetBranchStatus("Muon_eta", 0);
  outputtree->SetBranchStatus("Muon_ip3d", 0);
  outputtree->SetBranchStatus("Muon_mass", 0);
  outputtree->SetBranchStatus("Muon_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("Muon_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso04_all", 0);
  outputtree->SetBranchStatus("Muon_phi", 0);
  outputtree->SetBranchStatus("Muon_pt", 0);
  outputtree->SetBranchStatus("Muon_ptErr", 0);
  outputtree->SetBranchStatus("Muon_segmentComp", 0);
  outputtree->SetBranchStatus("Muon_sip3d", 0);
  outputtree->SetBranchStatus("Muon_mvaTTH", 0);
  outputtree->SetBranchStatus("Muon_charge", 0);
  outputtree->SetBranchStatus("Muon_jetIdx", 0);
  outputtree->SetBranchStatus("Muon_nStations", 0);
  outputtree->SetBranchStatus("Muon_nTrackerLayers", 0);
  outputtree->SetBranchStatus("Muon_pdgId", 0);
  outputtree->SetBranchStatus("Muon_tightCharge", 0);
  outputtree->SetBranchStatus("Muon_highPtId", 0);
  outputtree->SetBranchStatus("Muon_isPFcand", 0);
  outputtree->SetBranchStatus("Muon_mediumId", 0);
  outputtree->SetBranchStatus("Muon_softId", 0);
  outputtree->SetBranchStatus("Muon_tightId", 0);
  outputtree->SetBranchStatus("nPhoton", 0);
  outputtree->SetBranchStatus("Photon_eCorr", 0);
  outputtree->SetBranchStatus("Photon_energyErr", 0);
  outputtree->SetBranchStatus("Photon_eta", 0);
  outputtree->SetBranchStatus("Photon_hoe", 0);
  outputtree->SetBranchStatus("Photon_mass", 0);
  outputtree->SetBranchStatus("Photon_mvaID", 0);
  outputtree->SetBranchStatus("Photon_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("Photon_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("Photon_phi", 0);
  outputtree->SetBranchStatus("Photon_pt", 0);
  outputtree->SetBranchStatus("Photon_r9", 0);
  outputtree->SetBranchStatus("Photon_sieie", 0);
  outputtree->SetBranchStatus("Photon_charge", 0);
  outputtree->SetBranchStatus("Photon_cutBasedBitmap", 0);
  outputtree->SetBranchStatus("Photon_electronIdx", 0);
  outputtree->SetBranchStatus("Photon_jetIdx", 0);
  outputtree->SetBranchStatus("Photon_pdgId", 0);
  outputtree->SetBranchStatus("Photon_vidNestedWPBitmap", 0);
  outputtree->SetBranchStatus("Photon_electronVeto", 0);
  outputtree->SetBranchStatus("Photon_isScEtaEB", 0);
  outputtree->SetBranchStatus("Photon_isScEtaEE", 0);
  outputtree->SetBranchStatus("Photon_mvaID_WP80", 0);
  outputtree->SetBranchStatus("Photon_mvaID_WP90", 0);
  outputtree->SetBranchStatus("Photon_pixelSeed", 0);
  outputtree->SetBranchStatus("Pileup_nTrueInt", 0);
  outputtree->SetBranchStatus("Pileup_nPU", 0);
  outputtree->SetBranchStatus("Pileup_sumEOOT", 0);
  outputtree->SetBranchStatus("Pileup_sumLOOT", 0);
  outputtree->SetBranchStatus("PuppiMET_phi", 0);
  outputtree->SetBranchStatus("PuppiMET_pt", 0);
  outputtree->SetBranchStatus("PuppiMET_sumEt", 0);
  outputtree->SetBranchStatus("RawMET_phi", 0);
  outputtree->SetBranchStatus("RawMET_pt", 0);
  outputtree->SetBranchStatus("RawMET_sumEt", 0);
  outputtree->SetBranchStatus("fixedGridRhoFastjetAll", 0);
  outputtree->SetBranchStatus("fixedGridRhoFastjetCentralCalo", 0);
  outputtree->SetBranchStatus("fixedGridRhoFastjetCentralNeutral", 0);
  outputtree->SetBranchStatus("nGenDressedLepton", 0);
  outputtree->SetBranchStatus("GenDressedLepton_eta", 0);
  outputtree->SetBranchStatus("GenDressedLepton_mass", 0);
  outputtree->SetBranchStatus("GenDressedLepton_phi", 0);
  outputtree->SetBranchStatus("GenDressedLepton_pt", 0);
  outputtree->SetBranchStatus("GenDressedLepton_pdgId", 0);
  outputtree->SetBranchStatus("nSoftActivityJet", 0);
  outputtree->SetBranchStatus("SoftActivityJet_eta", 0);
  outputtree->SetBranchStatus("SoftActivityJet_phi", 0);
  outputtree->SetBranchStatus("SoftActivityJet_pt", 0);
  outputtree->SetBranchStatus("SoftActivityJetHT", 0);
  outputtree->SetBranchStatus("SoftActivityJetHT10", 0);
  outputtree->SetBranchStatus("SoftActivityJetHT2", 0);
  outputtree->SetBranchStatus("SoftActivityJetHT5", 0);
  outputtree->SetBranchStatus("SoftActivityJetNjets10", 0);
  outputtree->SetBranchStatus("SoftActivityJetNjets2", 0);
  outputtree->SetBranchStatus("SoftActivityJetNjets5", 0);
  outputtree->SetBranchStatus("nSubJet", 0);
  outputtree->SetBranchStatus("SubJet_btagCMVA", 0);
  outputtree->SetBranchStatus("SubJet_btagCSVV2", 0);
  outputtree->SetBranchStatus("SubJet_btagDeepB", 0);
  outputtree->SetBranchStatus("SubJet_eta", 0);
  outputtree->SetBranchStatus("SubJet_mass", 0);
  outputtree->SetBranchStatus("SubJet_n2b1", 0);
  outputtree->SetBranchStatus("SubJet_n3b1", 0);
  outputtree->SetBranchStatus("SubJet_phi", 0);
  outputtree->SetBranchStatus("SubJet_pt", 0);
  outputtree->SetBranchStatus("SubJet_tau1", 0);
  outputtree->SetBranchStatus("SubJet_tau2", 0);
  outputtree->SetBranchStatus("SubJet_tau3", 0);
  outputtree->SetBranchStatus("SubJet_tau4", 0);
  outputtree->SetBranchStatus("nTau", 0);
  outputtree->SetBranchStatus("Tau_chargedIso", 0);
  outputtree->SetBranchStatus("Tau_dxy", 0);
  outputtree->SetBranchStatus("Tau_dz", 0);
  outputtree->SetBranchStatus("Tau_eta", 0);
  outputtree->SetBranchStatus("Tau_leadTkDeltaEta", 0);
  outputtree->SetBranchStatus("Tau_leadTkDeltaPhi", 0);
  outputtree->SetBranchStatus("Tau_leadTkPtOverTauPt", 0);
  outputtree->SetBranchStatus("Tau_mass", 0);
  outputtree->SetBranchStatus("Tau_neutralIso", 0);
  outputtree->SetBranchStatus("Tau_phi", 0);
  outputtree->SetBranchStatus("Tau_photonsOutsideSignalCone", 0);
  outputtree->SetBranchStatus("Tau_pt", 0);
  outputtree->SetBranchStatus("Tau_puCorr", 0);
  outputtree->SetBranchStatus("Tau_rawAntiEle", 0);
  outputtree->SetBranchStatus("Tau_rawIso", 0);
  outputtree->SetBranchStatus("Tau_rawIsodR03", 0);
  outputtree->SetBranchStatus("Tau_rawMVAnewDM2017v2", 0);
  outputtree->SetBranchStatus("Tau_rawMVAoldDM", 0);
  outputtree->SetBranchStatus("Tau_rawMVAoldDM2017v1", 0);
  outputtree->SetBranchStatus("Tau_rawMVAoldDM2017v2", 0);
  outputtree->SetBranchStatus("Tau_rawMVAoldDMdR032017v2", 0);
  outputtree->SetBranchStatus("Tau_charge", 0);
  outputtree->SetBranchStatus("Tau_decayMode", 0);
  outputtree->SetBranchStatus("Tau_jetIdx", 0);
  outputtree->SetBranchStatus("Tau_rawAntiEleCat", 0);
  outputtree->SetBranchStatus("Tau_idAntiEle", 0);
  outputtree->SetBranchStatus("Tau_idAntiMu", 0);
  outputtree->SetBranchStatus("Tau_idDecayMode", 0);
  outputtree->SetBranchStatus("Tau_idDecayModeNewDMs", 0);
  outputtree->SetBranchStatus("Tau_idMVAnewDM2017v2", 0);
  outputtree->SetBranchStatus("Tau_idMVAoldDM", 0);
  outputtree->SetBranchStatus("Tau_idMVAoldDM2017v1", 0);
  outputtree->SetBranchStatus("Tau_idMVAoldDM2017v2", 0);
  outputtree->SetBranchStatus("Tau_idMVAoldDMdR032017v2", 0);
  outputtree->SetBranchStatus("TkMET_phi", 0);
  outputtree->SetBranchStatus("TkMET_pt", 0);
  outputtree->SetBranchStatus("TkMET_sumEt", 0);
  outputtree->SetBranchStatus("nTrigObj", 0);
  outputtree->SetBranchStatus("TrigObj_pt", 0);
  outputtree->SetBranchStatus("TrigObj_eta", 0);
  outputtree->SetBranchStatus("TrigObj_phi", 0);
  outputtree->SetBranchStatus("TrigObj_l1pt", 0);
  outputtree->SetBranchStatus("TrigObj_l1pt_2", 0);
  outputtree->SetBranchStatus("TrigObj_l2pt", 0);
  outputtree->SetBranchStatus("TrigObj_id", 0);
  outputtree->SetBranchStatus("TrigObj_l1iso", 0);
  outputtree->SetBranchStatus("TrigObj_l1charge", 0);
  outputtree->SetBranchStatus("TrigObj_filterBits", 0);
  outputtree->SetBranchStatus("genTtbarId", 0);
  outputtree->SetBranchStatus("nOtherPV", 0);
  outputtree->SetBranchStatus("OtherPV_z", 0);
  outputtree->SetBranchStatus("PV_ndof", 0);
  outputtree->SetBranchStatus("PV_x", 0);
  outputtree->SetBranchStatus("PV_y", 0);
  outputtree->SetBranchStatus("PV_z", 0);
  outputtree->SetBranchStatus("PV_chi2", 0);
  outputtree->SetBranchStatus("PV_score", 0);
  outputtree->SetBranchStatus("PV_npvs", 0);
  outputtree->SetBranchStatus("PV_npvsGood", 0);
  outputtree->SetBranchStatus("nSV", 0);
  outputtree->SetBranchStatus("SV_dlen", 0);
  outputtree->SetBranchStatus("SV_dlenSig", 0);
  outputtree->SetBranchStatus("SV_pAngle", 0);
  outputtree->SetBranchStatus("Electron_genPartIdx", 0);
  outputtree->SetBranchStatus("Electron_genPartFlav", 0);
  outputtree->SetBranchStatus("GenJetAK8_partonFlavour", 0);
  outputtree->SetBranchStatus("GenJetAK8_hadronFlavour", 0);
  outputtree->SetBranchStatus("GenJet_partonFlavour", 0);
  outputtree->SetBranchStatus("GenJet_hadronFlavour", 0);
  outputtree->SetBranchStatus("Jet_genJetIdx", 0);
  outputtree->SetBranchStatus("Jet_hadronFlavour", 0);
  outputtree->SetBranchStatus("Jet_partonFlavour", 0);
  outputtree->SetBranchStatus("Muon_genPartIdx", 0);
  outputtree->SetBranchStatus("Muon_genPartFlav", 0);
  outputtree->SetBranchStatus("Photon_genPartIdx", 0);
  outputtree->SetBranchStatus("Photon_genPartFlav", 0);
  outputtree->SetBranchStatus("MET_fiducialGenPhi", 0);
  outputtree->SetBranchStatus("MET_fiducialGenPt", 0);
  outputtree->SetBranchStatus("Electron_cleanmask", 0);
  outputtree->SetBranchStatus("Jet_cleanmask", 0);
  outputtree->SetBranchStatus("Muon_cleanmask", 0);
  outputtree->SetBranchStatus("Photon_cleanmask", 0);
  outputtree->SetBranchStatus("Tau_cleanmask", 0);
  outputtree->SetBranchStatus("SV_chi2", 0);
  outputtree->SetBranchStatus("SV_eta", 0);
  outputtree->SetBranchStatus("SV_mass", 0);
  outputtree->SetBranchStatus("SV_ndof", 0);
  outputtree->SetBranchStatus("SV_phi", 0);
  outputtree->SetBranchStatus("SV_pt", 0);
  outputtree->SetBranchStatus("SV_x", 0);
  outputtree->SetBranchStatus("SV_y", 0);
  outputtree->SetBranchStatus("SV_z", 0);
  outputtree->SetBranchStatus("Tau_genPartIdx", 0);
  outputtree->SetBranchStatus("Tau_genPartFlav", 0);
  outputtree->SetBranchStatus("L1simulation_step", 0);
  outputtree->SetBranchStatus("HLTriggerFirstPath", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet360_TrimMass30", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet380_TrimMass30", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet400_TrimMass30", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet420_TrimMass30", 0);
  outputtree->SetBranchStatus("HLT_AK8PFHT750_TrimMass50", 0);
  outputtree->SetBranchStatus("HLT_AK8PFHT800_TrimMass50", 0);
  outputtree->SetBranchStatus("HLT_AK8PFHT850_TrimMass50", 0);
  outputtree->SetBranchStatus("HLT_AK8PFHT900_TrimMass50", 0);
  outputtree->SetBranchStatus("HLT_CaloJet500_NoJetID", 0);
  outputtree->SetBranchStatus("HLT_CaloJet550_NoJetID", 0);
  outputtree->SetBranchStatus("HLT_Trimuon5_3p5_2_Upsilon_Muon", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle25_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle27_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle33_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle24_eta2p1_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_Ele27_Ele37_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_Mu27_Ele37_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_Mu37_Ele27_CaloIdL_MW", 0);
  outputtree->SetBranchStatus("HLT_Mu37_TkMu27", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_3_Bs", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_3_Jpsi_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_JpsiTrk_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_Trk_Tau3mu", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_PsiPrimeTrk_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_Mass8_DZ_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu8_Mass8_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_Mu3_PFJet40", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_L2Mu2_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_L2Mu2_Upsilon", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track2_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track3p5_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track7_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track2_Upsilon", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track3p5_Upsilon", 0);
  outputtree->SetBranchStatus("HLT_Mu7p5_Track7_Upsilon", 0);
  outputtree->SetBranchStatus("HLT_DoublePhoton33_CaloIdL", 0);
  outputtree->SetBranchStatus("HLT_DoublePhoton70", 0);
  outputtree->SetBranchStatus("HLT_DoublePhoton85", 0);
  outputtree->SetBranchStatus("HLT_Ele20_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele20_WPLoose_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele20_eta2p1_WPLoose_Gsf", 0);
  outputtree->SetBranchStatus("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", 0);
  outputtree->SetBranchStatus("HLT_Ele27_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele35_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele35_WPTight_Gsf_L1EGMT", 0);
  outputtree->SetBranchStatus("HLT_Ele38_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele40_WPTight_Gsf", 0);
  outputtree->SetBranchStatus("HLT_Ele32_WPTight_Gsf_L1DoubleEG", 0);
  outputtree->SetBranchStatus("HLT_HT450_Beamspot", 0);
  outputtree->SetBranchStatus("HLT_HT300_Beamspot", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu20", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu27", 0);
  outputtree->SetBranchStatus("HLT_IsoMu30", 0);
  outputtree->SetBranchStatus("HLT_UncorrectedJetE30_NoBPTX", 0);
  outputtree->SetBranchStatus("HLT_UncorrectedJetE30_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_UncorrectedJetE60_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_UncorrectedJetE70_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_L1SingleMu18", 0);
  outputtree->SetBranchStatus("HLT_L1SingleMu25", 0);
  outputtree->SetBranchStatus("HLT_L2Mu10", 0);
  outputtree->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX", 0);
  outputtree->SetBranchStatus("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", 0);
  outputtree->SetBranchStatus("HLT_L2Mu50", 0);
  outputtree->SetBranchStatus("HLT_DoubleL2Mu50", 0);
  outputtree->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", 0);
  outputtree->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", 0);
  outputtree->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 0);
  outputtree->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", 0);
  outputtree->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", 0);
  outputtree->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", 0);
  outputtree->SetBranchStatus("HLT_Mu25_TkMu0_Onia", 0);
  outputtree->SetBranchStatus("HLT_Mu30_TkMu0_Onia", 0);
  outputtree->SetBranchStatus("HLT_Mu20_TkMu0_Phi", 0);
  outputtree->SetBranchStatus("HLT_Mu25_TkMu0_Phi", 0);
  outputtree->SetBranchStatus("HLT_Mu20", 0);
  outputtree->SetBranchStatus("HLT_Mu27", 0);
  outputtree->SetBranchStatus("HLT_Mu50", 0);
  outputtree->SetBranchStatus("HLT_Mu55", 0);
  outputtree->SetBranchStatus("HLT_OldMu100", 0);
  outputtree->SetBranchStatus("HLT_TkMu100", 0);
  outputtree->SetBranchStatus("HLT_DiPFJet15_NoCaloMatched", 0);
  outputtree->SetBranchStatus("HLT_DiPFJet25_NoCaloMatched", 0);
  outputtree->SetBranchStatus("HLT_DiPFJet15_FBEta3_NoCaloMatched", 0);
  outputtree->SetBranchStatus("HLT_DiPFJet25_FBEta3_NoCaloMatched", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve40", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve60", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve80", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve140", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve200", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve260", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve320", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve400", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve500", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve15_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve25_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve35_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve60_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve80_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve100_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve160_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve220_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_DiPFJetAve300_HFJEC", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet40", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet60", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet80", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet140", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet200", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet260", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet320", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet400", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet450", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet500", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet550", 0);
  outputtree->SetBranchStatus("HLT_PFJet40", 0);
  outputtree->SetBranchStatus("HLT_PFJet60", 0);
  outputtree->SetBranchStatus("HLT_PFJet80", 0);
  outputtree->SetBranchStatus("HLT_PFJet140", 0);
  outputtree->SetBranchStatus("HLT_PFJet200", 0);
  outputtree->SetBranchStatus("HLT_PFJet260", 0);
  outputtree->SetBranchStatus("HLT_PFJet320", 0);
  outputtree->SetBranchStatus("HLT_PFJet400", 0);
  outputtree->SetBranchStatus("HLT_PFJet450", 0);
  outputtree->SetBranchStatus("HLT_PFJet500", 0);
  outputtree->SetBranchStatus("HLT_PFJet550", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd40", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd60", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd80", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd140", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd200", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd260", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd320", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd400", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd450", 0);
  outputtree->SetBranchStatus("HLT_PFJetFwd500", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd40", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd60", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd80", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd140", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd200", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd260", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd320", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd400", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd450", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJetFwd500", 0);
  outputtree->SetBranchStatus("HLT_PFHT180", 0);
  outputtree->SetBranchStatus("HLT_PFHT250", 0);
  outputtree->SetBranchStatus("HLT_PFHT370", 0);
  outputtree->SetBranchStatus("HLT_PFHT430", 0);
  outputtree->SetBranchStatus("HLT_PFHT510", 0);
  outputtree->SetBranchStatus("HLT_PFHT590", 0);
  outputtree->SetBranchStatus("HLT_PFHT680", 0);
  outputtree->SetBranchStatus("HLT_PFHT780", 0);
  outputtree->SetBranchStatus("HLT_PFHT890", 0);
  outputtree->SetBranchStatus("HLT_PFHT1050", 0);
  outputtree->SetBranchStatus("HLT_PFHT500_PFMET100_PFMHT100_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFHT500_PFMET110_PFMHT110_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFHT700_PFMET85_PFMHT85_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFHT700_PFMET95_PFMHT95_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFHT800_PFMET75_PFMHT75_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFHT800_PFMET85_PFMHT85_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", 0);
  outputtree->SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", 0);
  outputtree->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1", 0);
  outputtree->SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1", 0);
  outputtree->SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1", 0);
  outputtree->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne110_PFMHT110_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne130_PFMHT130_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne140_PFMHT140_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", 0);
  outputtree->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", 0);
  outputtree->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", 0);
  outputtree->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", 0);
  outputtree->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", 0);
  outputtree->SetBranchStatus("HLT_L1ETMHadSeeds", 0);
  outputtree->SetBranchStatus("HLT_CaloMHT90", 0);
  outputtree->SetBranchStatus("HLT_CaloMET80_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET90_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET100_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET110_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET250_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET70_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET80_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET90_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET100_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET250_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET300_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_CaloMET350_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMET200_NotCleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMET200_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMET250_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMET300_HBHECleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMET200_HBHE_BeamHaloCleaned", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", 0);
  outputtree->SetBranchStatus("HLT_MET105_IsoTrk50", 0);
  outputtree->SetBranchStatus("HLT_MET120_IsoTrk50", 0);
  outputtree->SetBranchStatus("HLT_SingleJet30_Mu12_SinglePFJet40", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets40_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets100_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets200_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets350_CaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33", 0);
  outputtree->SetBranchStatus("HLT_Photon300_NoHE", 0);
  outputtree->SetBranchStatus("HLT_Mu8_TrkIsoVVL", 0);
  outputtree->SetBranchStatus("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", 0);
  outputtree->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 0);
  outputtree->SetBranchStatus("HLT_Mu17_TrkIsoVVL", 0);
  outputtree->SetBranchStatus("HLT_Mu19_TrkIsoVVL", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4DiJet20_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4DiJet40_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4DiJet70_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4DiJet110_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4DiJet170_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK4Jet300_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK8DiJet170_Mu5", 0);
  outputtree->SetBranchStatus("HLT_BTagMu_AK8Jet300_Mu5", 0);
  outputtree->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", 0);
  outputtree->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", 0);
  outputtree->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 0);
  outputtree->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu12_DoublePhoton20", 0);
  outputtree->SetBranchStatus("HLT_TriplePhoton_20_20_20_CaloIdLV2", 0);
  outputtree->SetBranchStatus("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", 0);
  outputtree->SetBranchStatus("HLT_TriplePhoton_30_30_10_CaloIdLV2", 0);
  outputtree->SetBranchStatus("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", 0);
  outputtree->SetBranchStatus("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", 0);
  outputtree->SetBranchStatus("HLT_Photon25", 0);
  outputtree->SetBranchStatus("HLT_Photon33", 0);
  outputtree->SetBranchStatus("HLT_Photon50", 0);
  outputtree->SetBranchStatus("HLT_Photon75", 0);
  outputtree->SetBranchStatus("HLT_Photon90", 0);
  outputtree->SetBranchStatus("HLT_Photon120", 0);
  outputtree->SetBranchStatus("HLT_Photon150", 0);
  outputtree->SetBranchStatus("HLT_Photon175", 0);
  outputtree->SetBranchStatus("HLT_Photon200", 0);
  outputtree->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM", 0);
  outputtree->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM", 0);
  outputtree->SetBranchStatus("HLT_Photon90_R9Id90_HE10_IsoM", 0);
  outputtree->SetBranchStatus("HLT_Photon120_R9Id90_HE10_IsoM", 0);
  outputtree->SetBranchStatus("HLT_Photon165_R9Id90_HE10_IsoM", 0);
  outputtree->SetBranchStatus("HLT_Photon90_CaloIdL_PFHT700", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi_L1_NoOS", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Jpsi3p5_Muon2", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_5", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5NoOS", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5er2p0", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_NoVertexing", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_5M", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass_L1_0er1p5R", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass_L1_0er1p5", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass_L1_4", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass_L1_4R", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_LowMass_L1_TM530", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_Muon_L1_TM0", 0);
  outputtree->SetBranchStatus("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", 0);
  outputtree->SetBranchStatus("HLT_TripleMu_5_3_3_Mass3p8to60_DZ", 0);
  outputtree->SetBranchStatus("HLT_TripleMu_10_5_5_DZ", 0);
  outputtree->SetBranchStatus("HLT_TripleMu_12_10_5", 0);
  outputtree->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", 0);
  outputtree->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", 0);
  outputtree->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", 0);
  outputtree->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_Jpsi_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_Jpsi_NoVertexing", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu4_JpsiTrkTrk_Displaced", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu43NoFiltersNoVtx", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu48NoFiltersNoVtx", 0);
  outputtree->SetBranchStatus("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", 0);
  outputtree->SetBranchStatus("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_L1_DM4", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", 0);
  outputtree->SetBranchStatus("HLT_HT425", 0);
  outputtree->SetBranchStatus("HLT_HT430_DisplacedDijet40_DisplacedTrack", 0);
  outputtree->SetBranchStatus("HLT_HT430_DisplacedDijet60_DisplacedTrack", 0);
  outputtree->SetBranchStatus("HLT_HT430_DisplacedDijet80_DisplacedTrack", 0);
  outputtree->SetBranchStatus("HLT_HT400_DisplacedDijet40_DisplacedTrack", 0);
  outputtree->SetBranchStatus("HLT_HT650_DisplacedDijet60_Inclusive", 0);
  outputtree->SetBranchStatus("HLT_HT550_DisplacedDijet80_Inclusive", 0);
  outputtree->SetBranchStatus("HLT_HT550_DisplacedDijet60_Inclusive", 0);
  outputtree->SetBranchStatus("HLT_HT650_DisplacedDijet80_Inclusive", 0);
  outputtree->SetBranchStatus("HLT_HT750_DisplacedDijet80_Inclusive", 0);
  outputtree->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET110", 0);
  outputtree->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET120", 0);
  outputtree->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET130", 0);
  outputtree->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET110", 0);
  outputtree->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET120", 0);
  outputtree->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET130", 0);
  outputtree->SetBranchStatus("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", 0);
  outputtree->SetBranchStatus("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", 0);
  outputtree->SetBranchStatus("HLT_Ele28_HighEta_SC20_Mass55", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_Photon23", 0);
  outputtree->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", 0);
  outputtree->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", 0);
  outputtree->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450", 0);
  outputtree->SetBranchStatus("HLT_Ele50_IsoVVVL_PFHT450", 0);
  outputtree->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT600", 0);
  outputtree->SetBranchStatus("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", 0);
  outputtree->SetBranchStatus("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", 0);
  outputtree->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", 0);
  outputtree->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", 0);
  outputtree->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450", 0);
  outputtree->SetBranchStatus("HLT_Mu50_IsoVVVL_PFHT450", 0);
  outputtree->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT600", 0);
  outputtree->SetBranchStatus("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", 0);
  outputtree->SetBranchStatus("HLT_Dimuon20_Jpsi_Barrel_Seagulls", 0);
  outputtree->SetBranchStatus("HLT_Dimuon10_Upsilon_Barrel_Seagulls", 0);
  outputtree->SetBranchStatus("HLT_Dimuon12_Upsilon_eta1p5", 0);
  outputtree->SetBranchStatus("HLT_Dimuon14_Phi_Barrel_Seagulls", 0);
  outputtree->SetBranchStatus("HLT_Dimuon18_PsiPrime", 0);
  outputtree->SetBranchStatus("HLT_Dimuon25_Jpsi", 0);
  outputtree->SetBranchStatus("HLT_Dimuon18_PsiPrime_noCorrL1", 0);
  outputtree->SetBranchStatus("HLT_Dimuon24_Upsilon_noCorrL1", 0);
  outputtree->SetBranchStatus("HLT_Dimuon24_Phi_noCorrL1", 0);
  outputtree->SetBranchStatus("HLT_Dimuon25_Jpsi_noCorrL1", 0);
  outputtree->SetBranchStatus("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", 0);
  outputtree->SetBranchStatus("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", 0);
  outputtree->SetBranchStatus("HLT_DoubleIsoMu20_eta2p1", 0);
  outputtree->SetBranchStatus("HLT_DoubleIsoMu24_eta2p1", 0);
  outputtree->SetBranchStatus("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", 0);
  outputtree->SetBranchStatus("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", 0);
  outputtree->SetBranchStatus("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", 0);
  outputtree->SetBranchStatus("HLT_Mu8", 0);
  outputtree->SetBranchStatus("HLT_Mu17", 0);
  outputtree->SetBranchStatus("HLT_Mu19", 0);
  outputtree->SetBranchStatus("HLT_Mu17_Photon30_IsoCaloId", 0);
  outputtree->SetBranchStatus("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", 0);
  outputtree->SetBranchStatus("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", 0);
  outputtree->SetBranchStatus("HLT_Ele115_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_Ele135_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_Ele145_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_Ele200_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_Ele250_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_Ele300_CaloIdVT_GsfTrkIdT", 0);
  outputtree->SetBranchStatus("HLT_PFHT300PT30_QuadPFJet_75_60_45_40", 0);
  outputtree->SetBranchStatus("HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0", 0);
  outputtree->SetBranchStatus("HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2", 0);
  outputtree->SetBranchStatus("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2", 0);
  outputtree->SetBranchStatus("HLT_PFHT380_SixPFJet32", 0);
  outputtree->SetBranchStatus("HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5", 0);
  outputtree->SetBranchStatus("HLT_PFHT430_SixPFJet40", 0);
  outputtree->SetBranchStatus("HLT_PFHT350", 0);
  outputtree->SetBranchStatus("HLT_PFHT350MinPFJet15", 0);
  outputtree->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL", 0);
  outputtree->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", 0);
  outputtree->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", 0);
  outputtree->SetBranchStatus("HLT_FullTrack_Multiplicity85", 0);
  outputtree->SetBranchStatus("HLT_FullTrack_Multiplicity100", 0);
  outputtree->SetBranchStatus("HLT_FullTrack_Multiplicity130", 0);
  outputtree->SetBranchStatus("HLT_FullTrack_Multiplicity155", 0);
  outputtree->SetBranchStatus("HLT_ECALHT800", 0);
  outputtree->SetBranchStatus("HLT_DiSC30_18_EIso_AND_HE_Mass70", 0);
  outputtree->SetBranchStatus("HLT_Physics", 0);
  outputtree->SetBranchStatus("HLT_Physics_part0", 0);
  outputtree->SetBranchStatus("HLT_Physics_part1", 0);
  outputtree->SetBranchStatus("HLT_Physics_part2", 0);
  outputtree->SetBranchStatus("HLT_Physics_part3", 0);
  outputtree->SetBranchStatus("HLT_Physics_part4", 0);
  outputtree->SetBranchStatus("HLT_Physics_part5", 0);
  outputtree->SetBranchStatus("HLT_Physics_part6", 0);
  outputtree->SetBranchStatus("HLT_Physics_part7", 0);
  outputtree->SetBranchStatus("HLT_Random", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part0", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part1", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part2", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part3", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part4", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part5", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part6", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_part7", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet30", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet40", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet50", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet80", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet100", 0);
  outputtree->SetBranchStatus("HLT_AK4CaloJet120", 0);
  outputtree->SetBranchStatus("HLT_AK4PFJet30", 0);
  outputtree->SetBranchStatus("HLT_AK4PFJet50", 0);
  outputtree->SetBranchStatus("HLT_AK4PFJet80", 0);
  outputtree->SetBranchStatus("HLT_AK4PFJet100", 0);
  outputtree->SetBranchStatus("HLT_AK4PFJet120", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton10_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton20_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton30_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton40_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton50_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_HISinglePhoton60_Eta3p1ForPPRef", 0);
  outputtree->SetBranchStatus("HLT_Photon20_HoverELoose", 0);
  outputtree->SetBranchStatus("HLT_Photon30_HoverELoose", 0);
  outputtree->SetBranchStatus("HLT_Photon40_HoverELoose", 0);
  outputtree->SetBranchStatus("HLT_Photon50_HoverELoose", 0);
  outputtree->SetBranchStatus("HLT_Photon60_HoverELoose", 0);
  outputtree->SetBranchStatus("HLT_EcalCalibration", 0);
  outputtree->SetBranchStatus("HLT_HcalCalibration", 0);
  outputtree->SetBranchStatus("HLT_L1UnpairedBunchBptxMinus", 0);
  outputtree->SetBranchStatus("HLT_L1UnpairedBunchBptxPlus", 0);
  outputtree->SetBranchStatus("HLT_L1NotBptxOR", 0);
  outputtree->SetBranchStatus("HLT_L1MinimumBiasHF_OR", 0);
  outputtree->SetBranchStatus("HLT_L1MinimumBiasHF0OR", 0);
  outputtree->SetBranchStatus("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", 0);
  outputtree->SetBranchStatus("HLT_HcalNZS", 0);
  outputtree->SetBranchStatus("HLT_HcalPhiSym", 0);
  outputtree->SetBranchStatus("HLT_HcalIsolatedbunch", 0);
  outputtree->SetBranchStatus("HLT_IsoTrackHB", 0);
  outputtree->SetBranchStatus("HLT_IsoTrackHE", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_FirstCollisionAfterAbortGap", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_IsolatedBunches", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_FirstCollisionInTrain", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_LastCollisionInTrain", 0);
  outputtree->SetBranchStatus("HLT_ZeroBias_FirstBXAfterTrain", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", 0);
  outputtree->SetBranchStatus("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", 0);
  outputtree->SetBranchStatus("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", 0);
  outputtree->SetBranchStatus("HLT_Rsq0p35", 0);
  outputtree->SetBranchStatus("HLT_Rsq0p40", 0);
  outputtree->SetBranchStatus("HLT_RsqMR300_Rsq0p09_MR200", 0);
  outputtree->SetBranchStatus("HLT_RsqMR320_Rsq0p09_MR200", 0);
  outputtree->SetBranchStatus("HLT_RsqMR300_Rsq0p09_MR200_4jet", 0);
  outputtree->SetBranchStatus("HLT_RsqMR320_Rsq0p09_MR200_4jet", 0);
  outputtree->SetBranchStatus("HLT_L1_DoubleJet30_Mass_Min400_Mu10", 0);
  outputtree->SetBranchStatus("HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1", 0);
  outputtree->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", 0);
  outputtree->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", 0);
  outputtree->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", 0);
  outputtree->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", 0);
  outputtree->SetBranchStatus("HLT_Mu18_Mu9_SameSign", 0);
  outputtree->SetBranchStatus("HLT_Mu18_Mu9_SameSign_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu18_Mu9", 0);
  outputtree->SetBranchStatus("HLT_Mu18_Mu9_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu20_Mu10_SameSign", 0);
  outputtree->SetBranchStatus("HLT_Mu20_Mu10_SameSign_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu20_Mu10", 0);
  outputtree->SetBranchStatus("HLT_Mu20_Mu10_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu23_Mu12_SameSign", 0);
  outputtree->SetBranchStatus("HLT_Mu23_Mu12_SameSign_DZ", 0);
  outputtree->SetBranchStatus("HLT_Mu23_Mu12", 0);
  outputtree->SetBranchStatus("HLT_Mu23_Mu12_DZ", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", 0);
  outputtree->SetBranchStatus("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", 0);
  outputtree->SetBranchStatus("HLT_TripleMu_5_3_3_Mass3p8to60_DCA", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet98_83_71_15", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet103_88_75_15", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet105_88_76_15", 0);
  outputtree->SetBranchStatus("HLT_QuadPFJet111_90_80_15", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet330_PFAK8BTagCSV_p17", 0);
  outputtree->SetBranchStatus("HLT_AK8PFJet330_PFAK8BTagCSV_p1", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", 0);
  outputtree->SetBranchStatus("HLTriggerFinalPath", 0);
  outputtree->SetBranchStatus("Flag_HBHENoiseFilter", 0);
  outputtree->SetBranchStatus("Flag_HBHENoiseIsoFilter", 0);
  outputtree->SetBranchStatus("Flag_CSCTightHaloFilter", 0);
  outputtree->SetBranchStatus("Flag_CSCTightHaloTrkMuUnvetoFilter", 0);
  outputtree->SetBranchStatus("Flag_CSCTightHalo2015Filter", 0);
  outputtree->SetBranchStatus("Flag_globalTightHalo2016Filter", 0);
  outputtree->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 0);
  outputtree->SetBranchStatus("Flag_HcalStripHaloFilter", 0);
  outputtree->SetBranchStatus("Flag_hcalLaserEventFilter", 0);
  outputtree->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 0);
  outputtree->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 0);
  outputtree->SetBranchStatus("Flag_ecalBadCalibFilter", 0);
  outputtree->SetBranchStatus("Flag_goodVertices", 0);
  outputtree->SetBranchStatus("Flag_eeBadScFilter", 0);
  outputtree->SetBranchStatus("Flag_ecalLaserCorrFilter", 0);
  outputtree->SetBranchStatus("Flag_trkPOGFilters", 0);
  outputtree->SetBranchStatus("Flag_chargedHadronTrackResolutionFilter", 0);
  outputtree->SetBranchStatus("Flag_muonBadTrackFilter", 0);
  outputtree->SetBranchStatus("Flag_BadChargedCandidateFilter", 0);
  outputtree->SetBranchStatus("Flag_BadPFMuonFilter", 0);
  outputtree->SetBranchStatus("Flag_BadChargedCandidateSummer16Filter", 0);
  outputtree->SetBranchStatus("Flag_BadPFMuonSummer16Filter", 0);
  outputtree->SetBranchStatus("Flag_trkPOG_manystripclus53X", 0);
  outputtree->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 0);
  outputtree->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 0);
  outputtree->SetBranchStatus("Flag_METFilters", 0);
  outputtree->SetBranchStatus("nLepton", 0);
  outputtree->SetBranchStatus("Lepton_pdgId", 0);
  outputtree->SetBranchStatus("Lepton_electronIdx", 0);
  outputtree->SetBranchStatus("Lepton_muonIdx", 0);
  outputtree->SetBranchStatus("Lepton_pt", 0);
  outputtree->SetBranchStatus("Lepton_eta", 0);
  outputtree->SetBranchStatus("Lepton_phi", 0);
  outputtree->SetBranchStatus("Lepton_eCorr", 0);
  outputtree->SetBranchStatus("nVetoLepton", 0);
  outputtree->SetBranchStatus("VetoLepton_pdgId", 0);
  outputtree->SetBranchStatus("VetoLepton_electronIdx", 0);
  outputtree->SetBranchStatus("VetoLepton_muonIdx", 0);
  outputtree->SetBranchStatus("VetoLepton_pt", 0);
  outputtree->SetBranchStatus("VetoLepton_eta", 0);
  outputtree->SetBranchStatus("VetoLepton_phi", 0);
  outputtree->SetBranchStatus("VetoLepton_eCorr", 0);
  outputtree->SetBranchStatus("nCleanJet", 0);
  outputtree->SetBranchStatus("CleanJet_jetIdx", 0);
  outputtree->SetBranchStatus("CleanJet_pt", 0);
  outputtree->SetBranchStatus("CleanJet_eta", 0);
  outputtree->SetBranchStatus("CleanJet_phi", 0);
  outputtree->SetBranchStatus("Lepton_isLoose", 0);
  outputtree->SetBranchStatus("Lepton_isVeto", 0);
  outputtree->SetBranchStatus("Lepton_isWgs", 0);
  outputtree->SetBranchStatus("dmZll_veto", 0);
  //outputtree->SetBranchStatus("Lepton_isTightElectron_mvaFall17Iso_WP90", 0);
  //outputtree->SetBranchStatus("Lepton_isTightElectron_mvaFall17Iso_WP90_SS", 0);
  outputtree->SetBranchStatus("Lepton_isTightMuon_cut_Tight80x_HWWW", 0);
  outputtree->SetBranchStatus("puWeight", 0);
  outputtree->SetBranchStatus("puWeightUp", 0);
  outputtree->SetBranchStatus("puWeightDown", 0);
  outputtree->SetBranchStatus("mll", 0);
  outputtree->SetBranchStatus("dphill", 0);
  outputtree->SetBranchStatus("yll", 0);
  outputtree->SetBranchStatus("ptll", 0);
  outputtree->SetBranchStatus("pt1", 0);
  outputtree->SetBranchStatus("pt2", 0);
  outputtree->SetBranchStatus("mth", 0);
  outputtree->SetBranchStatus("mcoll", 0);
  outputtree->SetBranchStatus("mcollWW", 0);
  outputtree->SetBranchStatus("mTi", 0);
  outputtree->SetBranchStatus("mTe", 0);
  outputtree->SetBranchStatus("choiMass", 0);
  outputtree->SetBranchStatus("mR", 0);
  outputtree->SetBranchStatus("mT2", 0);
  outputtree->SetBranchStatus("channel", 0);
  outputtree->SetBranchStatus("drll", 0);
  outputtree->SetBranchStatus("dphilljet", 0);
  outputtree->SetBranchStatus("dphilljetjet", 0);
  outputtree->SetBranchStatus("dphilljetjet_cut", 0);
  outputtree->SetBranchStatus("dphillmet", 0);
  outputtree->SetBranchStatus("dphilmet", 0);
  outputtree->SetBranchStatus("dphilmet1", 0);
  outputtree->SetBranchStatus("dphilmet2", 0);
  outputtree->SetBranchStatus("mtw1", 0);
  outputtree->SetBranchStatus("mtw2", 0);
  outputtree->SetBranchStatus("mjj", 0);
  outputtree->SetBranchStatus("detajj", 0);
  outputtree->SetBranchStatus("njet", 0);
  outputtree->SetBranchStatus("mllWgSt", 0);
  outputtree->SetBranchStatus("drllWgSt", 0);
  outputtree->SetBranchStatus("mllThird", 0);
  outputtree->SetBranchStatus("mllOneThree", 0);
  outputtree->SetBranchStatus("mllTwoThree", 0);
  outputtree->SetBranchStatus("drllOneThree", 0);
  outputtree->SetBranchStatus("drllTwoThree", 0);
  outputtree->SetBranchStatus("dphijet1met", 0);
  outputtree->SetBranchStatus("dphijet2met", 0);
  outputtree->SetBranchStatus("dphijjmet", 0);
  outputtree->SetBranchStatus("dphijjmet_cut", 0);
  outputtree->SetBranchStatus("dphilep1jet1", 0);
  outputtree->SetBranchStatus("dphilep1jet2", 0);
  outputtree->SetBranchStatus("dphilep2jet1", 0);
  outputtree->SetBranchStatus("dphilep2jet2", 0);
  outputtree->SetBranchStatus("ht", 0);
  outputtree->SetBranchStatus("vht_pt", 0);
  outputtree->SetBranchStatus("vht_phi", 0);
  outputtree->SetBranchStatus("projpfmet", 0);
  outputtree->SetBranchStatus("dphiltkmet", 0);
  outputtree->SetBranchStatus("projtkmet", 0);
  outputtree->SetBranchStatus("mpmet", 0);
  outputtree->SetBranchStatus("pTWW", 0);
  outputtree->SetBranchStatus("recoil", 0);
  outputtree->SetBranchStatus("jetpt1_cut", 0);
  outputtree->SetBranchStatus("jetpt2_cut", 0);
  outputtree->SetBranchStatus("dphilljet_cut", 0);
  outputtree->SetBranchStatus("dphijet1met_cut", 0);
  outputtree->SetBranchStatus("dphijet2met_cut", 0);
  outputtree->SetBranchStatus("PfMetDivSumMet", 0);
  outputtree->SetBranchStatus("upara", 0);
  outputtree->SetBranchStatus("uperp", 0);
  outputtree->SetBranchStatus("m2ljj20", 0);
  outputtree->SetBranchStatus("m2ljj30", 0);
  outputtree->SetBranchStatus("mllmin3l", 0);
  outputtree->SetBranchStatus("zveto_3l", 0);
  outputtree->SetBranchStatus("pt3", 0);
  outputtree->SetBranchStatus("eta1", 0);
  outputtree->SetBranchStatus("eta2", 0);
  outputtree->SetBranchStatus("eta3", 0);
  outputtree->SetBranchStatus("phi1", 0);
  outputtree->SetBranchStatus("phi2", 0);
  outputtree->SetBranchStatus("phi3", 0);
  outputtree->SetBranchStatus("drllmin3l", 0);
  outputtree->SetBranchStatus("njet_3l", 0);
  outputtree->SetBranchStatus("nbjet_3l", 0);
  outputtree->SetBranchStatus("chlll", 0);
  outputtree->SetBranchStatus("pfmet", 0);
  outputtree->SetBranchStatus("mlll", 0);
  outputtree->SetBranchStatus("flagOSSF", 0);
  outputtree->SetBranchStatus("mtwww", 0);
  outputtree->SetBranchStatus("mtw1_wh3l", 0);
  outputtree->SetBranchStatus("mtw2_wh3l", 0);
  outputtree->SetBranchStatus("mtw3_wh3l", 0);
  outputtree->SetBranchStatus("minmtw_wh3l", 0);
  outputtree->SetBranchStatus("mindphi_lmet", 0);
  outputtree->SetBranchStatus("dphilllmet", 0);
  outputtree->SetBranchStatus("ptlll", 0);
  outputtree->SetBranchStatus("pTWWW", 0);
  outputtree->SetBranchStatus("dphilmet1_wh3l", 0);
  outputtree->SetBranchStatus("dphilmet2_wh3l", 0);
  outputtree->SetBranchStatus("dphilmet3_wh3l", 0);
  outputtree->SetBranchStatus("ptbest", 0);
  outputtree->SetBranchStatus("pfmetPhi_zh4l", 0);
  outputtree->SetBranchStatus("z0Mass_zh4l", 0);
  outputtree->SetBranchStatus("z1Mass_zh4l", 0);
  outputtree->SetBranchStatus("zaMass_zh4l", 0);
  outputtree->SetBranchStatus("zbMass_zh4l", 0);
  outputtree->SetBranchStatus("flagZ1SF_zh4l", 0);
  outputtree->SetBranchStatus("z0DeltaPhi_zh4l", 0);
  outputtree->SetBranchStatus("z1DeltaPhi_zh4l", 0);
  outputtree->SetBranchStatus("zaDeltaPhi_zh4l", 0);
  outputtree->SetBranchStatus("zbDeltaPhi_zh4l", 0);
  outputtree->SetBranchStatus("minDeltaPhi_zh4l", 0);
  outputtree->SetBranchStatus("z0DeltaR_zh4l", 0);
  outputtree->SetBranchStatus("z1DeltaR_zh4l", 0);
  outputtree->SetBranchStatus("zaDeltaR_zh4l", 0);
  outputtree->SetBranchStatus("zbDeltaR_zh4l", 0);
  outputtree->SetBranchStatus("lep1Mt_zh4l", 0);
  outputtree->SetBranchStatus("lep2Mt_zh4l", 0);
  outputtree->SetBranchStatus("lep3Mt_zh4l", 0);
  outputtree->SetBranchStatus("lep4Mt_zh4l", 0);
  outputtree->SetBranchStatus("minMt_zh4l", 0);
  outputtree->SetBranchStatus("z1Mt_zh4l", 0);
  outputtree->SetBranchStatus("mllll_zh4l", 0);
  outputtree->SetBranchStatus("chllll_zh4l", 0);
  outputtree->SetBranchStatus("Jet_btagSF", 0);
  outputtree->SetBranchStatus("Jet_btagSF_up", 0);
  outputtree->SetBranchStatus("Jet_btagSF_down", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_jes", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_jes", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_lf", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_lf", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_hf", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_hf", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_hfstats1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_hfstats1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_hfstats2", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_hfstats2", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_lfstats1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_lfstats1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_lfstats2", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_lfstats2", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_cferr1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_cferr1", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_up_cferr2", 0);
  outputtree->SetBranchStatus("Jet_btagSF_shape_down_cferr2", 0);
  outputtree->SetBranchStatus("btagWeight", 0);
  outputtree->SetBranchStatus("btagWeight_up_jes", 0);
  outputtree->SetBranchStatus("btagWeight_down_jes", 0);
  outputtree->SetBranchStatus("btagWeight_up_lf", 0);
  outputtree->SetBranchStatus("btagWeight_down_lf", 0);
  outputtree->SetBranchStatus("btagWeight_up_hf", 0);
  outputtree->SetBranchStatus("btagWeight_down_hf", 0);
  outputtree->SetBranchStatus("btagWeight_up_hfstats1", 0);
  outputtree->SetBranchStatus("btagWeight_down_hfstats1", 0);
  outputtree->SetBranchStatus("btagWeight_up_hfstats2", 0);
  outputtree->SetBranchStatus("btagWeight_down_hfstats2", 0);
  outputtree->SetBranchStatus("btagWeight_up_lfstats1", 0);
  outputtree->SetBranchStatus("btagWeight_down_lfstats1", 0);
  outputtree->SetBranchStatus("btagWeight_up_lfstats2", 0);
  outputtree->SetBranchStatus("btagWeight_down_lfstats2", 0);
  outputtree->SetBranchStatus("btagWeight_up_cferr1", 0);
  outputtree->SetBranchStatus("btagWeight_down_cferr1", 0);
  outputtree->SetBranchStatus("btagWeight_up_cferr2", 0);
  outputtree->SetBranchStatus("btagWeight_down_cferr2", 0);
  outputtree->SetBranchStatus("baseW", 0);
  outputtree->SetBranchStatus("Xsec", 0);
  /*outputtree->SetBranchStatus("TriggerEmulator", 0);
  outputtree->SetBranchStatus("EMTFbug_veto", 0);
  outputtree->SetBranchStatus("run_period", 0);
  outputtree->SetBranchStatus("Trigger_sngEl", 0);
  outputtree->SetBranchStatus("Trigger_sngMu", 0);
  outputtree->SetBranchStatus("Trigger_dblEl", 0);
  outputtree->SetBranchStatus("Trigger_dblMu", 0);
  outputtree->SetBranchStatus("Trigger_ElMu", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_2l", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_2l_u", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_2l_d", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_3l", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_3l_u", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_3l_d", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_4l", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_4l_u", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_4l_d", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_sngEl", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_sngMu", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_dblEl", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_dblMu", 0);
  outputtree->SetBranchStatus("TriggerEffWeight_ElMu", 0);
  outputtree->SetBranchStatus("LepCut3l", 0);
  outputtree->SetBranchStatus("METFilter_DATA", 0);
  outputtree->SetBranchStatus("LepCut2l", 0);
  outputtree->SetBranchStatus("LepCut4l", 0);
  outputtree->SetBranchStatus("LepCut2lSS", 0);*/
  outputtree->SetBranchStatus("run", 0);
  outputtree->SetBranchStatus("luminosityBlock", 0);
  outputtree->SetBranchStatus("event", 0);
  outputtree->SetBranchStatus("CaloMET_phi", 0);
  outputtree->SetBranchStatus("CaloMET_pt", 0);
  outputtree->SetBranchStatus("CaloMET_sumEt", 0);
  outputtree->SetBranchStatus("nElectron", 0);
  outputtree->SetBranchStatus("Electron_deltaEtaSC", 0);
  outputtree->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 0);
  outputtree->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 0);
  outputtree->SetBranchStatus("Electron_dr03TkSumPt", 0);
  outputtree->SetBranchStatus("Electron_dxy", 0);
  outputtree->SetBranchStatus("Electron_dxyErr", 0);
  outputtree->SetBranchStatus("Electron_dz", 0);
  outputtree->SetBranchStatus("Electron_dzErr", 0);
  outputtree->SetBranchStatus("Electron_eCorr", 0);
  outputtree->SetBranchStatus("Electron_eInvMinusPInv", 0);
  outputtree->SetBranchStatus("Electron_energyErr", 0);
  outputtree->SetBranchStatus("Electron_eta", 0);
  outputtree->SetBranchStatus("Electron_hoe", 0);
  outputtree->SetBranchStatus("Electron_ip3d", 0);
  outputtree->SetBranchStatus("Electron_mass", 0);
  outputtree->SetBranchStatus("Electron_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("Electron_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso", 0);
  outputtree->SetBranchStatus("Electron_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("Electron_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("Electron_phi", 0);
  outputtree->SetBranchStatus("Electron_pt", 0);
  outputtree->SetBranchStatus("Electron_r9", 0);
  outputtree->SetBranchStatus("Electron_sieie", 0);
  outputtree->SetBranchStatus("Electron_sip3d", 0);
  outputtree->SetBranchStatus("Electron_mvaTTH", 0);
  outputtree->SetBranchStatus("Electron_charge", 0);
  outputtree->SetBranchStatus("Electron_cutBased", 0);
  outputtree->SetBranchStatus("Electron_jetIdx", 0);
  outputtree->SetBranchStatus("Electron_pdgId", 0);
  outputtree->SetBranchStatus("Electron_photonIdx", 0);
  outputtree->SetBranchStatus("Electron_tightCharge", 0);
  outputtree->SetBranchStatus("Electron_vidNestedWPBitmap", 0);
  outputtree->SetBranchStatus("Electron_convVeto", 0);
  outputtree->SetBranchStatus("Electron_cutBased_HEEP", 0);
  outputtree->SetBranchStatus("Electron_isPFcand", 0);
  outputtree->SetBranchStatus("Electron_lostHits", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP80", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WP90", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17Iso_WPL", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP80", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WP90", 0);
  outputtree->SetBranchStatus("Electron_mvaFall17noIso_WPL", 0);
  outputtree->SetBranchStatus("nFatJet", 0);
  outputtree->SetBranchStatus("FatJet_area", 0);
  outputtree->SetBranchStatus("FatJet_btagCMVA", 0);
  outputtree->SetBranchStatus("FatJet_btagCSVV2", 0);
  outputtree->SetBranchStatus("FatJet_btagDeepB", 0);
  outputtree->SetBranchStatus("FatJet_btagHbb", 0);
  outputtree->SetBranchStatus("FatJet_eta", 0);
  outputtree->SetBranchStatus("FatJet_mass", 0);
  outputtree->SetBranchStatus("FatJet_msoftdrop", 0);
  outputtree->SetBranchStatus("FatJet_n2b1", 0);
  outputtree->SetBranchStatus("FatJet_n3b1", 0);
  outputtree->SetBranchStatus("FatJet_phi", 0);
  outputtree->SetBranchStatus("FatJet_pt", 0);
  outputtree->SetBranchStatus("FatJet_tau1", 0);
  outputtree->SetBranchStatus("FatJet_tau2", 0);
  outputtree->SetBranchStatus("FatJet_tau3", 0);
  outputtree->SetBranchStatus("FatJet_tau4", 0);
  outputtree->SetBranchStatus("FatJet_jetId", 0);
  outputtree->SetBranchStatus("FatJet_subJetIdx1", 0);
  outputtree->SetBranchStatus("FatJet_subJetIdx2", 0);
  outputtree->SetBranchStatus("nGenJetAK8", 0);
  outputtree->SetBranchStatus("GenJetAK8_eta", 0);
  outputtree->SetBranchStatus("GenJetAK8_mass", 0);
  outputtree->SetBranchStatus("GenJetAK8_phi", 0);
  outputtree->SetBranchStatus("GenJetAK8_pt", 0);
  outputtree->SetBranchStatus("nGenJet", 0);
  outputtree->SetBranchStatus("GenJet_eta", 0);
  outputtree->SetBranchStatus("GenJet_mass", 0);
  outputtree->SetBranchStatus("GenJet_phi", 0);
  outputtree->SetBranchStatus("GenJet_pt", 0);
  outputtree->SetBranchStatus("nGenPart", 0);
  outputtree->SetBranchStatus("GenPart_eta", 0);
  outputtree->SetBranchStatus("GenPart_mass", 0);
  outputtree->SetBranchStatus("GenPart_phi", 0);
  outputtree->SetBranchStatus("GenPart_pt", 0);
  outputtree->SetBranchStatus("GenPart_genPartIdxMother", 0);
  outputtree->SetBranchStatus("GenPart_pdgId", 0);
  outputtree->SetBranchStatus("GenPart_status", 0);
  outputtree->SetBranchStatus("GenPart_statusFlags", 0);
  outputtree->SetBranchStatus("nSubGenJetAK8", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_eta", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_mass", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_phi", 0);
  outputtree->SetBranchStatus("SubGenJetAK8_pt", 0);
  outputtree->SetBranchStatus("Generator_binvar", 0);
  outputtree->SetBranchStatus("Generator_scalePDF", 0);
  outputtree->SetBranchStatus("Generator_weight", 0);
  outputtree->SetBranchStatus("Generator_x1", 0);
  outputtree->SetBranchStatus("Generator_x2", 0);
  outputtree->SetBranchStatus("Generator_xpdf1", 0);
  outputtree->SetBranchStatus("Generator_xpdf2", 0);
  outputtree->SetBranchStatus("Generator_id1", 0);
  outputtree->SetBranchStatus("Generator_id2", 0);
  outputtree->SetBranchStatus("nGenVisTau", 0);
  outputtree->SetBranchStatus("GenVisTau_eta", 0);
  outputtree->SetBranchStatus("GenVisTau_mass", 0);
  outputtree->SetBranchStatus("GenVisTau_phi", 0);
  outputtree->SetBranchStatus("GenVisTau_pt", 0);
  outputtree->SetBranchStatus("GenVisTau_charge", 0);
  outputtree->SetBranchStatus("GenVisTau_genPartIdxMother", 0);
  outputtree->SetBranchStatus("GenVisTau_status", 0);
  outputtree->SetBranchStatus("genWeight", 0);
  outputtree->SetBranchStatus("LHEWeight_originalXWGTUP", 0);
  outputtree->SetBranchStatus("nLHEPdfWeight", 0);
  outputtree->SetBranchStatus("LHEPdfWeight", 0);
  outputtree->SetBranchStatus("nLHEScaleWeight", 0);
  outputtree->SetBranchStatus("LHEScaleWeight", 0);
  outputtree->SetBranchStatus("nIsoTrack", 0);
  outputtree->SetBranchStatus("IsoTrack_dxy", 0);
  outputtree->SetBranchStatus("IsoTrack_dz", 0);
  outputtree->SetBranchStatus("IsoTrack_eta", 0);
  outputtree->SetBranchStatus("IsoTrack_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("IsoTrack_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("IsoTrack_phi", 0);
  outputtree->SetBranchStatus("IsoTrack_pt", 0);
  outputtree->SetBranchStatus("IsoTrack_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("IsoTrack_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("IsoTrack_pdgId", 0);
  outputtree->SetBranchStatus("IsoTrack_isHighPurityTrack", 0);
  outputtree->SetBranchStatus("IsoTrack_isPFcand", 0);
  outputtree->SetBranchStatus("nJet", 0);
  outputtree->SetBranchStatus("Jet_area", 0);
  outputtree->SetBranchStatus("Jet_btagCMVA", 0);
  outputtree->SetBranchStatus("Jet_btagCSVV2", 0);
  outputtree->SetBranchStatus("Jet_btagDeepB", 0);
  outputtree->SetBranchStatus("Jet_btagDeepC", 0);
  outputtree->SetBranchStatus("Jet_btagDeepFlavB", 0);
  outputtree->SetBranchStatus("Jet_chEmEF", 0);
  outputtree->SetBranchStatus("Jet_chHEF", 0);
  outputtree->SetBranchStatus("Jet_eta", 0);
  outputtree->SetBranchStatus("Jet_mass", 0);
  outputtree->SetBranchStatus("Jet_neEmEF", 0);
  outputtree->SetBranchStatus("Jet_neHEF", 0);
  outputtree->SetBranchStatus("Jet_phi", 0);
  outputtree->SetBranchStatus("Jet_pt", 0);
  outputtree->SetBranchStatus("Jet_qgl", 0);
  outputtree->SetBranchStatus("Jet_rawFactor", 0);
  outputtree->SetBranchStatus("Jet_bReg", 0);
  outputtree->SetBranchStatus("Jet_electronIdx1", 0);
  outputtree->SetBranchStatus("Jet_electronIdx2", 0);
  outputtree->SetBranchStatus("Jet_jetId", 0);
  outputtree->SetBranchStatus("Jet_muonIdx1", 0);
  outputtree->SetBranchStatus("Jet_muonIdx2", 0);
  outputtree->SetBranchStatus("Jet_nConstituents", 0);
  outputtree->SetBranchStatus("Jet_nElectrons", 0);
  outputtree->SetBranchStatus("Jet_nMuons", 0);
  outputtree->SetBranchStatus("Jet_puId", 0);
  outputtree->SetBranchStatus("LHE_HT", 0);
  outputtree->SetBranchStatus("LHE_HTIncoming", 0);
  outputtree->SetBranchStatus("LHE_Vpt", 0);
  outputtree->SetBranchStatus("LHE_Njets", 0);
  outputtree->SetBranchStatus("LHE_Nb", 0);
  outputtree->SetBranchStatus("LHE_Nc", 0);
  outputtree->SetBranchStatus("LHE_Nuds", 0);
  outputtree->SetBranchStatus("LHE_Nglu", 0);
  outputtree->SetBranchStatus("LHE_NpNLO", 0);
  outputtree->SetBranchStatus("LHE_NpLO", 0);
  outputtree->SetBranchStatus("nLHEPart", 0);
  outputtree->SetBranchStatus("LHEPart_pt", 0);
  outputtree->SetBranchStatus("LHEPart_eta", 0);
  outputtree->SetBranchStatus("LHEPart_phi", 0);
  outputtree->SetBranchStatus("LHEPart_mass", 0);
  outputtree->SetBranchStatus("LHEPart_pdgId", 0);
  outputtree->SetBranchStatus("GenMET_phi", 0);
  outputtree->SetBranchStatus("GenMET_pt", 0);
  outputtree->SetBranchStatus("MET_MetUnclustEnUpDeltaX", 0);
  outputtree->SetBranchStatus("MET_MetUnclustEnUpDeltaY", 0);
  outputtree->SetBranchStatus("MET_covXX", 0);
  outputtree->SetBranchStatus("MET_covXY", 0);
  outputtree->SetBranchStatus("MET_covYY", 0);
  outputtree->SetBranchStatus("MET_phi", 0);
  outputtree->SetBranchStatus("MET_pt", 0);
  outputtree->SetBranchStatus("MET_significance", 0);
  outputtree->SetBranchStatus("MET_sumEt", 0);
  outputtree->SetBranchStatus("nMuon", 0);
  outputtree->SetBranchStatus("Muon_dxy", 0);
  outputtree->SetBranchStatus("Muon_dxyErr", 0);
  outputtree->SetBranchStatus("Muon_dz", 0);
  outputtree->SetBranchStatus("Muon_dzErr", 0);
  outputtree->SetBranchStatus("Muon_eta", 0);
  outputtree->SetBranchStatus("Muon_ip3d", 0);
  outputtree->SetBranchStatus("Muon_mass", 0);
  outputtree->SetBranchStatus("Muon_miniPFRelIso_all", 0);
  outputtree->SetBranchStatus("Muon_miniPFRelIso_chg", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso03_all", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso03_chg", 0);
  outputtree->SetBranchStatus("Muon_pfRelIso04_all", 0);
  outputtree->SetBranchStatus("Muon_phi", 0);
  outputtree->SetBranchStatus("Muon_pt", 0);
  outputtree->SetBranchStatus("Muon_ptErr", 0);
  outputtree->SetBranchStatus("Muon_segmentComp", 0);
  outputtree->SetBranchStatus("Muon_sip3d", 0);
  outputtree->SetBranchStatus("Muon_mvaTTH", 0);
  outputtree->SetBranchStatus("Muon_charge", 0);
  outputtree->SetBranchStatus("Muon_jetIdx", 0);
  outputtree->SetBranchStatus("Muon_nStations", 0);
  outputtree->SetBranchStatus("Muon_nTrackerLayers", 0);
  outputtree->SetBranchStatus("Muon_pdgId", 0);
  outputtree->SetBranchStatus("Muon_tightCharge", 0);
}
#endif // #ifdef SkimMeBaby_cxx

