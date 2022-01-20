#ifndef OfflineHelper_h
#define OfflineHelper_h
/** \class OfflineProducerHelper
 *
 *  Collection of functions to analyze HTauTau primary trees, without CMSSW dependencies.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */
 
#include <TString.h>
#include "TMath.h"
#include "HTauTauTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <utility>

using namespace std;
 
class OfflineProducerHelper {
 public:
  enum particleType {
  MUON = 0,
  ELECTRON = 1,
  TAU =2
  };

  enum pairType {
    MuHad  = 0,
    EHad   = 1,
    HadHad = 2,
    MuMu   = 3,
    EE     = 4,
    EMu    = 5,
    EEPrompt = 6, // prompt Z->ee/mumu decays
    MuMuPrompt = 7,
    Other  = 8 // for e.g. h->bb
  };

  OfflineProducerHelper();
  int FindTriggerNumber(TString triggername);
  bool IsTriggerFired(int triggerbit, int triggerNumber);
  bool IsTriggerFired(int triggerbit, TString triggerName){return IsTriggerFired(triggerbit, FindTriggerNumber(triggerName));}
  int printFiredPaths(int);
  bool isMuon(int type){if(type == MUON)return true; else return false;}
  bool isElectron(int type){if(type == ELECTRON)return true; else return false;}
  bool isTau(int type){if(type == TAU)return true; else return false;}
  int getPairType (int type1, int type2); // return pair type
  bool checkBit (int word, int bitpos); // check bit "bitpos" in a word
  
  // whatApply: use "OSCharge" (appplies on pairs only)
  // whatApply: use "All", "Iso", "pTMin", "etaMax", "againstEle", "againstMu", "Vertex"; separate various arguments with a semicolon
  // is contains "All" it will override all the other settings; additional parameters are not considered (have no effect) 
  // a selection is applied by default if no parameter is specified
  bool pairPassBaseline (HTauTauTree* tree, int iPair, TString whatApply = "All");
  bool eleBaseline (HTauTauTree* tree, int iDau, float ptMin, float relIso, TString whatApply = "All"); // return true if leptons passes the baseline selections
  bool muBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, float relIso, TString whatApply = "All");
  bool tauBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, TString whatApply = "All");
  bool tightEleMVAID (float BDT, float fSCeta); // compute tight ele MVA id WP, but isBDT in ntuples has been fixed --> this will be soon deprecated
                                                // approx!! I call it using lepton eta and not superCluster eta
  
  /*
  // separate check of various requirements applied in baseline
  bool combRelIso (HTauTauTree* tree, int iDau, float iso); // for mu, ele
  bool combIsodBetaRaw3Hits (HTauTauTree* tree, int iDau, float iso); // for tau
  bool passEleID (HTauTauTree* tree, int iDau); // tight ele ID
  bool passMuID (HTauTauTree* tree, int iDau); // medium mu ID
  bool passTauID (HTauTauTree* tree, int iDau); // tau pog ID (decay Mode Finding OLD || NEW)
  bool passTauAntiEle (HTauTauTree* tree, int iDau, int againstEleWP);
  bool passTauAntiMu  (HTauTauTree* tree, int iDau, int againstMuWP); 
  */

  int MCHiggsTauTauDecayMode (HTauTauTree* tree); // find the MC decay of the Higgs to tau in the event

  TLorentzVector buildDauP4 (HTauTauTree* tree, int iDau); // build daughter 4 vector
  TLorentzVector buildMothP4 (HTauTauTree* tree, int iMoth); // build pair 4 vector
  bool getBestJets (HTauTauTree* tree, int& jet1, int& jet2, int strategy); // select jets, possibly two b jets, returns true if found, else false
  
  ~OfflineProducerHelper(){}

 private:
  static const int nTriggers =19;
  TString triggerlist[nTriggers];
  


};

OfflineProducerHelper::OfflineProducerHelper(){
  TString tmptrigger[nTriggers]={
    "IsoMu17_eta2p1_LooseIsoPFTau20",
    "IsoMu17_eta2p1",
    "IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg",
    "IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1",
    "IsoMu24_eta2p1_IterTrk01",
    "IsoMu24_eta2p1_IterTrk02",
    "IsoMu24_eta2p1_IterTrk02_LooseIsoPFTau20",
    "Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "Ele32_eta2p1_WP85_Gsf",
    "Ele32_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "LooseIsoPFTau50_Trk30_eta2p1_MET120",
    "IsoMu16_eta2p1_CaloMET30_LooseIsoPFTau50_Trk30_eta2p1",
    "IsoMu16_eta2p1_CaloMET30",
    "Mu16_eta2p1_CaloMET30",
    "LooseIsoPFTau50_Trk30_eta2p1",
    "DoubleIsoMu17_eta2p1",
    "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
    "Ele27_eta2p1_WP85_Gsf_LooseIsoPFTau20",
    "Ele27_eta2p1_WP85_Gsf"
  };
  for(int i=0;i<nTriggers;i++){
    triggerlist[i]=tmptrigger[i];
    triggerlist[i].Prepend("HLT_");
    triggerlist[i].Append("_v1");
  }
}

int OfflineProducerHelper::FindTriggerNumber(TString triggername){
  for(int it=0;it<nTriggers;it++){ 	
  	if(triggerlist[it].CompareTo(triggername.Data())==0)return it;
  	else {
  	    TString newName=triggername.Data();
  	    newName.Prepend("HLT_");
  	    newName.Append("_v1");
  	    if(triggerlist[it].CompareTo(newName.Data())==0)return it;
  	}
  }
  return -1;
}

bool OfflineProducerHelper::IsTriggerFired(int triggerbit, int triggernumber){ 
  if(triggernumber>=0 && triggernumber<nTriggers) return triggerbit & (1 << triggernumber);
  return false;
}

int OfflineProducerHelper::printFiredPaths(int triggerbit){
  int nFired = 0;
  for(int it=0;it<nTriggers;it++){ 	
  	if(IsTriggerFired(triggerbit,it)) {
  	  printf("%s\n",triggerlist[it].Data());
  	  nFired++;
  	  }
  }
  return nFired;
}

bool OfflineProducerHelper::checkBit (int word, int bitpos)
{
    bool res = word & (1 << bitpos);
    return res;
}

int OfflineProducerHelper::getPairType (int type1, int type2)
{
    int nmu = 0;
    int nele = 0;
    int ntau = 0;
    
    if (isMuon (type1) )     nmu++;
    if (isElectron (type1) ) nele++;
    if (isTau (type1) )      ntau++;

    if (isMuon (type2) )     nmu++;
    if (isElectron (type2) ) nele++;
    if (isTau (type2) )      ntau++;

    if (nmu == 1 && nele == 0 && ntau == 1) return (int) pairType::MuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) return (int) pairType::EHad;
    if (nmu == 0 && nele == 0 && ntau == 2) return (int) pairType::HadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) return (int) pairType::MuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) return (int) pairType::EE;
    if (nmu == 1 && nele == 1 && ntau == 0) return (int) pairType::EMu;
    
    return -1;

}

bool OfflineProducerHelper::pairPassBaseline (HTauTauTree* tree, int iPair, TString whatApply)
{
    int dau1index = tree->indexDau1->at(iPair);
    int dau2index = tree->indexDau2->at(iPair);
    int type1 = tree->particleType->at(dau1index);
    int type2 = tree->particleType->at(dau2index);
    int pairType = getPairType (type1, type2);

    bool isOS = tree->isOSCand->at(iPair);
    if ((whatApply.Contains("OScharge") || whatApply.Contains("All"))  && !isOS) return false; // do not even check the rest if requiring the charge

    // pairs are always ordered as: e mu | e tau | mu tau  (e < mu < tau)
    // if same type of particle, highest pt one is the first
    if (pairType == MuHad)
    {
        bool leg1 = muBaseline (tree, dau1index, 18., 2.1, 0.1, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 0, 1, 1.5, whatApply);
        return (leg1 && leg2);
    }

    if (pairType == EHad)
    {
        bool leg1 = eleBaseline (tree, dau1index, 23., 0.1, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 3, 0, 1.5, whatApply);
        return (leg1 && leg2);
    }

    // ordered by pT and not by most isolated, but baseline asked in sync is the same...
    if (pairType == HadHad)
    {
        bool leg1 = tauBaseline (tree, dau1index, 45., 2.1, 0, 0, 1.0, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 45., 2.1, 0, 0, 1.0, whatApply);
        return (leg1 && leg2);
    }

    if (pairType == EMu)
    {
        bool leg1 = eleBaseline (tree, dau1index, 13., 0.15, whatApply);
        bool leg2 = muBaseline (tree, dau2index, 9., 2.4, 0.15, whatApply);
        return (leg1 && leg2);
    }
    
    // e e, mu mu missing for the moment...
    if (pairType == EE) return false;
    if (pairType == MuMu) return false;
    
    return false;
        
}


bool OfflineProducerHelper::eleBaseline (HTauTauTree* tree, int iDau, float ptMin, float relIso, TString whatApply)
{ 
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
 
    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }

    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    //bool idS = checkBit (tree->daughters_iseleCUT->at(iDau), 3) || byp_idS; // 3 is TIGHT ele id CUT BASED
    bool idS = tree->daughters_iseleBDT->at(iDau) || byp_idS; // use it in ntuples produced after 11 June 2015, contains tight WP bool  
    //bool idS = tightEleMVAID (tree->discriminator->at(iDau), TMath::Abs(p4.Eta())) || byp_idS; // APPROX! Using lepton eta and not super cluster eta, discriminator contains ele BDT  
    bool isoS = (tree->combreliso->at(iDau) < relIso) || byp_isoS;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < 2.5) || byp_etaS;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
    
}

bool OfflineProducerHelper::muBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, float relIso, TString whatApply)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    int discr = tree->daughters_muonID->at(iDau);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }
        
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool idS = checkBit (discr, 2) || byp_idS; // bit 2 is MEDIUM mu id
    bool isoS = (tree->combreliso->at(iDau) < relIso) || byp_isoS;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
}

// againstEleWP: 0 = VLoose, 1 = Loose, 2 = Medium, 3 = Tight, 4 = VTight [all are MVA discr]
// againstMuWP: 0 = Loose, 1 = Tight
bool OfflineProducerHelper::tauBaseline (HTauTauTree* tree, int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, TString whatApply)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_dmfS  = false;
    bool byp_agEleS = false;
    bool byp_agMuS = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All"))
    {
      byp_vertexS = byp_dmfS = byp_agEleS = byp_agMuS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_dmfS = false; 
      if (whatApply.Contains("againstEle"))  byp_agEleS = false; 
      if (whatApply.Contains("againstMu"))   byp_agMuS = false;       
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }


    if (againstEleWP < 0 || againstEleWP > 4) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstEleWP must be between 0 and 4 --> using 0" << endl;
        againstEleWP = 0;
    } 

    if (againstMuWP < 0 || againstMuWP > 1) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstMuWP must be between 0 and 1 --> using 0" << endl;
        againstMuWP = 0;
    }
    
    int agEleVal = 0;
    int agMuVal = 0;
    
    // ag ele:
    if (againstEleWP == 0)      agEleVal = tree->daughters_againstElectronVLooseMVA5->at(iDau);
    else if (againstEleWP == 1) agEleVal = tree->daughters_againstElectronLooseMVA5->at(iDau);
    else if (againstEleWP == 2) agEleVal = tree->daughters_againstElectronMediumMVA5->at(iDau);
    else if (againstEleWP == 3) agEleVal = tree->daughters_againstElectronTightMVA5->at(iDau);
    else if (againstEleWP == 4) agEleVal = tree->daughters_againstElectronVTightMVA5->at(iDau);   

    // ag mu:
    if (againstMuWP == 0)      agMuVal = tree->daughters_againstMuonLoose3->at(iDau);
    else if (againstMuWP == 1) agMuVal = tree->daughters_againstMuonTight3->at(iDau);

    bool dmfS = (tree->daughters_decayModeFindingOldDMs->at(iDau) == 1 || tree->daughters_decayModeFindingNewDMs->at(iDau) == 1) || byp_dmfS;
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool agEleS = (agEleVal == 1) || byp_agEleS; 
    bool agMuS  = (agMuVal == 1) || byp_agMuS; 
    bool isoS = (tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(iDau) < isoRaw3Hits) || byp_isoS;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;

    bool totalS = (dmfS && vertexS && agEleS && agMuS && isoS && ptS && etaS);
    return totalS;    
}

// branch isBDT has been updated, this will soon be unnecessary (and eventually obsolete)
bool OfflineProducerHelper::tightEleMVAID (float BDT, float fSCeta)
{
    // POG_MVA_ID_Run2_NonTrig_Tight PHYS14
    if (fSCeta < 0) fSCeta *= -1;
    bool isBDT = false;
    if (fSCeta <0.8) isBDT = (BDT>0.73);
    else if (fSCeta < 1.479) isBDT = (BDT>0.57);
    else isBDT = (BDT >0.05);
    
    return isBDT;
}

TLorentzVector OfflineProducerHelper::buildDauP4 (HTauTauTree* tree, int iDau)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

TLorentzVector OfflineProducerHelper::buildMothP4 (HTauTauTree* tree, int iMoth)
{
    float px = tree->mothers_px->at(iMoth);
    float py = tree->mothers_py->at(iMoth);
    float pz = tree->mothers_pz->at(iMoth);
    float e =  tree->mothers_e->at(iMoth);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

int OfflineProducerHelper::MCHiggsTauTauDecayMode (HTauTauTree* tree)
{
    int decay = -1; // good decays go from 0 to 7, see enum
    for (int i = 0; i < tree->genpart_HZDecayMode->size(); i++)
    {
        int val = tree->genpart_HZDecayMode->at(i);
        if (val >= 0 && val <= 7)
        {
            decay = val;
            break; // I don;t expect more than 1 H to tau tau per event
        }
    }
    return decay; 
}


bool OfflineProducerHelper::getBestJets (HTauTauTree* tree, int& jet1, int& jet2, int strategy)
{
    jet1 = jet2 = -1;
    switch (strategy)
    {
        case(0): // two with highest b score
        {
            int njets = tree->bCSVscore->size();
            if (njets < 2) return false;
            std::vector<std::pair<float, int>> scores;
            for (int i = 0; i < njets; i++) scores.push_back (std::make_pair(tree->bCSVscore->at(i), i));
            std::sort (scores.begin(), scores.end()); // are sorted according to the first index, i.e. the CSV score
            jet1 = (scores.at(njets-1)).second; //leading
            jet2 = (scores.at(njets-2)).second; // subleading
            return true;                            
            break;
        }
        
        default:
        {
            std::cout << "B jet selection strategy" << strategy << " not implemented" << std::endl;
            return false;
        }
    }
}


#endif
