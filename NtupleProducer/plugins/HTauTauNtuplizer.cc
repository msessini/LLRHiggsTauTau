/** \class HTauTauNtuplizer
 *
 *  No description available.
 *
 *  $Date: 2014/10/31 10:08:19 $
 *  $Revision: 1.00 $
 *  \author G. Ortona (LLR) and L. Cadamuro (LLR)
 */

#define EPSIL 1.e-5

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>
#include <bitset>
#include <any>
//#include <XYZTLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "HiggsCPinTauDecays/TauRefit/interface/RefitVertex.h"
//#include "HiggsCPinTauDecays/TauRefit/plugins/AdvancedRefitVertexProducer.h"
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
//#include <HiggsCPinTauDecays/TauRefit/plugins/AdvancedRefitVertexProducer.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>

#include <LLRHiggsTauTau/NtupleProducer/interface/SysHelper.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CutSet.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LeptonIsoHelper.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/FinalStates.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/MCHistoryTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CPTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/PUReweight.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/VBFCandidateJetSelector.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/bitops.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/Fisher.h>
//#include <LLRHiggsTauTau/NtupleProducer/interface/HTauTauConfigHelper.h>
//#include "HZZ4lNtupleFactory.h"
#include <LLRHiggsTauTau/NtupleProducer/interface/PhotonFwd.h>
#include "LLRHiggsTauTau/NtupleProducer/interface/triggerhelper.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/PUReweight.h"
//#include "LLRHiggsTauTau/NtupleProducer/Utils/OfflineProducerHelper.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/ParticleType.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/GenFlags.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/GenHelper.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/TrackParticle.h"

#include "LLRHiggsTauTau/NtupleProducer/interface/TauDecay_CMSSW.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "LLRHiggsTauTau/NtupleProducer/interface//DataMCType.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/PDGInfo.h"


#include "Geometry/Records/interface/HcalParametersRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/GeometryObjects/interface/HcalParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "TLorentzVector.h"
#include "TMatrixTSym.h"


#include "TVectorD.h"
#include "TVector3.h"
#include "boost/functional/hash.hpp"

#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/SFProvider.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "LLRHiggsTauTau/NtupleProducer/interface/PileUp.h"

#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

namespace {
  //   bool writePhotons = false;  // Write photons in the tree.
  //bool writeJets = true;     // Write jets in the tree.
  //bool writeFatJets = true;
  //bool writeSoftLep = false;
  //bool writeL1 = true;
  bool DEBUG = true;
  //int DETAIL=1;
}

using namespace std;
using namespace edm;
using namespace reco;

// Map for JEC uncertainty sources
typedef std::map<std::string, std::unique_ptr<JetCorrectionUncertainty>> myJECMap;

// class declaration

class HTauTauNtuplizer : public edm::EDAnalyzer {
public:
  /// Constructor
  explicit HTauTauNtuplizer(const edm::ParameterSet&);

  /// Destructor
  virtual ~HTauTauNtuplizer();

private:
  //----edm control---
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void Initialize();
  Int_t FindCandIndex(const reco::Candidate&, Int_t iCand);
  //----To implement here-----
  //virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  //virtual void FillPhoton(const pat::Photon& photon);
  int FillJet(const edm::View<pat::Jet>* jet,const edm::View<pat::Jet>* jetDown, const edm::View<pat::Jet>* jetUp, const edm::Event&, edm::EventSetup const&, JetCorrectionUncertainty*, myJECMap*, myJECMap*,const reco::Candidate* cand1,const reco::Candidate* cand2,int ipair);
  void FillFatJet(const edm::View<pat::Jet>* fatjets, const edm::Event&);
  void FillSoftLeptons(const edm::View<reco::Candidate> *dauhandler, const edm::Event& event, const edm::EventSetup& setup, bool theFSR, const edm::View<pat::Jet>* jets, const BXVector<l1t::Tau>* l1taus);
  void FillGenInfo(const edm::Event&);
  void FillGenJetInfo(const edm::Event&);

  void fillMCTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  bool EleVeto(const reco::Candidate* cand);
  bool MuVeto(const reco::Candidate* cand);
  bool DiEle(const reco::Candidate* cand1,const reco::Candidate* cand2);
  bool DiMuon(const reco::Candidate* cand1,const reco::Candidate* cand2);
  TLorentzVector getVisMomentumNoLep(const reco::Candidate* genL/*,std::vector<std::vector<std::vector<double> > > MCTauandProd_p4,std::vector<std::vector<int> > MCTauandProd_pdgid*/);
  int GenMatch(const reco::Candidate* genL, const edm::Event& event);
  math::XYZTLorentzVector P4Corrected(math::XYZTLorentzVector p4, int genmatch, int DM);
  int GetMatchedGen (const reco::Candidate* genL, const edm::Event& event); // return the index of the associated gen particle in the filtered gen collection, in not existing return -1
  //int CreateFlagsWord (const pat::GenericParticle* part); // build int with each bit containing some boolean flags
  void FillL1Obj(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event); // chia
  static bool CompareLegs(const reco::Candidate *, const reco::Candidate *);
  float ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy);
  float ComputePUPPIMT (math::XYZTLorentzVector visP4,float PUPPImetPhi, float PUPPImetPt);
  float ComputeMTtot (math::XYZTLorentzVector visP41, math::XYZTLorentzVector visP42,float METx, float METy);
  static bool ComparePairsbyPt(pat::Jet i, pat::Jet j);
  static bool ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j);
  static bool CHECK_BIT(unsigned long long var, int pos);
  static bool isGoodGenParticle(const reco::GenParticle &GenPar);

  bool refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup);
  bool findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup);
  TVector3 getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
		  const reco::Track *aTrack, const GlobalPoint & aPoint);

  // ----------member data ---------------------------
  //std::map <int, int> genFlagPosMap_; // to convert from input to output enum format for H/Z decays

  int trigincr=0;
  int filtersincr=0;

  //Configs
  TString theFileName;
  bool theFSR;
  Bool_t theisMC;
  Bool_t IsEmbed;
  bool do_MCSummary_;
  bool do_MCComplete_;
  Bool_t doCPVariables;
  Bool_t computeQGVar;
  string theJECName;
  Int_t theYear;
  Int_t dataMCtype;
  //Trigger
  vector<int> indexOfPath;
  vector<string> foundPaths;
  edm::InputTag processName;
  HLTConfigProvider hltConfig_;
  //Output Objects
  //TTree *myTree;//->See from ntuplefactory in zz4l
  TH1F *hCounter;
  TH1F *hTauIDs;
  TH1F *hYear;
  triggerhelper* myTriggerHelper;
  SysHelper* mySysHelper;

  //PUReweight reweight;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBits_;
  edm::EDGetTokenT<vector<Vertex>> theVtxTag;
  edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> theSecVtxTag; //FRA
  edm::EDGetTokenT<double> theRhoTag;
  edm::EDGetTokenT<double> theRhoMiniRelIsoTag;
  edm::EDGetTokenT<double> theRhoForJERTag;
  edm::EDGetTokenT<vector<PileupSummaryInfo> > thePUTag;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > thePFCandTag;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theCandTag;
  edm::EDGetTokenT<edm::View<pat::Jet> > theJetTag;
  edm::EDGetTokenT<edm::View<pat::Jet> > theSmearedJetTag;
  edm::EDGetTokenT<edm::View<pat::Jet> > theSmearedJetDownTag;
  edm::EDGetTokenT<edm::View<pat::Jet> > theSmearedJetUpTag;
  edm::EDGetTokenT<edm::View<pat::Jet> > theFatJetTag;
  edm::EDGetTokenT<edm::ValueMap<float> > theQGTaggerTag;
  edm::EDGetTokenT<edm::View<reco::Candidate> > theLepTag;
  edm::EDGetTokenT<LHEEventProduct> theLHETag;
  edm::EDGetTokenT<GenEventInfoProduct> theGenTag;
  edm::EDGetTokenT<pat::METCollection> thePUPPIMetTag;
  edm::EDGetTokenT<math::Error<2>::type> thePUPPIMETCovTag;
  edm::EDGetTokenT<double> thePUPPIMETSignifTag;
  edm::EDGetTokenT<edm::View<pat::GenericParticle> > theGenericTag;
  edm::EDGetTokenT<edm::View<reco::GenJet> > theGenJetTag;
  edm::EDGetTokenT<edm::MergeableCounter> theTotTag;
  edm::EDGetTokenT<edm::MergeableCounter> thePassTag;
  edm::EDGetTokenT<LHEEventProduct> theLHEPTag;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotTag;
  edm::EDGetTokenT<BXVector<l1t::Tau> > theL1TauTag;
  edm::EDGetTokenT<BXVector<l1t::Jet> > theL1JetTag;
  edm::EDGetTokenT<GenLumiInfoHeader> genLumiHeaderTag;
  edm::EDGetTokenT< bool >ecalBadCalibFilterUpdate_token;
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > ThePrunedGenTag_;

  edm::EDGetTokenT<pat::TauCollection> theNewTauTag;
  edm::EDGetTokenT<edm::View<reco::Candidate> >theNewTauTagCandidate;

  edm::EDGetTokenT<RefitVertexCollection > RefitVtxBSTag;
  edm::EDGetTokenT<RefitVertexCollection > RefitVtxNoBSTag;
  edm::EDGetTokenT<double> theTauSpinnerWTEventag;
  edm::EDGetTokenT<double> theTauSpinnerWTOddtag;
  edm::EDGetTokenT<double> theTauSpinnerWTMMtag;

  //flags
  //static const int nOutVars =14;
  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  //PUReweight reweight;
  //counters
  Int_t Nevt_Gen;
  Int_t Nevt_PassTrigger;
  Int_t Npairs;
  Int_t SelectedPairs;

  //Event Output variables
  ULong64_t _indexevents;
  ULong64_t _runNumber;
  Int_t _lumi;
  Int_t _year;
  Long64_t _triggerbit;
  Int_t _metfilterbit;
  //Int_t _NBadMu;
  Bool_t _passecalBadCalibFilterUpdate;
  Float_t _prefiringweight;
  Float_t _prefiringweightup;
  Float_t _prefiringweightdown;

  Float_t _PUPPIMETsignif;
  Double_t _MC_weight;
  Float_t _nominal_wt;
  std::vector<Double_t> _TheoreticalPSUnc;
  Double_t _aMCatNLOweight;
  Float_t _lheHt;
  Int_t   _lheNOutPartons;
  Int_t   _lheNOutB;
  Int_t   _lheNOutC;
  Int_t   _DataMC_Type;
  bool Event_isRealData;
  Int_t   _PUNumInteractions;
  Float_t _PUReweight;
  Double_t _MC_weight_scale_muF0p5;
  Double_t _MC_weight_scale_muF2;
  Double_t _MC_weight_scale_muR0p5;
  Double_t _MC_weight_scale_muR2;

  Float_t _pv_x=0, _pv_y=0, _pv_z=0;
  Float_t _pvGen_x=0, _pvGen_y=0, _pvGen_z=0;
  Float_t _pvRefit_x=0, _pvRefit_y=0, _pvRefit_z=0;

  std::vector<Float_t> _RefitPVBS_x,  _RefitPVBS_y, _RefitPVBS_z;
  std::vector<vector<vector<double>>> _RefitPVBS_Cov;
  std::vector<Float_t> _RefitPVBS_xError,  _RefitPVBS_yError, _RefitPVBS_zError;
  std::vector<Float_t> _RefitPVNoBS_x,  _RefitPVNoBS_y, _RefitPVNoBS_z;
  std::vector<Float_t> _RefitPVNoBS_xError,  _RefitPVNoBS_yError, _RefitPVNoBS_zError;

  std::vector<size_t> _VertexHashNoBS1, _VertexHashNoBS2;
  std::vector<size_t> _VertexHashNoBSTracksRemovedOld1, _VertexHashNoBSTracksRemovedOld2;
  std::vector<size_t> _VertexHashBS1, _VertexHashBS2;
  std::vector<size_t> _VertexHashBSTracksRemovedOld1, _VertexHashBSTracksRemovedOld2;
  std::vector<size_t> _LeptonHash;

  std::vector<double > _pvRefit_cov;

  bool _isRefitPV=false;
  unsigned int  DataMC_Type_idx;
  // pairs
  std::vector<Float_t> _mothers_px;
  std::vector<Float_t> _mothers_py;
  std::vector<Float_t> _mothers_pz;
  std::vector<Float_t> _mothers_e;
  std::vector<Long64_t> _mothers_trgSeparateMatch; // are the two legs matched to different HLT objs?
                                                   // stored bitwise for HLT paths as done for daughters_trgMatch


  std::vector<int> Muon_charge;
  std::vector<int> Muon_trackCharge;
  std::vector<int> Muon_pdgid;
  std::vector<double> Muon_B;
  std::vector<double> Muon_M;
  std::vector<std::vector<double> > Muon_par;
  std::vector<std::vector<double> > Muon_cov;
  std::vector<TrackParticle> MuonTrack;


  std::vector<std::vector<double> > PFTau_Track_par;
  std::vector<std::vector<double> > PFTau_Track_cov;
  std::vector<int>  PFTau_Track_charge;
  std::vector<int>  PFTau_Track_pdgid;
  std::vector<double>  PFTau_Track_B;
  std::vector<double>  PFTau_Track_M;

  std::vector<int> a1_charge;
  std::vector<int> a1_pdgid;
  std::vector<double> a1_B;
  std::vector<double> a1_M;
  std::vector<std::vector<double> > a1_par;
  std::vector<std::vector<double> > a1_cov;


  // reco leptons
  //std::vector<TLorentzVector> _daughters;
  std::vector<string> _trigger_name;
  std::vector<Int_t> _trigger_accept;
  std::vector<Float_t> _daughters_px;
  std::vector<Float_t> _daughters_py;
  std::vector<Float_t> _daughters_pz;
  std::vector<Float_t> _daughters_e;

  std::vector<Float_t> _daughters_charged_px;
  std::vector<Float_t> _daughters_charged_py;
  std::vector<Float_t> _daughters_charged_pz;
  std::vector<Float_t> _daughters_charged_e;

  std::vector<Float_t> _daughters_neutral_px;
  std::vector<Float_t> _daughters_neutral_py;
  std::vector<Float_t> _daughters_neutral_pz;
  std::vector<Float_t> _daughters_neutral_e;

  std::vector<Float_t> _daughters_vx;
  std::vector<Float_t> _daughters_vy;
  std::vector<Float_t> _daughters_vz;

  std::vector<Int_t> _daughters_hasTES;
  std::vector<Float_t> _daughters_px_TauUp;
  std::vector<Float_t> _daughters_py_TauUp;
  std::vector<Float_t> _daughters_pz_TauUp;
  std::vector<Float_t> _daughters_e_TauUp;
  std::vector<Float_t> _daughters_px_TauDown;
  std::vector<Float_t> _daughters_py_TauDown;
  std::vector<Float_t> _daughters_pz_TauDown;
  std::vector<Float_t> _daughters_e_TauDown;
  std::vector<Int_t> _daughters_hasEES;
  std::vector<Float_t> _daughters_px_EleUp;
  std::vector<Float_t> _daughters_py_EleUp;
  std::vector<Float_t> _daughters_pz_EleUp;
  std::vector<Float_t> _daughters_e_EleUp;
  std::vector<Float_t> _daughters_px_EleDown;
  std::vector<Float_t> _daughters_py_EleDown;
  std::vector<Float_t> _daughters_pz_EleDown;
  std::vector<Float_t> _daughters_e_EleDown;
  std::vector<Int_t> _daughters_genindex;
  std::vector<Int_t> _daughters_charge;
  std::vector<Int_t> _daughters_isTauMatched;

  std::vector<const reco::Candidate*> _softLeptons;

  std::vector<Float_t> _genpart_px;
  std::vector<Float_t> _genpart_py;
  std::vector<Float_t> _genpart_pz;
  std::vector<Float_t> _genpart_e;

  std::vector<Float_t> _genpart_pca_x;
  std::vector<Float_t> _genpart_pca_y;
  std::vector<Float_t> _genpart_pca_z;

  std::vector<Int_t> _genpart_pdg;
  std::vector<Int_t> _genpart_status;

  // Signal particles Z, W, H0, Hpm
  std::vector<std::vector<double> > MCSignalParticle_p4;
  std::vector<int> MCSignalParticle_pdgid;
  std::vector<std::vector<int> > MCSignalParticle_childpdgid;
  std::vector<int> MCSignalParticle_charge;
  std::vector<std::vector<double> > MCSignalParticle_Poca;
  std::vector<std::vector<unsigned int> > MCSignalParticle_Tauidx;

  // MC Tau Info
  std::vector<std::vector<std::vector<double> > > MCTauandProd_p4;
  std::vector<TLorentzVector > MCTauProdVisible_P4;
  std::vector<std::vector<std::vector<double> > > MCTauandProd_Vertex;
  std::vector<std::vector<int> > MCTauandProd_pdgid;
  std::vector<std::vector<unsigned int> > MCTauandProd_midx;
  std::vector<std::vector<int> > MCTauandProd_charge;
  std::vector<unsigned int> MCTau_JAK;
  std::vector<unsigned int> MCTau_DecayBitMask;
  std::vector<std::vector<float> > MC_p4;
  std::vector<int> MC_pdgid;
  std::vector<std::vector<int> > MC_childpdgid;
  std::vector<std::vector<int> > MC_childidx;
  std::vector<int> MC_charge;
  std::vector<int> MC_midx;
  std::vector<int> MC_status;

  //std::vector<Int_t> _genpart_mothInd;
  std::vector<Int_t> _genpart_HMothInd;
  std::vector<Int_t> _genpart_MSSMHMothInd;
  std::vector<Int_t> _genpart_TopMothInd;
  std::vector<Int_t> _genpart_TauMothInd;
  std::vector<Int_t> _genpart_ZMothInd;
  std::vector<Int_t> _genpart_WMothInd;
  std::vector<Int_t> _genpart_bMothInd;
  std::vector<Int_t> _genpart_HZDecayMode;
  std::vector<Int_t> _genpart_WDecayMode;
  std::vector<Int_t> _genpart_TopDecayMode;
  std::vector<Int_t> _genpart_TauGenDecayMode;
  std::vector<Int_t> _genpart_TauGenDetailedDecayMode;

  std::vector<Int_t> _genpart_flags; // vector of bit flags bout gen info

  // gen jets
  std::vector<Float_t> _genjet_px;
  std::vector<Float_t> _genjet_py;
  std::vector<Float_t> _genjet_pz;
  std::vector<Float_t> _genjet_e;
  std::vector<Int_t> _genjet_partonFlavour; // from matched pat::Jet
  std::vector<Int_t> _genjet_hadronFlavour; // (eh yes, because it is not accessible easily from gen jets)

  //L1 taus
  std::vector<Float_t> _L1_tauEt;
  std::vector<Float_t> _L1_tauEta;
  std::vector<Float_t> _L1_tauPhi;
  std::vector<short int> _L1_tauIso;

  //L1 jets
  std::vector<Float_t> _L1_jetEt;
  std::vector<Float_t> _L1_jetEta;
  std::vector<Float_t> _L1_jetPhi;

  //std::vector<math::XYZTLorentzVector> _daughter2;

  //Mothers output variables
  std::vector<Int_t> _indexDau1;
  std::vector<Int_t> _indexDau2;
  std::vector<Float_t> _daughters_HLTpt;
  //  std::vector<Bool_t>  _daughters_isL1IsoTau28Matched;
  std::vector<Float_t>  _daughters_highestEt_L1IsoTauMatched;
  //std::vector<Int_t> _genDaughters;
  std::vector<Bool_t> _isOSCand;

  std::vector<Float_t> _metx;
  std::vector<Float_t> _mety;

  std::vector<Float_t> _mTDau1;
  std::vector<Float_t> _mTDau2;

  //Leptons variables
  std::vector<Int_t> _pdgdau;
  std::vector<Int_t> _particleType;//0=muon, 1=e, 2=tau
  std::vector<Float_t> _combreliso;
  std::vector<Float_t> _combreliso03;
  //std::vector<Float_t> _discriminator;//BDT for ele, discriminator for tau,
  std::vector<Int_t> _daughters_muonID; //bitwise (bit 0 loose, 1 soft , 2 medium, 3 tight, 4 highPT 5 tight_noVtx)
  std::vector<Int_t> _daughters_typeOfMuon; //bitwise, 0=PF, 1=Global, 2=Tracker
  std::vector<Float_t> _dxy;
  std::vector<Float_t> _dz;
  std::vector<Float_t> _dxy_innerTrack;
  std::vector<Float_t> _dz_innerTrack;
  std::vector<Float_t> _daughters_rel_error_trackpt;
  std::vector<Float_t> _SIP;
  //std::vector<bool> _daughters_iseleBDT; //isBDT for ele
  std::vector<bool> _daughters_iseleWPLoose; //isBDT for ele
  std::vector<bool> _daughters_iseleWP80; //isBDT for ele
  std::vector<bool> _daughters_iseleWP90; //isBDT for ele
  std::vector<bool> _daughters_iseleNoIsoWPLoose; //isBDT for ele no Iso
  std::vector<bool> _daughters_iseleNoIsoWP80; //isBDT for ele no Iso
  std::vector<bool> _daughters_iseleNoIsoWP90; //isBDT for ele no Iso

  std::vector<bool> _daughters_iseleCutBased; //isBDT for ele no Iso
  std::vector<Float_t> _daughters_eleMVAnt; //isBDT for ele
  //std::vector<Float_t> _daughters_eleMVAntNoIso;
  std::vector<Float_t> _daughters_eleMVA_HZZ; //isBDT for ele
  std::vector<bool> _daughters_passConversionVeto; //isBDT for ele
  std::vector<int>  _daughters_eleMissingHits;
  //std::vector<int>  _daughters_eleMissingLostHits;
  std::vector<bool>  _daughters_iseleChargeConsistent;
  //std::vector<int> _daughters_iseleCUT; //CUT ID for ele (0=veto,1=loose,2=medium,3=tight)
  std::vector<Int_t> _decayType;//for taus only
  std::vector<Int_t> _genmatch;//for taus only
  std::vector<Long64_t> _daughters_tauID; //bitwise. check h_tauID for histogram list
  static const int ntauIds = 61;
  TString tauIDStrings[ntauIds] = {
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    "againstMuonLoose3",
    "againstMuonTight3",
    "againstElectronVLooseMVA6",
    "againstElectronLooseMVA6",
    "againstElectronMediumMVA6",
    "againstElectronTightMVA6",
    "againstElectronVTightMVA6",
    "byVLooseIsolationMVArun2v1DBoldDMwLT",
    "byLooseIsolationMVArun2v1DBoldDMwLT",
    "byMediumIsolationMVArun2v1DBoldDMwLT",
    "byTightIsolationMVArun2v1DBoldDMwLT",
    "byVTightIsolationMVArun2v1DBoldDMwLT",
    "byVLooseIsolationMVArun2v1DBnewDMwLT",
    "byLooseIsolationMVArun2v1DBnewDMwLT",
    "byMediumIsolationMVArun2v1DBnewDMwLT",
    "byTightIsolationMVArun2v1DBnewDMwLT",
    "byVTightIsolationMVArun2v1DBnewDMwLT",
    "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
    "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
    "byTightIsolationMVArun2v1DBdR03oldDMwLT",
    "byVTightIsolationMVArun2v1DBdR03oldDMwLT",
    "byVLooseIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    "byVVVLooseDeepTau2017v2p1VSjet",
    "byVVLooseDeepTau2017v2p1VSjet",
    "byVLooseDeepTau2017v2p1VSjet",
    "byLooseDeepTau2017v2p1VSjet",
    "byMediumDeepTau2017v2p1VSjet",
    "byTightDeepTau2017v2p1VSjet",
    "byVTightDeepTau2017v2p1VSjet",
    "byVVTightDeepTau2017v2p1VSjet",
    "byVVVLooseDeepTau2017v2p1VSe",
    "byVVLooseDeepTau2017v2p1VSe",
    "byVLooseDeepTau2017v2p1VSe",
    "byLooseDeepTau2017v2p1VSe",
    "byMediumDeepTau2017v2p1VSe",
    "byTightDeepTau2017v2p1VSe",
    "byVTightDeepTau2017v2p1VSe",
    "byVVTightDeepTau2017v2p1VSe",
    "byVLooseDeepTau2017v2p1VSmu",
    "byLooseDeepTau2017v2p1VSmu",
    "byMediumDeepTau2017v2p1VSmu",
    "byTightDeepTau2017v2p1VSmu"
  };
  std::vector<Float_t> _daughters_IetaIeta;
  std::vector<Float_t> _daughters_full5x5_IetaIeta;
  std::vector<Float_t> _daughters_hOverE;
  std::vector<Float_t> _daughters_deltaEtaSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_deltaPhiSuperClusterTrackAtVtx;
  std::vector<Float_t> _daughters_IoEmIoP;
  std::vector<Float_t> _daughters_IoEmIoP_ttH;
  //  std::vector<Float_t> _daughters_SCeta;
  std::vector<Float_t> _daughters_depositR03_tracker;
  std::vector<Float_t> _daughters_depositR03_ecal;
  std::vector<Float_t> _daughters_depositR03_hcal;
  std::vector<Int_t> _daughters_decayModeFindingOldDMs;

  //std::vector<Int_t> _MVADM2016v1;
  std::vector<Int_t> _MVADM2017v1;
  //std::vector<Int_t> _MVADM2017v1InTauFiller;


  std::vector<Float_t> _daughters_footprintCorrection;
  std::vector<Float_t> _daughters_neutralIsoPtSumWeight;
  std::vector<Float_t> _daughters_photonPtSumOutsideSignalCone;

  std::vector<Int_t> _daughters_decayModeFindingNewDMs;
  std::vector<Float_t> _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  std::vector<Float_t> _daughters_byIsolationMVArun2v1DBoldDMwLTraw;
  std::vector<Float_t> _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017; //FRA
  std::vector<Float_t> _daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017; //FRA
  std::vector<Float_t> _daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017; //FRA
  std::vector<Int_t> _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017; //FRA
  std::vector<Float_t> _daughters_byDeepTau2017v2p1VSjetraw;
  std::vector<Float_t> _daughters_byDeepTau2017v2p1VSeraw;
  std::vector<Float_t> _daughters_byDeepTau2017v2p1VSmuraw;
  std::vector<Float_t> _daughters_chargedIsoPtSum;
  std::vector<Float_t> _daughters_neutralIsoPtSum;
  std::vector<Float_t> _daughters_puCorrPtSum;
  std::vector<Int_t> _daughters_numChargedParticlesSignalCone;
  std::vector<Int_t> _daughters_numNeutralHadronsSignalCone;
  std::vector<Int_t> _daughters_numPhotonsSignalCone;
  std::vector<Int_t> _daughters_numParticlesSignalCone;
  std::vector<Int_t> _daughters_numChargedParticlesIsoCone;
  std::vector<Int_t> _daughters_numNeutralHadronsIsoCone;
  std::vector<Int_t> _daughters_numPhotonsIsoCone;
  std::vector<Int_t> _daughters_numParticlesIsoCone;
  std::vector<Float_t> _daughters_leadChargedParticlePt;
  std::vector<Float_t> _daughters_trackRefPt;
  std::vector<Long64_t> _daughters_trgMatched;
  std::vector<Long64_t> _daughters_FilterFired;
  std::vector<Long64_t> _daughters_isGoodTriggerType;
  std::vector<Long64_t> _daughters_L3FilterFired;
  std::vector<Long64_t> _daughters_L3FilterFiredLast;


  std::vector<std::vector<double > > _PFTauSVPos;
  std::vector<std::vector<double > > _PFTauSVCov;


  std::vector<double>  PFTauGEOMFlightLenght;
  std::vector<double>  PFTauGEOMFlightLenghtSignificance;



  std::vector<std::vector<std::vector<double > > > _PFTauPionsP4;
  std::vector<std::vector<std::vector<double > > > _PFTauRefitPionsP4;
  std::vector<std::vector<double > > _PFTauPionsCharge;
  std::vector<std::vector<double > > _PFTauRefitPionsCharge;
  std::vector<float> PFTauTrack_deltaR;
  std::vector<std::vector<double > > PFTauLeadTrackLV;
  std::vector<std::vector<double > > _PFTauSVChi2NDofMatchingQuality;


  std::vector<std::vector<double> > PFTau_a1_lvp;
  std::vector<std::vector<double> > PFTau_a1_cov;
  std::vector<int>  PFTau_a1_charge;
  std::vector<int>  PFTau_a1_pdgid;
  std::vector<double>  PFTau_a1_B;
  std::vector<double>  PFTau_a1_M;
  std::vector<LorentzVectorParticle> A1LVP;

  std::vector<Float_t> TauFLSignificance;



  //  std::vector<Int_t> _daughters_jetNDauChargedMVASel;
  std::vector<Float_t> _daughters_miniRelIsoCharged;
  std::vector<Float_t> _daughters_miniRelIsoNeutral;
  std::vector<Float_t> _daughters_jetPtRel;
  std::vector<Float_t> _daughters_jetPtRatio;
  std::vector<Float_t> _daughters_jetBTagCSV;
  std::vector<Float_t> _daughters_jetBTagDeepCSV;
  std::vector<Float_t> _daughters_jetBTagDeepFlavor;

  std::vector<Float_t> _daughters_pca_x;
  std::vector<Float_t> _daughters_pca_y;
  std::vector<Float_t> _daughters_pca_z;

  std::vector<Float_t> _daughters_pcaRefitPV_x;
  std::vector<Float_t> _daughters_pcaRefitPV_y;
  std::vector<Float_t> _daughters_pcaRefitPV_z;

  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_x;
  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_y;
  std::vector<std::vector<float>> _daughters_pcaRefitPVBS_z;

  std::vector<Float_t> _daughters_pcaGenPV_x;
  std::vector<Float_t> _daughters_pcaGenPV_y;
  std::vector<Float_t> _daughters_pcaGenPV_z;


  //Jets variables
  //Int_t _numberOfJets;
  //std::vector<TLorentzVector> _jets;
  std::vector<Long64_t> _jets_VBFleadFilterMatch;    //FRA
  std::vector<Long64_t> _jets_VBFsubleadFilterMatch; //FRA
  std::vector<Float_t> _jets_px;
  std::vector<Float_t> _jets_py;
  std::vector<Float_t> _jets_pz;
  std::vector<Float_t> _jets_e;
  //  std::vector<Float_t> _jets_rawPt;
  std::vector<Float_t> _jets_area;
  std::vector<Float_t> _jets_mT;
  std::vector<Float_t> _jets_vtxPt;
  std::vector<Float_t> _jets_vtxMass;
  std::vector<Float_t> _jets_vtx3dL;
  std::vector<Float_t> _jets_vtxNtrk;
  std::vector<Float_t> _jets_vtx3deL;
  std::vector<Float_t> _jets_leadTrackPt;
  std::vector<Float_t> _jets_leptonPtRel;
  std::vector<Float_t> _jets_leptonPt;
  std::vector<Float_t> _jets_leptonDeltaR;
  std::vector<Float_t> _jets_chEmEF;
  std::vector<Float_t> _jets_chHEF;
  std::vector<Float_t> _jets_nEmEF;
  std::vector<Float_t> _jets_nHEF;
  std::vector<Float_t> _jets_MUF;
  std::vector<Int_t>   _jets_neMult;
  std::vector<Int_t>   _jets_chMult;
  std::vector<Float_t> _jets_jecUnc;
  std::vector<Float_t> _pileupjetidMVA;

  std::vector<Float_t> _jetsDown_px;
  std::vector<Float_t> _jetsDown_py;
  std::vector<Float_t> _jetsDown_pz;
  std::vector<Float_t> _jetsDown_e;
  std::vector<Float_t> _jetsDown_area;
  std::vector<Float_t> _jetsDown_mT;
  std::vector<Float_t> _jetsDown_leadTrackPt;
  std::vector<Float_t> _jetsDown_leptonPtRel;
  std::vector<Float_t> _jetsDown_leptonPt;
  std::vector<Float_t> _jetsDown_leptonDeltaR;
  std::vector<Float_t> _jetsDown_chEmEF;
  std::vector<Float_t> _jetsDown_chHEF;
  std::vector<Float_t> _jetsDown_nEmEF;
  std::vector<Float_t> _jetsDown_nHEF;
  std::vector<Float_t> _jetsDown_MUF;
  std::vector<Int_t>   _jetsDown_neMult;
  std::vector<Int_t>   _jetsDown_chMult;
  std::vector<Float_t> _pileupjetidMVADown;

  std::vector<Float_t> _jetsUp_px;
  std::vector<Float_t> _jetsUp_py;
  std::vector<Float_t> _jetsUp_pz;
  std::vector<Float_t> _jetsUp_e;
  std::vector<Float_t> _jetsUp_area;
  std::vector<Float_t> _jetsUp_mT;
  std::vector<Float_t> _jetsUp_leadTrackPt;
  std::vector<Float_t> _jetsUp_leptonPtRel;
  std::vector<Float_t> _jetsUp_leptonPt;
  std::vector<Float_t> _jetsUp_leptonDeltaR;
  std::vector<Float_t> _jetsUp_chEmEF;
  std::vector<Float_t> _jetsUp_chHEF;
  std::vector<Float_t> _jetsUp_nEmEF;
  std::vector<Float_t> _jetsUp_nHEF;
  std::vector<Float_t> _jetsUp_MUF;
  std::vector<Int_t>   _jetsUp_neMult;
  std::vector<Int_t>   _jetsUp_chMult;
  std::vector<Float_t> _pileupjetidMVAUp;

  // JEC uncertainty sources
  std::vector<Float_t> _jets_jetUnc_AbsoluteFlavMap_up; // up variations
  std::vector<Float_t> _jets_jetUnc_AbsoluteMPFBias_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteSample_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteScale_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteStat_up;
  std::vector<Float_t> _jets_jetUnc_FlavorQCD_up;
  std::vector<Float_t> _jets_jetUnc_Fragmentation_up;
  std::vector<Float_t> _jets_jetUnc_PileUpDataMC_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtBB_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC1_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC2_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtHF_up;
  std::vector<Float_t> _jets_jetUnc_PileUpPtRef_up;
  std::vector<Float_t> _jets_jetUnc_RelativeBal_up;
  std::vector<Float_t> _jets_jetUnc_RelativeFSR_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC1_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC2_up;
  std::vector<Float_t> _jets_jetUnc_RelativeJERHF_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtBB_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC1_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC2_up;
  std::vector<Float_t> _jets_jetUnc_RelativePtHF_up;
  std::vector<Float_t> _jets_jetUnc_RelativeSample_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatEC_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatFSR_up;
  std::vector<Float_t> _jets_jetUnc_RelativeStatHF_up;
  std::vector<Float_t> _jets_jetUnc_SinglePionECAL_up;
  std::vector<Float_t> _jets_jetUnc_SinglePionHCAL_up;
  std::vector<Float_t> _jets_jetUnc_TimePtEta_up;
  std::vector<Float_t> _jets_jetUnc_AbsoluteFlavMap_dw; // down variations
  std::vector<Float_t> _jets_jetUnc_AbsoluteMPFBias_dw;
  std::vector<Float_t> _jets_jetUnc_AbsoluteSample_dw;
  std::vector<Float_t> _jets_jetUnc_AbsoluteScale_dw;
  std::vector<Float_t> _jets_jetUnc_AbsoluteStat_dw;
  std::vector<Float_t> _jets_jetUnc_FlavorQCD_dw;
  std::vector<Float_t> _jets_jetUnc_Fragmentation_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpDataMC_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtBB_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC1_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtEC2_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtHF_dw;
  std::vector<Float_t> _jets_jetUnc_PileUpPtRef_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeBal_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeFSR_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC1_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJEREC2_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeJERHF_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtBB_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC1_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtEC2_dw;
  std::vector<Float_t> _jets_jetUnc_RelativePtHF_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeSample_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatEC_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatFSR_dw;
  std::vector<Float_t> _jets_jetUnc_RelativeStatHF_dw;
  std::vector<Float_t> _jets_jetUnc_SinglePionECAL_dw;
  std::vector<Float_t> _jets_jetUnc_SinglePionHCAL_dw;
  std::vector<Float_t> _jets_jetUnc_TimePtEta_dw;

  myJECMap jecSourceUncProviders;

  // JEC uncertainty sources Regrouped
  std::vector<Float_t> _jets_jetUncRegrouped_FlavorQCD_up;              // up variations
  std::vector<Float_t> _jets_jetUncRegrouped_RelativeBal_up;
  std::vector<Float_t> _jets_jetUncRegrouped_HF_up;
  std::vector<Float_t> _jets_jetUncRegrouped_BBEC1_up;
  std::vector<Float_t> _jets_jetUncRegrouped_EC2_up;
  std::vector<Float_t> _jets_jetUncRegrouped_Absolute_up;
  std::vector<Float_t> _jets_jetUncRegrouped_BBEC1_YEAR_up;
  std::vector<Float_t> _jets_jetUncRegrouped_EC2_YEAR_up;
  std::vector<Float_t> _jets_jetUncRegrouped_Absolute_YEAR_up;
  std::vector<Float_t> _jets_jetUncRegrouped_HF_YEAR_up;
  std::vector<Float_t> _jets_jetUncRegrouped_RelativeSample_YEAR_up;
  std::vector<Float_t> _jets_jetUncRegrouped_Total_up;
  std::vector<Float_t> _jets_jetUncRegrouped_FlavorQCD_dw;              // down variations
  std::vector<Float_t> _jets_jetUncRegrouped_RelativeBal_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_HF_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_BBEC1_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_EC2_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_Absolute_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_BBEC1_YEAR_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_EC2_YEAR_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_Absolute_YEAR_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_HF_YEAR_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_RelativeSample_YEAR_dw;
  std::vector<Float_t> _jets_jetUncRegrouped_Total_dw;
  myJECMap jecSourceUncRegroupedProviders;

  std::vector<std::string> m_jec_sources_2016 = {
    "AbsoluteFlavMap",
    "AbsoluteMPFBias",
    "AbsoluteScale",
    "AbsoluteStat",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeBal",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeSample",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta" };
  std::vector<std::string> m_jec_sources_2017 = {
    "AbsoluteFlavMap",
    "AbsoluteMPFBias",
    "AbsoluteScale",
    "AbsoluteStat",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeBal",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeSample",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta" };
  std::vector<std::string> m_jec_sources_2018 = {
    "AbsoluteFlavMap",
    "AbsoluteMPFBias",
    "AbsoluteSample",  //FRA: new for 2018 data
    "AbsoluteScale",
    "AbsoluteStat",
    "FlavorQCD",
    "Fragmentation",
    "PileUpDataMC",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
    "PileUpPtRef",
    "RelativeBal",
    "RelativeFSR",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeSample",
    "RelativeStatEC",
    "RelativeStatFSR",
    "RelativeStatHF",
    "SinglePionECAL",
    "SinglePionHCAL",
    "TimePtEta" };
  std::map<std::string, std::vector<Float_t>> _SourceUncVal_up;
  std::map<std::string, std::vector<Float_t>> _SourceUncVal_dw;

  std::vector<std::string> m_jec_sources_regrouped_2016 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_2016",
    "EC2_2016",
    "Absolute_2016",
    "HF_2016",
    "RelativeSample_2016",
    "Total"
  };
  std::vector<std::string> m_jec_sources_regrouped_2017 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_2017",
    "EC2_2017",
    "Absolute_2017",
    "HF_2017",
    "RelativeSample_2017",
    "Total"
  };
  std::vector<std::string> m_jec_sources_regrouped_2018 = {
    "FlavorQCD",
    "RelativeBal",
    "HF",
    "BBEC1",
    "EC2",
    "Absolute",
    "BBEC1_2018",
    "EC2_2018",
    "Absolute_2018",
    "HF_2018",
    "RelativeSample_2018",
    "Total"
  };

  std::map<std::string, std::vector<Float_t>> _SourceUncValRegrouped_up;
  std::map<std::string, std::vector<Float_t>> _SourceUncValRegrouped_dw;

  //std::map<std::string, Double_t> _TheoreticalScaleUnc;
  Double_t _TheoreticalScaleUncTab[18];

  std::vector<Float_t> _jets_QGdiscr;

  std::vector<Float_t> _ak8jets_px;
  std::vector<Float_t> _ak8jets_py;
  std::vector<Float_t> _ak8jets_pz;
  std::vector<Float_t> _ak8jets_e;
  std::vector<Float_t> _ak8jets_SoftDropMass;
  std::vector<Float_t> _ak8jets_PrunedMass;
  std::vector<Float_t> _ak8jets_TrimmedMass;
  std::vector<Float_t> _ak8jets_FilteredMass;
  std::vector<Float_t> _ak8jets_tau1; // subjettiness
  std::vector<Float_t> _ak8jets_tau2; // subjettiness
  std::vector<Float_t> _ak8jets_tau3; // subjettiness
  std::vector<Float_t> _ak8jets_tau4; // subjettiness
  std::vector<Float_t> _ak8jets_CSV; // CSV score
  std::vector<Float_t> _ak8jets_deepCSV_probb; // CSV score
  std::vector<Float_t> _ak8jets_deepCSV_probbb; // CSV score
  std::vector<Float_t> _ak8jets_deepFlavor_probb; // Flavor score
  std::vector<Float_t> _ak8jets_deepFlavor_probbb; // Flavor score
  std::vector<Float_t> _ak8jets_deepFlavor_problepb; // Flavor score
  std::vector<Int_t>   _ak8jets_nsubjets;

  // subjets of ak8 -- store ALL subjets, and link them with an idx to the ak8 jet vectors
  std::vector<Float_t> _subjets_px;
  std::vector<Float_t> _subjets_py;
  std::vector<Float_t> _subjets_pz;
  std::vector<Float_t> _subjets_e;
  std::vector<Float_t> _subjets_CSV;
  std::vector<Float_t> _subjets_deepCSV_probb;
  std::vector<Float_t> _subjets_deepCSV_probbb;
  std::vector<Float_t> _subjets_deepFlavor_probb;
  std::vector<Float_t> _subjets_deepFlavor_probbb;
  std::vector<Float_t> _subjets_deepFlavor_problepb;
  std::vector<Int_t>   _subjets_ak8MotherIdx;

  std::vector<Int_t> _jets_Flavour; // parton flavour
  std::vector<Int_t> _jets_HadronFlavour; // hadron flavour
  std::vector<Int_t> _jets_genjetIndex; // index of matched gen jet in genjet vector
  std::vector<Int_t> _jetsDown_Flavour; // parton flavour
  std::vector<Int_t> _jetsDown_HadronFlavour; // hadron flavour
  std::vector<Int_t> _jetsDown_genjetIndex; // index of matched gen jet in genjet vector
  std::vector<Int_t> _jetsUp_Flavour; // parton flavour
  std::vector<Int_t> _jetsUp_HadronFlavour; // hadron flavour
  std::vector<Int_t> _jetsUp_genjetIndex; // index of matched gen jet in genjet vector
  std::vector<Float_t> _bdiscr;
  std::vector<Float_t> _bdiscr2; //CSVv2
  std::vector<Float_t> _bdiscr3;

  std::vector<Float_t> _bdiscr4; //DeepCSV_probb
  std::vector<Float_t> _bdiscr5; //DeepCSV_probbb
  std::vector<Float_t> _bdiscr6; //DeepCSV_probudsg
  std::vector<Float_t> _bdiscr7; //DeepCSV_probc
  std::vector<Float_t> _bdiscr8; //DeepCSV_probcc

  std::vector<Float_t> _bdiscr9;  //DeepFlavor_probb
  std::vector<Float_t> _bdiscr10; //DeepFlavor_probbb
  std::vector<Float_t> _bdiscr11; //DeepFlavor_problepb
  std::vector<Float_t> _bdiscr12; //DeepFlavor_probc
  std::vector<Float_t> _bdiscr13; //DeepFlavor_probuds
  std::vector<Float_t> _bdiscr14; //DeepFlavor_probg

  std::vector<Float_t> _bdiscrDown;
  std::vector<Float_t> _bdiscr2Down; //CSVv2
  std::vector<Float_t> _bdiscr3Down;

  std::vector<Float_t> _bdiscr4Down; //DeepCSV_probb
  std::vector<Float_t> _bdiscr5Down; //DeepCSV_probbb
  std::vector<Float_t> _bdiscr6Down; //DeepCSV_probudsg
  std::vector<Float_t> _bdiscr7Down; //DeepCSV_probc
  std::vector<Float_t> _bdiscr8Down; //DeepCSV_probcc

  std::vector<Float_t> _bdiscr9Down;  //DeepFlavor_probb
  std::vector<Float_t> _bdiscr10Down; //DeepFlavor_probbb
  std::vector<Float_t> _bdiscr11Down; //DeepFlavor_problepb
  std::vector<Float_t> _bdiscr12Down; //DeepFlavor_probc
  std::vector<Float_t> _bdiscr13Down; //DeepFlavor_probuds
  std::vector<Float_t> _bdiscr14Down; //DeepFlavor_probg

  std::vector<Float_t> _bdiscrUp;
  std::vector<Float_t> _bdiscr2Up; //CSVv2
  std::vector<Float_t> _bdiscr3Up;

  std::vector<Float_t> _bdiscr4Up; //DeepCSV_probb
  std::vector<Float_t> _bdiscr5Up; //DeepCSV_probbb
  std::vector<Float_t> _bdiscr6Up; //DeepCSV_probudsg
  std::vector<Float_t> _bdiscr7Up; //DeepCSV_probc
  std::vector<Float_t> _bdiscr8Up; //DeepCSV_probcc

  std::vector<Float_t> _bdiscr9Up;  //DeepFlavor_probb
  std::vector<Float_t> _bdiscr10Up; //DeepFlavor_probbb
  std::vector<Float_t> _bdiscr11Up; //DeepFlavor_problepb
  std::vector<Float_t> _bdiscr12Up; //DeepFlavor_probc
  std::vector<Float_t> _bdiscr13Up; //DeepFlavor_probuds
  std::vector<Float_t> _bdiscr14Up; //DeepFlavor_probg
  //std::vector<Int_t> _jetID; //1=loose, 2=tight, 3=tightlepveto

  std::vector<Bool_t> _looseJetID;
  std::vector<Bool_t> _tightJetID;
  std::vector<Bool_t> _tightLepVetoJetID;
  std::vector<Float_t> _jetrawf;

  std::vector<Bool_t> _looseJetIDDown;
  std::vector<Bool_t> _tightJetIDDown;
  std::vector<Bool_t> _tightLepVetoJetIDDown;
  std::vector<Float_t> _jetrawfDown;

  std::vector<Bool_t> _looseJetIDUp;
  std::vector<Bool_t> _tightJetIDUp;
  std::vector<Bool_t> _tightLepVetoJetIDUp;
  std::vector<Float_t> _jetrawfUp;

  std::vector<Float_t> _jets_JER; // Jet Energy Resolution

  // SUSY info
  TString _susyModel;

  std::vector<std::vector<int> > vTrgMatchedToDau_idx;

  //Synchro

  //Bool_t dilepton_veto=false;
  Bool_t eleveto=false;
  Bool_t muonveto=false;
  std::vector<Bool_t> trg_singleelectron;
  std::vector<Bool_t> trg_singlemuon;
  std::vector<Bool_t> trg_singletau;
  std::vector<Bool_t> trg_muonelectron;
  std::vector<Bool_t> trg_mutaucross;
  std::vector<Bool_t> trg_doubletau;

  std::vector<Float_t> pt_1;
  std::vector<Float_t> phi_1;
  std::vector<Float_t> eta_1;
  std::vector<Float_t> m_1;
  std::vector<Float_t> q_1;
  std::vector<Float_t> d0_1;
  std::vector<Float_t> dz_1;
  std::vector<Float_t> mt_1;
  //std::vector<Float_t> pfmt_1;
  std::vector<Float_t> puppimt_1;
  std::vector<Float_t> iso_1;
  std::vector<Float_t> gen_match_1;
  std::vector<Bool_t> againstElectronLooseMVA6_1;
  std::vector<Bool_t> againstElectronMediumMVA6_1;
  std::vector<Bool_t> againstElectronTightMVA6_1;
  std::vector<Bool_t> againstElectronVLooseMVA6_1;
  std::vector<Bool_t> againstElectronVTightMVA6_1;
  std::vector<Bool_t> againstMuonLoose3_1;
  std::vector<Bool_t> againstMuonTight3_1;
  std::vector<Float_t> byIsolationMVA3oldDMwLTraw_1;
  string cmsswBase = (getenv("CMSSW_BASE"));

  TauIDSFTool *tauSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSmu","VLoose");
  TauIDSFTool *tauSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSmu","VLoose");
  TauIDSFTool *tauSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSmu","VLoose");

  std::vector<Float_t> idisoweight_1;
  std::vector<Float_t> antieweight_1;
  std::vector<Float_t> antimuweight_1;

  std::vector<Float_t> pt_2;
  std::vector<Float_t> phi_2;
  std::vector<Float_t> eta_2;
  std::vector<Float_t> m_2;
  std::vector<Float_t> q_2;
  std::vector<Float_t> d0_2;
  std::vector<Float_t> dz_2;
  std::vector<Float_t> mt_2;
  //std::vector<Float_t> pfmt_2;
  std::vector<Float_t> puppimt_2;
  std::vector<Float_t> iso_2;
  std::vector<Float_t> gen_match_2;
  std::vector<Bool_t> againstElectronLooseMVA6_2;
  std::vector<Bool_t> againstElectronMediumMVA6_2;
  std::vector<Bool_t> againstElectronTightMVA6_2;
  std::vector<Bool_t> againstElectronVLooseMVA6_2;
  std::vector<Bool_t> againstElectronVTightMVA6_2;
  std::vector<Bool_t> againstMuonLoose3_2;
  std::vector<Bool_t> againstMuonTight3_2;
  std::vector<Float_t> byIsolationMVA3oldDMwLTraw_2;

  std::vector<Float_t> idisoweight_2;
  std::vector<Float_t> antieweight_2;
  std::vector<Float_t> antimuweight_2;

  std::vector<Float_t> pt_tt;
  std::vector<Float_t> pt_vis;
  std::vector<Float_t> mt_tot;
  std::vector<Float_t> m_vis;

  Float_t _PUPPImet=-99;
  Float_t _PUPPImetphi=-99;
  std::vector<Float_t> pzetavis;
  std::vector<Float_t> pzetamiss;
  //std::vector<Float_t> pfpzetamiss;
  std::vector<Float_t> puppipzetamiss;

  Float_t _PUPPIMETCov00=-99;
  Float_t _PUPPIMETCov01=-99;
  Float_t _PUPPIMETCov10=-99;
  Float_t _PUPPIMETCov11=-99;
  Float_t _puppimet_ex_JetEnUp=-99;
  Float_t _puppimet_ey_JetEnUp=-99;
  Float_t _puppimet_ex_JetEnDown=-99;
  Float_t _puppimet_ey_JetEnDown=-99;
  Float_t _puppimet_ex_UnclusteredEnUp=-99;
  Float_t _puppimet_ey_UnclusteredEnUp=-99;
  Float_t _puppimet_ex_UnclusteredEnDown=-99;
  Float_t _puppimet_ey_UnclusteredEnDown=-99;
  Float_t _puppimet_ex_JetResUp=-99;
  Float_t _puppimet_ey_JetResUp=-99;
  Float_t _puppimet_ex_JetResDown=-99;
  Float_t _puppimet_ey_JetResDown=-99;

  std::vector<Float_t> mjj;
  std::vector<Float_t> jdeta;
  std::vector<Int_t> njetingap;
  std::vector<Int_t> njetingap20;
  std::vector<Float_t> jdphi;
  std::vector<Float_t> dijetpt;
  std::vector<Float_t> dijetphi;
  std::vector<Float_t> ptvis;

  std::vector<Int_t> nbtag;
  std::vector<Int_t> njets;
  std::vector<Int_t> njetspt20;
  std::vector<Float_t> jpt_1;
  std::vector<Float_t> jeta_1;
  std::vector<Float_t> jphi_1;
  std::vector<Float_t> jcsv_1;
  std::vector<Float_t> jpt_2;
  std::vector<Float_t> jeta_2;
  std::vector<Float_t> jphi_2;
  std::vector<Float_t> jcsv_2;
  std::vector<Float_t> bpt_1;
  std::vector<Float_t> beta_1;
  std::vector<Float_t> bphi_1;
  std::vector<Float_t> bcsv_1;
  std::vector<Float_t> bpt_2;
  std::vector<Float_t> beta_2;
  std::vector<Float_t> bphi_2;
  std::vector<Float_t> bcsv_2;

  std::vector<Float_t> mjjDown;
  std::vector<Float_t> jdetaDown;
  std::vector<Int_t> njetingapDown;
  std::vector<Int_t> njetingap20Down;
  std::vector<Float_t> jdphiDown;
  std::vector<Float_t> dijetptDown;
  std::vector<Float_t> dijetphiDown;
  std::vector<Float_t> ptvisDown;

  std::vector<Int_t> nbtagDown;
  std::vector<Int_t> njetsDown;
  std::vector<Int_t> njetspt20Down;
  std::vector<Float_t> jptDown_1;
  std::vector<Float_t> jetaDown_1;
  std::vector<Float_t> jphiDown_1;
  std::vector<Float_t> jcsvDown_1;
  std::vector<Float_t> jptDown_2;
  std::vector<Float_t> jetaDown_2;
  std::vector<Float_t> jphiDown_2;
  std::vector<Float_t> jcsvDown_2;
  std::vector<Float_t> bptDown_1;
  std::vector<Float_t> betaDown_1;
  std::vector<Float_t> bphiDown_1;
  std::vector<Float_t> bcsvDown_1;
  std::vector<Float_t> bptDown_2;
  std::vector<Float_t> betaDown_2;
  std::vector<Float_t> bphiDown_2;
  std::vector<Float_t> bcsvDown_2;

  std::vector<Float_t> mjjUp;
  std::vector<Float_t> jdetaUp;
  std::vector<Int_t> njetingapUp;
  std::vector<Int_t> njetingap20Up;
  std::vector<Float_t> jdphiUp;
  std::vector<Float_t> dijetptUp;
  std::vector<Float_t> dijetphiUp;
  std::vector<Float_t> ptvisUp;

  std::vector<Int_t> nbtagUp;
  std::vector<Int_t> njetsUp;
  std::vector<Int_t> njetspt20Up;
  std::vector<Float_t> jptUp_1;
  std::vector<Float_t> jetaUp_1;
  std::vector<Float_t> jphiUp_1;
  std::vector<Float_t> jcsvUp_1;
  std::vector<Float_t> jptUp_2;
  std::vector<Float_t> jetaUp_2;
  std::vector<Float_t> jphiUp_2;
  std::vector<Float_t> jcsvUp_2;
  std::vector<Float_t> bptUp_1;
  std::vector<Float_t> betaUp_1;
  std::vector<Float_t> bphiUp_1;
  std::vector<Float_t> bcsvUp_1;
  std::vector<Float_t> bptUp_2;
  std::vector<Float_t> betaUp_2;
  std::vector<Float_t> bphiUp_2;
  std::vector<Float_t> bcsvUp_2;

  Float_t puweight=1;

  Int_t _nup=-99;

  std::vector<Float_t> weight;

  std::vector<Float_t> jpfid_1;
  std::vector<Float_t> jpuid_1;
  std::vector<Float_t> jpfid_2;
  std::vector<Float_t> jpuid_2;
  std::vector<Float_t> bpfid_1;
  std::vector<Float_t> bpuid_1;
  std::vector<Float_t> bpfid_2;
  std::vector<Float_t> bpuid_2;

  Int_t _npv=-99;
  Float_t _npu=-99;
  Float_t _rho=-99;

  std::vector<Float_t> pt_sv;
  std::vector<Float_t> eta_sv;
  std::vector<Float_t> phi_sv;
  std::vector<Float_t> met_sv;

  std::vector<Float_t> byDeepTau2017v2p1VSjetraw_1;
  std::vector<Float_t> byVVVLooseDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byVVLooseDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byTightDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byVTightDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> byVVTightDeepTau2017v2p1VSjet_1;

  std::vector<Float_t> byDeepTau2017v2p1VSeleraw_1;
  std::vector<Float_t> byVVVLooseDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byVVLooseDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byTightDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byVTightDeepTau2017v2p1VSe_1;
  std::vector<Float_t> byVVTightDeepTau2017v2p1VSe_1;

  std::vector<Float_t> byDeepTau2017v2p1VSmuraw_1;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSmu_1;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSmu_1;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSmu_1;
  std::vector<Float_t> byTightDeepTau2017v2p1VSmu_1;

  std::vector<Float_t> byDeepTau2017v2p1VSjetraw_2;
  std::vector<Float_t> byVVVLooseDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byVVLooseDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byTightDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byVTightDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> byVVTightDeepTau2017v2p1VSjet_2;

  std::vector<Float_t> byDeepTau2017v2p1VSeleraw_2;
  std::vector<Float_t> byVVVLooseDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byVVLooseDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byTightDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byVTightDeepTau2017v2p1VSe_2;
  std::vector<Float_t> byVVTightDeepTau2017v2p1VSe_2;

  std::vector<Float_t> byDeepTau2017v2p1VSmuraw_2;
  std::vector<Float_t> byVLooseDeepTau2017v2p1VSmu_2;
  std::vector<Float_t> byLooseDeepTau2017v2p1VSmu_2;
  std::vector<Float_t> byMediumDeepTau2017v2p1VSmu_2;
  std::vector<Float_t> byTightDeepTau2017v2p1VSmu_2;

  float pvx;
  float pvy;
  float pvz;
  std::vector<std::vector<double>> pvcov;
  std::vector<float> pcax;
  std::vector<float> pcay;
  std::vector<float> pcaz;

  std::vector<Float_t> dm_1;
  std::vector<Float_t> dmMVA_1;

  std::vector<Float_t> dm_2;
  std::vector<Float_t> dmMVA_2;

  std::vector<Float_t> _chpt_1;
  std::vector<Float_t> _chphi_1;
  std::vector<Float_t> _cheta_1;
  std::vector<Float_t> _chm_1;

  std::vector<Float_t> _chpt_2;
  std::vector<Float_t> _chphi_2;
  std::vector<Float_t> _cheta_2;
  std::vector<Float_t> _chm_2;

  std::vector<Float_t> _npt_1;
  std::vector<Float_t> _nphi_1;
  std::vector<Float_t> _neta_1;
  std::vector<Float_t> _nm_1;

  std::vector<Float_t> _npt_2;
  std::vector<Float_t> _nphi_2;
  std::vector<Float_t> _neta_2;
  std::vector<Float_t> _nm_2;

  std::vector<Float_t> svx_1;
  std::vector<Float_t> svy_1;
  std::vector<Float_t> svz_1;

  std::vector<Float_t> svx_2;
  std::vector<Float_t> svy_2;
  std::vector<Float_t> svz_2;

  PileUp * PUofficial = new PileUp();

  std::vector<math::XYZTLorentzVector> LeptonP4;
  std::vector<Int_t> tau1IndexVect;
  std::vector<Int_t> tau2IndexVect;

  // Output trees
  TTree *Nominal;
  TTree *TESUp;
  TTree *TESDown;
  TTree *MESUp;
  TTree *MESDown;
  TTree *JERUp;
  TTree *JERDown;
  TTree *METResoUp;
  TTree *METResoDown;
  TTree *METScaleUp;
  TTree *METScaleDown;
  TTree *FlavorQCDUp;
  TTree *FlavorQCDDown;
  TTree *RelativeBalUp;
  TTree *RelativeBalDown;
  TTree *HFUp;
  TTree *HFDown;
  TTree *BBEC1Up;
  TTree *BBEC1Down;
  TTree *EC2Up;
  TTree *EC2Down;
  TTree *AbsoluteUp;
  TTree *AbsoluteDown;
  TTree *BBEC1_YEARUp;
  TTree *BBEC1_YEARDown;
  TTree *EC2_YEARUp;
  TTree *EC2_YEARDown;
  TTree *Absolute_YEARUp;
  TTree *Absolute_YEARDown;
  TTree *HF_YEARUp;
  TTree *HF_YEARDown;
  TTree *RelativeSample_YEARUp;
  TTree *RelativeSample_YEARDown;

  std::map<std::string, TTree*> SystematicsMap {
    {"Nominal",  Nominal},
    {"TESUp",  TESUp},
    {"TESDown",  TESDown}, 
    {"MESUp",  MESUp}, 
    {"MESDown",  MESDown}, 
    {"JERUp",  JERUp}, 
    {"JERDown",  JERDown}, 
    {"METResoUp",  METResoUp}, 
    {"METResoDown",  METResoDown}, 
    {"METScaleUp",  METScaleUp}, 
    {"METScaleDown",  METScaleDown}, 
    {"FlavorQCDUp",  FlavorQCDUp}, 
    {"FlavorQCDDown",  FlavorQCDDown}, 
    {"RelativeBalUp",  RelativeBalUp}, 
    {"RelativeBalDown",  RelativeBalDown}, 
    {"HFUp",  HFUp}, 
    {"HFDown",  HFDown}, 
    {"BBEC1Up",  BBEC1Up}, 
    {"BBEC1Down",  BBEC1Down}, 
    {"EC2Up",  EC2Up}, 
    {"EC2Down",  EC2Down}, 
    {"AbsoluteUp",  AbsoluteUp}, 
    {"AbsoluteDown",  AbsoluteDown}, 
    {"BBEC1_YEARUp",  BBEC1_YEARUp}, 
    {"BBEC1_YEARDown",  BBEC1_YEARDown}, 
    {"EC2_YEARUp",  EC2_YEARUp}, 
    {"EC2_YEARDown",  EC2_YEARDown}, 
    {"Absolute_YEARUp",  Absolute_YEARUp}, 
    {"Absolute_YEARDown",  Absolute_YEARDown}, 
    {"HF_YEARUp",  HF_YEARUp}, 
    {"HF_YEARDown",  HF_YEARDown}, 
    {"RelativeSample_YEARUp",  RelativeSample_YEARUp}, 
    {"RelativeSample_YEARDown",  RelativeSample_YEARDown}
  };

};

const int HTauTauNtuplizer::ntauIds; // definition of static member

// ----Constructor and Destructor -----
HTauTauNtuplizer::HTauTauNtuplizer(const edm::ParameterSet& pset) : //reweight(),
  triggerObjects_      (consumes<pat::TriggerObjectStandAloneCollection> (pset.getParameter<edm::InputTag>("triggerSet"))),
  //l1ExtraIsoTau_       (consumes<vector<l1extra::L1JetParticle>>         (pset.getParameter<edm::InputTag>("l1extraIsoTau"))) ,
  triggerBits_         (consumes<edm::TriggerResults>                    (pset.getParameter<edm::InputTag>("triggerResultsLabel"))),
  metFilterBits_       (consumes<edm::TriggerResults>                    (pset.getParameter<edm::InputTag>("metFilters"))),
  theVtxTag            (consumes<vector<Vertex>>                         (pset.getParameter<edm::InputTag>("vtxCollection"))),
  theSecVtxTag         (consumes<edm::View<reco::VertexCompositePtrCandidate>> (pset.getParameter<edm::InputTag>("secVtxCollection"))), //FRA
  theRhoTag            (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoCollection"))),
  theRhoMiniRelIsoTag  (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoMiniRelIsoCollection"))),
  theRhoForJERTag      (consumes<double>                                 (pset.getParameter<edm::InputTag>("rhoForJER"))), //FRA
  thePUTag             (consumes<vector<PileupSummaryInfo>>              (pset.getParameter<edm::InputTag>("puCollection"))),
  thePFCandTag         (consumes<edm::View<pat::PackedCandidate>>        (pset.getParameter<edm::InputTag>("PFCandCollection"))),
  theCandTag           (consumes<edm::View<pat::CompositeCandidate>>     (pset.getParameter<edm::InputTag>("candCollection"))),
  theJetTag            (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("jetCollection"))),
  theSmearedJetTag     (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("SmearedJets"))),
  theSmearedJetDownTag (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("SmearedJetsDown"))),
  theSmearedJetUpTag   (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("SmearedJetsUp"))),
  theFatJetTag         (consumes<edm::View<pat::Jet>>                    (pset.getParameter<edm::InputTag>("ak8jetCollection"))),
  theQGTaggerTag       (consumes<edm::ValueMap<float>>                   (pset.getParameter<edm::InputTag>("QGTagger"))),
  theLepTag            (consumes<edm::View<reco::Candidate>>             (pset.getParameter<edm::InputTag>("lepCollection"))),
  theLHETag            (consumes<LHEEventProduct>                        (pset.getParameter<edm::InputTag>("lheCollection"))),
  theGenTag            (consumes<GenEventInfoProduct>                    (pset.getParameter<edm::InputTag>("genCollection"))),
//  theMetTag            (consumes<pat::METCollection>                     (pset.getParameter<edm::InputTag>("metCollection"))),
// theMetERTag          (consumes<pat::METCollection>                     (pset.getParameter<edm::InputTag>("metERCollection"))),
  thePUPPIMetTag       (consumes<pat::METCollection>                     (pset.getParameter<edm::InputTag>("PUPPImetCollection"))),
// thePFMETCovTag       (consumes<math::Error<2>::type>                   (pset.getParameter<edm::InputTag>("srcPFMETCov"))),
//thePFMETSignifTag    (consumes<double>                                 (pset.getParameter<edm::InputTag>("srcPFMETSignificance"))),
  thePUPPIMETCovTag    (consumes<math::Error<2>::type>                   (pset.getParameter<edm::InputTag>("srcPUPPIMETCov"))),
  thePUPPIMETSignifTag (consumes<double>                                 (pset.getParameter<edm::InputTag>("srcPUPPIMETSignificance"))),
  theGenericTag        (consumes<edm::View<pat::GenericParticle>>        (pset.getParameter<edm::InputTag>("genericCollection"))),
  theGenJetTag         (consumes<edm::View<reco::GenJet>>                (pset.getParameter<edm::InputTag>("genjetCollection"))),
  theTotTag            (consumes<edm::MergeableCounter, edm::InLumi>     (pset.getParameter<edm::InputTag>("totCollection"))),
  thePassTag           (consumes<edm::MergeableCounter, edm::InLumi>     (pset.getParameter<edm::InputTag>("passCollection"))),
  theLHEPTag           (consumes<LHEEventProduct>                        (pset.getParameter<edm::InputTag>("lhepCollection"))),
  beamSpotTag          (consumes<reco::BeamSpot>                         (pset.getParameter<edm::InputTag>("beamSpot"))),
  theL1TauTag          (consumes<BXVector<l1t::Tau>>                     (pset.getParameter<edm::InputTag>("stage2TauCollection"))),
  theL1JetTag          (consumes<BXVector<l1t::Jet>>                     (pset.getParameter<edm::InputTag>("stage2JetCollection"))),
//								    theNBadMuTag         (consumes<int>                                    (pset.getParameter<edm::InputTag>("nBadMu"))),
  genLumiHeaderTag     (consumes<GenLumiInfoHeader, edm::InLumi>         (pset.getParameter<edm::InputTag>("genLumiHeaderTag"))),
  ecalBadCalibFilterUpdate_token  (consumes< bool >                      (pset.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))),
  prefweight_token     (consumes< double >                               (pset.getParameter<edm::InputTag>("L1prefireProb"))),
  prefweightup_token   (consumes< double >                               (pset.getParameter<edm::InputTag>("L1prefireProbUp"))),
  prefweightdown_token (consumes< double >                               (pset.getParameter<edm::InputTag>("L1prefireProbDown"))),
  ThePrunedGenTag_     (consumes<edm::View<reco::GenParticle> >          (pset.getParameter<edm::InputTag>("PrunedGenCollection"))),
  theNewTauTag         (consumes<pat::TauCollection>                     (pset.getParameter<edm::InputTag>("TausNewIDCollection"))),
  theNewTauTagCandidate(consumes<edm::View<reco::Candidate> >                   (pset.getParameter<edm::InputTag>("TausNewIDCollection"))),
  RefitVtxBSTag        (consumes<RefitVertexCollection>  (pset.getParameter<edm::InputTag>("RefitVtxBSCollection"))),
  RefitVtxNoBSTag      (consumes<RefitVertexCollection>  (pset.getParameter<edm::InputTag>("RefitVtxNoBSCollection"))),
  theTauSpinnerWTEventag   (consumes<double> (pset.getParameter<edm::InputTag>("tauSpinnerWTEven"))),
  theTauSpinnerWTOddtag   (consumes<double> (pset.getParameter<edm::InputTag>("tauSpinnerWTOdd"))),
  theTauSpinnerWTMMtag   (consumes<double> (pset.getParameter<edm::InputTag>("tauSpinnerWTMM")))

{
  theFileName = pset.getUntrackedParameter<string>("fileName");
  theFSR = pset.getParameter<bool>("applyFSR");
  theisMC = pset.getParameter<bool>("IsMC");
  IsEmbed=pset.getParameter<bool>("IsEmbed");
  do_MCSummary_=pset.getParameter<bool>("do_MCSummary");
  do_MCComplete_=pset.getParameter<bool>("do_MCComplete");
  doCPVariables = pset.getParameter<bool>("doCPVariables");
  computeQGVar = pset.getParameter<bool>("computeQGVar");
  theJECName = pset.getUntrackedParameter<string>("JECset");
  theYear = pset.getParameter<int>("year");
  // theUseNoHFPFMet = pset.getParameter<bool>("useNOHFMet");
  //writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  //sampleName = pset.getParameter<string>("sampleName");
  Nevt_Gen=0;
  Nevt_PassTrigger = 0;
  Npairs=0;
  SelectedPairs=0;

  ///// TRIGGER
  processName= pset.getParameter<edm::InputTag>("triggerResultsLabel");
  std::vector<edm::ParameterSet> HLTList = pset.getParameter <std::vector<edm::ParameterSet> > ("triggerList");

  DataMCType DMT;
  dataMCtype = DMT.GetType();
  mySysHelper = new SysHelper(theYear,dataMCtype);
  myTriggerHelper = new triggerhelper();// (HLTList);

  if (HLTList.size() > 64) cout << endl << "** HTauTauNtuplizer : Warning : trigger list size exceeds 64, not enough bits in Long64_t type to store all" << endl;
  if (HLTList.size() > 64) cout << "** HTauTauNtuplizer : Warning : trigger list size: " << HLTList.size() << endl << endl;

  for (std::vector<edm::ParameterSet>::const_iterator iPSet = HLTList.begin();iPSet != HLTList.end(); ++iPSet) {
    const std::string& hlt = iPSet->getParameter<std::string>("HLT");
    const std::vector<std::string>& path1 = iPSet->getParameter<std::vector<std::string>>("path1");
    const std::vector<std::string>& path2 = iPSet->getParameter<std::vector<std::string>>("path2");
    const std::vector<std::string>& path3 = iPSet->getParameter<std::vector<std::string>>("path3"); //FRA
    const std::vector<std::string>& path4 = iPSet->getParameter<std::vector<std::string>>("path4"); //FRA
    const int& leg1 = iPSet->getParameter<int>("leg1");
    const int& leg2 = iPSet->getParameter<int>("leg2");
    myTriggerHelper->addTriggerMap(hlt,path1,path2,path3,path4,leg1,leg2); //FRA
    //myTriggerHelper->addTriggerMap(hlt,path1,path2,path3,path4,leg1,leg2, pt1, pt2); //FRA
    // Build the map
    //myTriggerHelper->addTriggerMap(hlt,path1,path2,path3,path4,leg1,leg2, pt1, pt2); //FRA
  }

  Initialize();
  _susyModel = ""; // not to initialize at every event, as updated ni new lumi block
}

HTauTauNtuplizer::~HTauTauNtuplizer(){
  delete myTriggerHelper;
  delete mySysHelper;
}
//

void HTauTauNtuplizer::Initialize(){

  //_mothers.clear();
  _mothers_px.clear();
  _mothers_py.clear();
  _mothers_pz.clear();
  _mothers_e.clear();
  _mothers_trgSeparateMatch.clear();

  //_daughters.clear();
  if(DEBUG){
    _trigger_name.clear();
    _trigger_accept.clear();
  }
  _daughters_px.clear();
  _daughters_py.clear();
  _daughters_pz.clear();
  _daughters_e.clear();

  _daughters_charged_px.clear();
  _daughters_charged_py.clear();
  _daughters_charged_pz.clear();
  _daughters_charged_e.clear();

  _daughters_neutral_px.clear();
  _daughters_neutral_py.clear();
  _daughters_neutral_pz.clear();
  _daughters_neutral_e.clear();

  _daughters_vx.clear();
  _daughters_vy.clear();
  _daughters_vz.clear();

  _daughters_hasTES.clear();
  _daughters_px_TauUp.clear();
  _daughters_py_TauUp.clear();
  _daughters_pz_TauUp.clear();
  _daughters_e_TauUp.clear();
  _daughters_px_TauDown.clear();
  _daughters_py_TauDown.clear();
  _daughters_pz_TauDown.clear();
  _daughters_e_TauDown.clear();
  _daughters_hasEES.clear();
  _daughters_px_EleUp.clear();
  _daughters_py_EleUp.clear();
  _daughters_pz_EleUp.clear();
  _daughters_e_EleUp.clear();
  _daughters_px_EleDown.clear();
  _daughters_py_EleDown.clear();
  _daughters_pz_EleDown.clear();
  _daughters_e_EleDown.clear();
  _daughters_charge.clear();
  _daughters_isTauMatched.clear();
  _daughters_genindex.clear();
  _daughters_IetaIeta.clear();
  _daughters_full5x5_IetaIeta.clear();
  _daughters_hOverE.clear();
  _daughters_deltaEtaSuperClusterTrackAtVtx.clear();
  _daughters_deltaPhiSuperClusterTrackAtVtx.clear();
  _daughters_IoEmIoP.clear();
  _daughters_IoEmIoP_ttH.clear();
  //  _daughters_SCeta.clear();
  _daughters_depositR03_tracker.clear();
  _daughters_depositR03_ecal.clear();
  _daughters_depositR03_hcal.clear();
  _daughters_decayModeFindingOldDMs.clear();
  _daughters_decayModeFindingNewDMs.clear();
  _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();

  //_daughters_byIsolationMVArun2v1DBoldDMwLTraw.clear();
  _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017.clear();      //FRA
  //_daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017.clear();      //FRA
  //_daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017.clear(); //FRA
  _daughters_byDeepTau2017v2p1VSjetraw.clear();
  _daughters_byDeepTau2017v2p1VSeraw.clear();
  _daughters_byDeepTau2017v2p1VSmuraw.clear();
  _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear(); //FRA

  _daughters_byIsolationMVArun2v1DBoldDMwLTraw.clear();
  _daughters_byDeepTau2017v2p1VSjetraw.clear();
  _daughters_byDeepTau2017v2p1VSeraw.clear();
  _daughters_byDeepTau2017v2p1VSmuraw.clear();

  //_MVADM2016v1.clear();
  _MVADM2017v1.clear();
  //_MVADM2017v1InTauFiller.clear();

  _daughters_footprintCorrection.clear();
  _daughters_neutralIsoPtSumWeight.clear();
  _daughters_photonPtSumOutsideSignalCone.clear();

  _daughters_chargedIsoPtSum.clear();
  _daughters_neutralIsoPtSum.clear();
  _daughters_puCorrPtSum.clear();
  _daughters_numChargedParticlesSignalCone.clear();
  _daughters_numNeutralHadronsSignalCone.clear();
  _daughters_numPhotonsSignalCone.clear();
  _daughters_numParticlesSignalCone.clear();
  _daughters_numChargedParticlesIsoCone.clear();
  _daughters_numNeutralHadronsIsoCone.clear();
  _daughters_numPhotonsIsoCone.clear();
  _daughters_numParticlesIsoCone.clear();
  _daughters_leadChargedParticlePt.clear();
  _daughters_trackRefPt.clear();
  //  _daughters_LFtrigger.clear();
  //_daughters_L3trigger.clear();
  _daughters_trgMatched.clear();
  _daughters_FilterFired.clear();
  _daughters_isGoodTriggerType.clear();
  _daughters_L3FilterFired.clear();
  _daughters_L3FilterFiredLast.clear();
  _PFTauSVPos.clear();
  _PFTauSVCov.clear();

  _PFTauPionsP4.clear();
  _PFTauRefitPionsP4.clear();
  _PFTauPionsCharge.clear();
  _PFTauRefitPionsCharge.clear();
  _PFTauSVChi2NDofMatchingQuality.clear();


  //  _daughters_jetNDauChargedMVASel.clear();
  _daughters_miniRelIsoCharged.clear();
  _daughters_miniRelIsoNeutral.clear();
  _daughters_jetPtRel.clear();
  _daughters_jetPtRatio.clear();
  _daughters_jetBTagCSV.clear();
  _daughters_jetBTagDeepCSV.clear();
  _daughters_jetBTagDeepFlavor.clear();

  //  _daughters_iseleBDT.clear();
  _daughters_iseleWPLoose.clear();
  _daughters_iseleWP80.clear();
  _daughters_iseleWP90.clear();
  _daughters_iseleNoIsoWPLoose.clear();
  _daughters_iseleNoIsoWP80.clear();
  _daughters_iseleNoIsoWP90.clear();
  _daughters_iseleCutBased.clear();

  //  _daughters_eleMVAntNoIso.clear();
  _daughters_eleMVAnt.clear();
  _daughters_eleMVA_HZZ.clear();
  _daughters_passConversionVeto.clear();
  _daughters_eleMissingHits.clear();
  //_daughters_eleMissingLostHits.clear();
  _daughters_iseleChargeConsistent.clear();
  //_daughters_iseleCUT.clear();
  //_daughter2.clear();
  _softLeptons.clear();
  //_genDaughters.clear();

  _daughters_pca_x.clear();
  _daughters_pca_y.clear();
  _daughters_pca_z.clear();

  _daughters_pcaRefitPV_x.clear();
  _daughters_pcaRefitPV_y.clear();
  _daughters_pcaRefitPV_z.clear();

  _daughters_pcaRefitPVBS_x.clear();
  _daughters_pcaRefitPVBS_y.clear();
  _daughters_pcaRefitPVBS_z.clear();

  _daughters_pcaGenPV_x.clear();
  _daughters_pcaGenPV_y.clear();
  _daughters_pcaGenPV_z.clear();

  _genpart_px.clear();
  _genpart_py.clear();
  _genpart_pz.clear();
  _genpart_e.clear();

  _genpart_pca_x.clear();
  _genpart_pca_y.clear();
  _genpart_pca_z.clear();

  _genpart_pdg.clear();
  _genpart_status.clear();
  //_genpart_mothInd.clear();
  _genpart_HMothInd.clear();
  _genpart_MSSMHMothInd.clear();
  _genpart_TopMothInd.clear();
  _genpart_TauMothInd.clear();
  _genpart_ZMothInd.clear();
  _genpart_WMothInd.clear();
  _genpart_bMothInd.clear();
  _genpart_HZDecayMode.clear();
  _genpart_TopDecayMode.clear();
  _genpart_WDecayMode.clear();
  _genpart_TauGenDecayMode.clear();
  _genpart_TauGenDetailedDecayMode.clear();
  _genpart_flags.clear();

  _L1_tauEt.clear();
  _L1_tauEta.clear();
  _L1_tauPhi.clear();
  _L1_tauIso.clear();

  _L1_jetEt.clear();
  _L1_jetEta.clear();
  _L1_jetPhi.clear();

  _genjet_px.clear();
  _genjet_py.clear();
  _genjet_pz.clear();
  _genjet_e.clear();
  _genjet_partonFlavour.clear();
  _genjet_hadronFlavour.clear();

  _indexDau1.clear();
  _indexDau2.clear();
  _pdgdau.clear();
  _isOSCand.clear();
  _daughters_HLTpt.clear();
  _daughters_highestEt_L1IsoTauMatched.clear();
  _metx.clear();
  _mety.clear();
  _mTDau1.clear();
  _mTDau2.clear();
  _particleType.clear();
  _daughters_typeOfMuon.clear();
  _daughters_muonID.clear();
  _dxy.clear();
  _dz.clear();
  _dxy_innerTrack.clear();
  _dz_innerTrack.clear();
  _daughters_rel_error_trackpt.clear();
  _SIP.clear();
  _decayType.clear();
  _genmatch.clear();
  _daughters_tauID.clear();
  _combreliso.clear();
  _combreliso03.clear();
  _indexevents=0;
  _runNumber=0;
  _lumi=0;
  _year=0;
  _passecalBadCalibFilterUpdate=false;
  _prefiringweight = 1.;
  _prefiringweightup = 1.;
  _prefiringweightdown = 1.;
  _triggerbit=0;
  _metfilterbit=0;
  _MC_weight=0.;
  _nominal_wt=0.;
  _TheoreticalPSUnc.clear();
  _DataMC_Type=0;
  DataMC_Type_idx=0;
  Event_isRealData=0;
  //npv=0;
  _lheHt=0;
  _lheNOutPartons=0;
  _lheNOutB=0;
  _lheNOutC=0;
  //npu=0.;
  _PUNumInteractions=0;
  _MC_weight_scale_muF0p5=0.;
  _MC_weight_scale_muF2=0.;
  _MC_weight_scale_muR0p5=0.;
  _MC_weight_scale_muR2=0.;

  if (do_MCComplete_ || IsEmbed) {
    MC_p4.clear();
    MC_pdgid.clear();
    MC_charge.clear();
    MC_midx.clear();
    MC_status.clear();
    MC_childpdgid.clear();
    MC_childidx.clear();
  }

  if (do_MCSummary_ || IsEmbed) {
    MCSignalParticle_p4.clear();
    MCSignalParticle_pdgid.clear();
    MCSignalParticle_charge.clear();
    MCSignalParticle_Poca.clear();
    MCSignalParticle_Tauidx.clear();
    MCTauandProd_p4.clear();
    MCTauProdVisible_P4.clear();
    MCTauandProd_Vertex.clear();
    MCTauandProd_pdgid.clear();
    MCTauandProd_midx.clear();
    MCTauandProd_charge.clear();
    MCTau_JAK.clear();
    MCTau_DecayBitMask.clear();
    MC_childpdgid.clear();
  }

  //_jets.clear();
  _jets_VBFleadFilterMatch.clear();    //FRA
  _jets_VBFsubleadFilterMatch.clear(); //FRA
  _jets_px.clear();
  _jets_py.clear();
  _jets_pz.clear();
  _jets_e.clear();
  //_jets_rawPt.clear();
  _jets_area.clear();
  _jets_mT.clear();
  _jets_vtxPt.clear();
  _jets_vtxMass.clear();
  _jets_vtx3dL.clear();
  _jets_vtxNtrk.clear();
  _jets_vtx3deL.clear();
  _jets_leadTrackPt.clear();
  _jets_leptonPtRel.clear();
  _jets_leptonPt.clear();
  _jets_leptonDeltaR.clear();
  _jets_chEmEF.clear();
  _jets_chHEF.clear();
  _jets_nEmEF.clear();
  _jets_nHEF.clear();
  _jets_MUF.clear();
  _jets_neMult.clear();
  _jets_chMult.clear();
  _jets_Flavour.clear();
  _jets_HadronFlavour.clear();
  _jets_genjetIndex.clear();
  _jets_jecUnc.clear();
  _pileupjetidMVA.clear();

  // JEC uncertainties sources Regrouped
  _jets_jetUncRegrouped_FlavorQCD_up.clear();  // up variations
  _jets_jetUncRegrouped_RelativeBal_up.clear();
  _jets_jetUncRegrouped_HF_up.clear();
  _jets_jetUncRegrouped_BBEC1_up.clear();
  _jets_jetUncRegrouped_EC2_up.clear();
  _jets_jetUncRegrouped_Absolute_up.clear();
  _jets_jetUncRegrouped_BBEC1_YEAR_up.clear();
  _jets_jetUncRegrouped_EC2_YEAR_up.clear();
  _jets_jetUncRegrouped_Absolute_YEAR_up.clear();
  _jets_jetUncRegrouped_HF_YEAR_up.clear();
  _jets_jetUncRegrouped_RelativeSample_YEAR_up.clear();
  _jets_jetUncRegrouped_Total_up.clear();
  _jets_jetUncRegrouped_FlavorQCD_dw.clear(); // down variations
  _jets_jetUncRegrouped_RelativeBal_dw.clear();
  _jets_jetUncRegrouped_HF_dw.clear();
  _jets_jetUncRegrouped_BBEC1_dw.clear();
  _jets_jetUncRegrouped_EC2_dw.clear();
  _jets_jetUncRegrouped_Absolute_dw.clear();
  _jets_jetUncRegrouped_BBEC1_YEAR_dw.clear();
  _jets_jetUncRegrouped_EC2_YEAR_dw.clear();
  _jets_jetUncRegrouped_Absolute_YEAR_dw.clear();
  _jets_jetUncRegrouped_HF_YEAR_dw.clear();
  _jets_jetUncRegrouped_RelativeSample_YEAR_dw.clear();
  _jets_jetUncRegrouped_Total_dw.clear();
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncValRegrouped_up.begin(); it!=_SourceUncValRegrouped_up.end(); ++it)
    {
      it->second.clear();
    }
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncValRegrouped_dw.begin(); it!=_SourceUncValRegrouped_dw.end(); ++it)
    {
      it->second.clear();
    }

  //  _TheoreticalScaleUnc.clear();

  for (int Idx = 0; Idx < 18; Idx++) {
    _TheoreticalScaleUncTab[Idx] = -999;
  }

  _jetsDown_px.clear();
  _jetsDown_py.clear();
  _jetsDown_pz.clear();
  _jetsDown_e.clear();
  _jetsDown_area.clear();
  _jetsDown_mT.clear();
  _jetsDown_leadTrackPt.clear();
  _jetsDown_leptonPtRel.clear();
  _jetsDown_leptonPt.clear();
  _jetsDown_leptonDeltaR.clear();
  _jetsDown_chEmEF.clear();
  _jetsDown_chHEF.clear();
  _jetsDown_nEmEF.clear();
  _jetsDown_nHEF.clear();
  _jetsDown_MUF.clear();
  _jetsDown_neMult.clear();
  _jetsDown_chMult.clear();
  _jetsDown_Flavour.clear();
  _jetsDown_HadronFlavour.clear();
  _jetsDown_genjetIndex.clear();
  _pileupjetidMVADown.clear();

  _jetsUp_px.clear();
  _jetsUp_py.clear();
  _jetsUp_pz.clear();
  _jetsUp_e.clear();
  _jetsUp_area.clear();
  _jetsUp_mT.clear();
  _jetsUp_leadTrackPt.clear();
  _jetsUp_leptonPtRel.clear();
  _jetsUp_leptonPt.clear();
  _jetsUp_leptonDeltaR.clear();
  _jetsUp_chEmEF.clear();
  _jetsUp_chHEF.clear();
  _jetsUp_nEmEF.clear();
  _jetsUp_nHEF.clear();
  _jetsUp_MUF.clear();
  _jetsUp_neMult.clear();
  _jetsUp_chMult.clear();
  _jetsUp_Flavour.clear();
  _jetsUp_HadronFlavour.clear();
  _jetsUp_genjetIndex.clear();
  _pileupjetidMVAUp.clear();

  // JEC uncertainty sources
  _jets_jetUnc_AbsoluteFlavMap_up.clear(); //up variations
  _jets_jetUnc_AbsoluteMPFBias_up.clear();
  _jets_jetUnc_AbsoluteSample_up.clear();
  _jets_jetUnc_AbsoluteScale_up.clear();
  _jets_jetUnc_AbsoluteStat_up.clear();
  _jets_jetUnc_FlavorQCD_up.clear();
  _jets_jetUnc_Fragmentation_up.clear();
  _jets_jetUnc_PileUpDataMC_up.clear();
  _jets_jetUnc_PileUpPtBB_up.clear();
  _jets_jetUnc_PileUpPtEC1_up.clear();
  _jets_jetUnc_PileUpPtEC2_up.clear();
  _jets_jetUnc_PileUpPtHF_up.clear();
  _jets_jetUnc_PileUpPtRef_up.clear();
  _jets_jetUnc_RelativeBal_up.clear();
  _jets_jetUnc_RelativeFSR_up.clear();
  _jets_jetUnc_RelativeJEREC1_up.clear();
  _jets_jetUnc_RelativeJEREC2_up.clear();
  _jets_jetUnc_RelativeJERHF_up.clear();
  _jets_jetUnc_RelativePtBB_up.clear();
  _jets_jetUnc_RelativePtEC1_up.clear();
  _jets_jetUnc_RelativePtEC2_up.clear();
  _jets_jetUnc_RelativePtHF_up.clear();
  _jets_jetUnc_RelativeSample_up.clear();
  _jets_jetUnc_RelativeStatEC_up.clear();
  _jets_jetUnc_RelativeStatFSR_up.clear();
  _jets_jetUnc_RelativeStatHF_up.clear();
  _jets_jetUnc_SinglePionECAL_up.clear();
  _jets_jetUnc_SinglePionHCAL_up.clear();
  _jets_jetUnc_TimePtEta_up.clear();
  _jets_jetUnc_AbsoluteFlavMap_dw.clear(); // down variations
  _jets_jetUnc_AbsoluteMPFBias_dw.clear();
  _jets_jetUnc_AbsoluteSample_dw.clear();
  _jets_jetUnc_AbsoluteScale_dw.clear();
  _jets_jetUnc_AbsoluteStat_dw.clear();
  _jets_jetUnc_FlavorQCD_dw.clear();
  _jets_jetUnc_Fragmentation_dw.clear();
  _jets_jetUnc_PileUpDataMC_dw.clear();
  _jets_jetUnc_PileUpPtBB_dw.clear();
  _jets_jetUnc_PileUpPtEC1_dw.clear();
  _jets_jetUnc_PileUpPtEC2_dw.clear();
  _jets_jetUnc_PileUpPtHF_dw.clear();
  _jets_jetUnc_PileUpPtRef_dw.clear();
  _jets_jetUnc_RelativeBal_dw.clear();
  _jets_jetUnc_RelativeFSR_dw.clear();
  _jets_jetUnc_RelativeJEREC1_dw.clear();
  _jets_jetUnc_RelativeJEREC2_dw.clear();
  _jets_jetUnc_RelativeJERHF_dw.clear();
  _jets_jetUnc_RelativePtBB_dw.clear();
  _jets_jetUnc_RelativePtEC1_dw.clear();
  _jets_jetUnc_RelativePtEC2_dw.clear();
  _jets_jetUnc_RelativePtHF_dw.clear();
  _jets_jetUnc_RelativeSample_dw.clear();
  _jets_jetUnc_RelativeStatEC_dw.clear();
  _jets_jetUnc_RelativeStatFSR_dw.clear();
  _jets_jetUnc_RelativeStatHF_dw.clear();
  _jets_jetUnc_SinglePionECAL_dw.clear();
  _jets_jetUnc_SinglePionHCAL_dw.clear();
  _jets_jetUnc_TimePtEta_dw.clear();
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncVal_up.begin(); it!=_SourceUncVal_up.end(); ++it)
    {
      it->second.clear();
    }
  for (std::map<std::string, std::vector<Float_t> >::iterator it=_SourceUncVal_dw.begin(); it!=_SourceUncVal_dw.end(); ++it)
    {
      it->second.clear();
    }

  _jets_QGdiscr.clear();
  //_numberOfJets=0;
  _bdiscr.clear();
  _bdiscr2.clear();
  _bdiscr3.clear();
  _bdiscr4.clear();
  _bdiscr5.clear();
  _bdiscr6.clear();
  _bdiscr7.clear();
  _bdiscr8.clear();
  _bdiscr9.clear();
  _bdiscr10.clear();
  _bdiscr11.clear();
  _bdiscr12.clear();
  _bdiscr13.clear();
  _bdiscr14.clear();

  _bdiscrDown.clear();
  _bdiscr2Down.clear();
  _bdiscr3Down.clear();
  _bdiscr4Down.clear();
  _bdiscr5Down.clear();
  _bdiscr6Down.clear();
  _bdiscr7Down.clear();
  _bdiscr8Down.clear();
  _bdiscr9Down.clear();
  _bdiscr10Down.clear();
  _bdiscr11Down.clear();
  _bdiscr12Down.clear();
  _bdiscr13Down.clear();
  _bdiscr14Down.clear();

  _bdiscrUp.clear();
  _bdiscr2Up.clear();
  _bdiscr3Up.clear();
  _bdiscr4Up.clear();
  _bdiscr5Up.clear();
  _bdiscr6Up.clear();
  _bdiscr7Up.clear();
  _bdiscr8Up.clear();
  _bdiscr9Up.clear();
  _bdiscr10Up.clear();
  _bdiscr11Up.clear();
  _bdiscr12Up.clear();
  _bdiscr13Up.clear();
  _bdiscr14Up.clear();

  //_jetID.clear();
  _looseJetID .clear();
  _tightJetID .clear();
  _tightLepVetoJetID .clear();
  _jetrawf.clear();
  _jets_JER.clear(); // Jet Energy Resolution

  _looseJetIDDown.clear();
  _tightJetIDDown.clear();
  _tightLepVetoJetIDDown.clear();
  _jetrawfDown.clear();

  _looseJetIDUp.clear();
  _tightJetIDUp.clear();
  _tightLepVetoJetIDUp.clear();
  _jetrawfUp.clear();

  _ak8jets_px.clear();
  _ak8jets_py.clear();
  _ak8jets_pz.clear();
  _ak8jets_e.clear();
  _ak8jets_SoftDropMass.clear();
  _ak8jets_PrunedMass.clear();
  _ak8jets_TrimmedMass.clear();
  _ak8jets_FilteredMass.clear();
  _ak8jets_tau1.clear();
  _ak8jets_tau2.clear();
  _ak8jets_tau3.clear();
  _ak8jets_tau4.clear();
  _ak8jets_CSV.clear();
  _ak8jets_deepCSV_probb.clear();
  _ak8jets_deepCSV_probbb.clear();
  _ak8jets_deepFlavor_probb.clear();
  _ak8jets_deepFlavor_probbb.clear();
  _ak8jets_deepFlavor_problepb.clear();
  _ak8jets_nsubjets.clear();

  _subjets_px.clear();
  _subjets_py.clear();
  _subjets_pz.clear();
  _subjets_e.clear();
  _subjets_CSV.clear();
  _subjets_deepCSV_probb.clear();
  _subjets_deepCSV_probbb.clear();
  _subjets_deepFlavor_probb.clear();
  _subjets_deepFlavor_probbb.clear();
  _subjets_deepFlavor_problepb.clear();
  _subjets_ak8MotherIdx.clear();

  _RefitPVBS_xError.clear();
  _RefitPVBS_yError.clear();
  _RefitPVBS_zError.clear();

  _RefitPVBS_x.clear();
  _RefitPVBS_y.clear();
  _RefitPVBS_z.clear();

  _RefitPVBS_Cov.clear();

  _RefitPVNoBS_x.clear();
  _RefitPVNoBS_y.clear();
  _RefitPVNoBS_z.clear();

  _RefitPVNoBS_xError.clear();
  _RefitPVNoBS_yError.clear();
  _RefitPVNoBS_zError.clear();

  _VertexHashNoBS1.clear();
  _VertexHashNoBS2.clear();
  _VertexHashBS1.clear();
  _VertexHashBS2.clear();
  _VertexHashNoBSTracksRemovedOld1.clear();
  _VertexHashNoBSTracksRemovedOld2.clear();
  _VertexHashBSTracksRemovedOld1.clear();
  _VertexHashBSTracksRemovedOld2.clear();

  _LeptonHash.clear();

  _pvRefit_cov.clear();
  //_genH_px.clear();
  //_genH_py.clear();
  //_genH_pz.clear();
  //_genH_e.clear();

  // not a tree var, but has to be filled once per daughter - reset here
  vTrgMatchedToDau_idx.clear();


  PFTau_a1_lvp.clear();
  PFTau_a1_cov.clear();
  PFTau_a1_charge.clear();
  PFTau_a1_pdgid.clear();
  PFTau_a1_B.clear();
  PFTau_a1_M.clear();
  A1LVP.clear();

  Muon_trackCharge.clear();
  Muon_pdgid.clear();
  Muon_B.clear();
  Muon_M.clear();
  Muon_par.clear();
  Muon_cov.clear();
  MuonTrack.clear();

  PFTau_Track_par.clear();
  PFTau_Track_cov.clear();
  PFTau_Track_charge.clear();
  PFTau_Track_pdgid.clear();
  PFTau_Track_B.clear();
  PFTau_Track_M.clear();
  PFTauLeadTrackLV.clear();
  PFTauTrack_deltaR.clear();
  TauFLSignificance.clear();
  PFTauGEOMFlightLenght.clear();
  PFTauGEOMFlightLenghtSignificance.clear();


  // //Synchro

  //dilepton_veto=false;
  eleveto=false;
  muonveto=false;
  trg_singleelectron.clear();
  trg_singlemuon.clear();
  trg_singletau.clear();
  trg_muonelectron.clear();
  trg_mutaucross.clear();
  trg_doubletau.clear();

  pt_1.clear();
  phi_1.clear();
  eta_1.clear();
  m_1.clear();
  q_1.clear();
  d0_1.clear();
  dz_1.clear();
  mt_1.clear();
  //pfmt_1.clear();
  puppimt_1.clear();
  iso_1.clear();
  gen_match_1.clear();
  againstElectronLooseMVA6_1.clear();
  againstElectronMediumMVA6_1.clear();
  againstElectronTightMVA6_1.clear();
  againstElectronVLooseMVA6_1.clear();
  againstElectronVTightMVA6_1.clear();
  againstMuonLoose3_1.clear();
  againstMuonTight3_1.clear();
  byIsolationMVA3oldDMwLTraw_1.clear();

  idisoweight_1.clear();
  antieweight_1.clear();
  antimuweight_1.clear();

  pt_2.clear();
  phi_2.clear();
  eta_2.clear();
  m_2.clear();
  q_2.clear();
  d0_2.clear();
  dz_2.clear();
  mt_2.clear();
  //pfmt_2.clear();
  puppimt_2.clear();
  iso_2.clear();
  gen_match_2.clear();
  againstElectronLooseMVA6_2.clear();
  againstElectronMediumMVA6_2.clear();
  againstElectronTightMVA6_2.clear();
  againstElectronVLooseMVA6_2.clear();
  againstElectronVTightMVA6_2.clear();
  againstMuonLoose3_2.clear();
  againstMuonTight3_2.clear();
  byIsolationMVA3oldDMwLTraw_2.clear();

  idisoweight_2.clear();
  antieweight_2.clear();
  antimuweight_2.clear();

  pt_tt.clear();
  pt_vis.clear();
  mt_tot.clear();
  m_vis.clear();

  _PUPPImet=-99;
  _PUPPImetphi=-99;
  pzetavis.clear();
  pzetamiss.clear();
  //pfpzetamiss.clear();
  puppipzetamiss.clear();
  // _PFMETCov00=-99;
  // _PFMETCov01=-99;
  // _PFMETCov10=-99;
  // _PFMETCov11=-99;
  _PUPPIMETCov00=-99;
  _PUPPIMETCov01=-99;
  _PUPPIMETCov10=-99;
  _PUPPIMETCov11=-99;
  _puppimet_ex_JetEnUp=-99;
  _puppimet_ey_JetEnUp=-99;
  _puppimet_ex_JetEnDown=-99;
  _puppimet_ey_JetEnDown=-99;
  _puppimet_ex_UnclusteredEnUp=-99;
  _puppimet_ey_UnclusteredEnUp=-99;
  _puppimet_ex_UnclusteredEnDown=-99;
  _puppimet_ey_UnclusteredEnDown=-99;
  _puppimet_ex_JetResUp=-99;
  _puppimet_ey_JetResUp=-99;
  _puppimet_ex_JetResDown=-99;
  _puppimet_ey_JetResDown=-99;
  mjj.clear();
  jdeta.clear();
  njetingap.clear();
  njetingap20.clear();
  jdphi.clear();
  dijetpt.clear();
  dijetphi.clear();
  ptvis.clear();

  nbtag.clear();
  njets.clear();
  njetspt20.clear();
  jpt_1.clear();
  jeta_1.clear();
  jphi_1.clear();
  jcsv_1.clear();
  jpt_2.clear();
  jeta_2.clear();
  jphi_2.clear();
  jcsv_2.clear();
  bpt_1.clear();
  beta_1.clear();
  bphi_1.clear();
  bcsv_1.clear();
  bpt_2.clear();
  beta_2.clear();
  bphi_2.clear();
  bcsv_2.clear();

  mjjDown.clear();
  jdetaDown.clear();
  njetingapDown.clear();
  njetingap20Down.clear();
  jdphiDown.clear();
  dijetptDown.clear();
  dijetphiDown.clear();
  ptvisDown.clear();

  nbtagDown.clear();
  njetsDown.clear();
  njetspt20Down.clear();
  jptDown_1.clear();
  jetaDown_1.clear();
  jphiDown_1.clear();
  jcsvDown_1.clear();
  jptDown_2.clear();
  jetaDown_2.clear();
  jphiDown_2.clear();
  jcsvDown_2.clear();
  bptDown_1.clear();
  betaDown_1.clear();
  bphiDown_1.clear();
  bcsvDown_1.clear();
  bptDown_2.clear();
  betaDown_2.clear();
  bphiDown_2.clear();
  bcsvDown_2.clear();

  mjjUp.clear();
  jdetaUp.clear();
  njetingapUp.clear();
  njetingap20Up.clear();
  jdphiUp.clear();
  dijetptUp.clear();
  dijetphiUp.clear();
  ptvisUp.clear();

  nbtagUp.clear();
  njetsUp.clear();
  njetspt20Up.clear();
  jptUp_1.clear();
  jetaUp_1.clear();
  jphiUp_1.clear();
  jcsvUp_1.clear();
  jptUp_2.clear();
  jetaUp_2.clear();
  jphiUp_2.clear();
  jcsvUp_2.clear();
  bptUp_1.clear();
  betaUp_1.clear();
  bphiUp_1.clear();
  bcsvUp_1.clear();
  bptUp_2.clear();
  betaUp_2.clear();
  bphiUp_2.clear();
  bcsvUp_2.clear();

  puweight=1;

  _nup=-99;

  weight.clear();

  jpfid_1.clear();
  jpuid_1.clear();
  jpfid_2.clear();
  jpuid_2.clear();
  bpfid_1.clear();
  bpuid_1.clear();
  bpfid_2.clear();
  bpuid_2.clear();

  _npv=-99;
  _npu=-99;
  _rho=-99;

  pt_sv.clear();
  eta_sv.clear();
  phi_sv.clear();
  met_sv.clear();

  byDeepTau2017v2p1VSjetraw_1.clear();
  byVVVLooseDeepTau2017v2p1VSjet_1.clear();
  byVVLooseDeepTau2017v2p1VSjet_1.clear();
  byVLooseDeepTau2017v2p1VSjet_1.clear();
  byLooseDeepTau2017v2p1VSjet_1.clear();
  byMediumDeepTau2017v2p1VSjet_1.clear();
  byTightDeepTau2017v2p1VSjet_1.clear();
  byVTightDeepTau2017v2p1VSjet_1.clear();
  byVVTightDeepTau2017v2p1VSjet_1.clear();

  byDeepTau2017v2p1VSeleraw_1.clear();
  byVVVLooseDeepTau2017v2p1VSe_1.clear();
  byVVLooseDeepTau2017v2p1VSe_1.clear();
  byVLooseDeepTau2017v2p1VSe_1.clear();
  byLooseDeepTau2017v2p1VSe_1.clear();
  byMediumDeepTau2017v2p1VSe_1.clear();
  byTightDeepTau2017v2p1VSe_1.clear();
  byVTightDeepTau2017v2p1VSe_1.clear();
  byVVTightDeepTau2017v2p1VSe_1.clear();

  byDeepTau2017v2p1VSmuraw_1.clear();
  byVLooseDeepTau2017v2p1VSmu_1.clear();
  byLooseDeepTau2017v2p1VSmu_1.clear();
  byMediumDeepTau2017v2p1VSmu_1.clear();
  byTightDeepTau2017v2p1VSmu_1.clear();

  byDeepTau2017v2p1VSjetraw_2.clear();
  byVVVLooseDeepTau2017v2p1VSjet_2.clear();
  byVVLooseDeepTau2017v2p1VSjet_2.clear();
  byVLooseDeepTau2017v2p1VSjet_2.clear();
  byLooseDeepTau2017v2p1VSjet_2.clear();
  byMediumDeepTau2017v2p1VSjet_2.clear();
  byTightDeepTau2017v2p1VSjet_2.clear();
  byVTightDeepTau2017v2p1VSjet_2.clear();
  byVVTightDeepTau2017v2p1VSjet_2.clear();

  byDeepTau2017v2p1VSeleraw_2.clear();
  byVVVLooseDeepTau2017v2p1VSe_2.clear();
  byVVLooseDeepTau2017v2p1VSe_2.clear();
  byVLooseDeepTau2017v2p1VSe_2.clear();
  byLooseDeepTau2017v2p1VSe_2.clear();
  byMediumDeepTau2017v2p1VSe_2.clear();
  byTightDeepTau2017v2p1VSe_2.clear();
  byVTightDeepTau2017v2p1VSe_2.clear();
  byVVTightDeepTau2017v2p1VSe_2.clear();

  byDeepTau2017v2p1VSmuraw_2.clear();
  byVLooseDeepTau2017v2p1VSmu_2.clear();
  byLooseDeepTau2017v2p1VSmu_2.clear();
  byMediumDeepTau2017v2p1VSmu_2.clear();
  byTightDeepTau2017v2p1VSmu_2.clear();

  pvx = 0;
  pvy = 0;
  pvz = 0;
  pvcov.clear();
  pcax.clear();
  pcay.clear();
  pcaz.clear();

  dm_1.clear();
  dmMVA_1.clear();

  dm_2.clear();
  dmMVA_2.clear();

  _chpt_1.clear();
  _chphi_1.clear();
  _cheta_1.clear();
  _chm_1.clear();

  _chpt_2.clear();
  _chphi_2.clear();
  _cheta_2.clear();
  _chm_2.clear();

  _npt_1.clear();
  _nphi_1.clear();
  _neta_1.clear();
  _nm_1.clear();

  _npt_2.clear();
  _nphi_2.clear();
  _neta_2.clear();
  _nm_2.clear();

  svx_1.clear();
  svy_1.clear();
  svz_1.clear();

  svx_2.clear();
  svy_2.clear();
  svz_2.clear();

  //tauspinnerH=-99;
  //tauspinnerA=-99;
  //tauspinnerMaxMix=-99;

  LeptonP4.clear();
  tau1IndexVect.clear();
  tau2IndexVect.clear();

}

void HTauTauNtuplizer::beginJob(){

  edm::Service<TFileService> fs;
  //myTree = fs->make<TTree>("HTauTauTree","HTauTauTree");
  for(std::map<std::string, TTree*>::iterator iMap = SystematicsMap.begin(); iMap != SystematicsMap.end(); ++iMap) {
    iMap->second = fs->make<TTree>((iMap->first).c_str(), (iMap->first).c_str());
    bool isNominal = false;
    if(iMap->first == "Nominal") isNominal = true;
    mySysHelper->MakeBranches(iMap->second, isNominal);
  }

  int nbins=3+(myTriggerHelper->GetNTriggers());
  hCounter = fs->make<TH1F>("Counters","Counters",nbins,0,nbins);
  hTauIDs = fs->make<TH1F>("TauIDs","TauIDs",ntauIds,0,ntauIds);
  hYear = fs->make<TH1F>("Year","Year",50,2000,2050);

  //Branches
  //myTree->Branch("DataMC_Type_idx" ,&_DataMC_Type);
  //  myTree->Branch("evt",&_indexevents,"EventNumber/l");
  //  myTree->Branch("run",&_runNumber,"RunNumber/I");
  //  myTree->Branch("lumi",&_lumi,"lumi/I");
 // myTree->Branch("year",&_year,"year/I");
  //  myTree->Branch("NBadMu",&_NBadMu,"NBadMu/I");
 // myTree->Branch("passecalBadCalibFilterUpdate",&_passecalBadCalibFilterUpdate);
 // myTree->Branch("prefiringweight",&_prefiringweight,"prefiringweight/F");
 // myTree->Branch("prefiringweightup",&_prefiringweightup,"prefiringweightup/F");
 // myTree->Branch("prefiringweightdown",&_prefiringweightdown,"prefiringweightdown/F");
 // myTree->Branch("triggerbit",&_triggerbit,"triggerbit/L");
 // myTree->Branch("metfilterbit",&_metfilterbit,"metfilterbit/I");
  //myTree->Branch("met",&_met,"met/F");
  //myTree->Branch("met_er",&_met_er,"met_er/F");
  //myTree->Branch("met_er_phi",&_met_er_phi,"met_er_phi/F");
  /*if(DETAIL>=1){
    myTree->Branch("daughters_IetaIeta",&_daughters_IetaIeta);
    myTree->Branch("daughters_full5x5_IetaIeta",&_daughters_full5x5_IetaIeta);
    myTree->Branch("daughters_hOverE",&_daughters_hOverE);
    myTree->Branch("daughters_deltaEtaSuperClusterTrackAtVtx",&_daughters_deltaEtaSuperClusterTrackAtVtx);
    myTree->Branch("daughters_deltaPhiSuperClusterTrackAtVtx",&_daughters_deltaPhiSuperClusterTrackAtVtx);
    myTree->Branch("daughters_IoEmIoP",&_daughters_IoEmIoP);
    myTree->Branch("daughters_IoEmIoP_ttH",&_daughters_IoEmIoP_ttH);
  }*/
  // myTree->Branch("PFMETCov00",&_PFMETCov00,"PFMETCov00/F");
  // myTree->Branch("PFMETCov01",&_PFMETCov01,"PFMETCov01/F");
  // myTree->Branch("PFMETCov10",&_PFMETCov10,"PFMETCov10/F");
  // myTree->Branch("PFMETCov11",&_PFMETCov11,"PFMETCov11/F");
  // myTree->Branch("PFMETsignif", &_PFMETsignif, "PFMETsignif/F");

  //  myTree->Branch("npv",&_npv,"npv/I");
  //  myTree->Branch("npu",&_npu,"npu/F");
  //  //  myTree->Branch("PUReweight",&_PUReweight,"PUReweight/F");
  //  myTree->Branch("rho",&_rho,"rho/F");

  /*myTree->Branch("mothers_px",&_mothers_px);
  myTree->Branch("mothers_py",&_mothers_py);
  myTree->Branch("mothers_pz",&_mothers_pz);
  myTree->Branch("mothers_e",&_mothers_e);
  myTree->Branch("mothers_trgSeparateMatch", &_mothers_trgSeparateMatch);
  if(DEBUG){
    myTree->Branch("trigger_name",&_trigger_name);
    myTree->Branch("trigger_accept",&_trigger_accept);
  }
  myTree->Branch("daughters_px",&_daughters_px);
  myTree->Branch("daughters_py",&_daughters_py);
  myTree->Branch("daughters_pz",&_daughters_pz);
  myTree->Branch("daughters_e",&_daughters_e);
  myTree->Branch("daughters_charge",&_daughters_charge);
  if(doCPVariables){
    myTree->Branch("daughters_charged_px",&_daughters_charged_px);
    myTree->Branch("daughters_charged_py",&_daughters_charged_py);
    myTree->Branch("daughters_charged_pz",&_daughters_charged_pz);
    myTree->Branch("daughters_charged_e",&_daughters_charged_e);

    myTree->Branch("daughters_neutral_px",&_daughters_neutral_px);
    myTree->Branch("daughters_neutral_py",&_daughters_neutral_py);
    myTree->Branch("daughters_neutral_pz",&_daughters_neutral_pz);
    myTree->Branch("daughters_neutral_e",&_daughters_neutral_e);

    myTree->Branch("daughters_vx",&_daughters_vx);
    myTree->Branch("daughters_vy",&_daughters_vy);
    myTree->Branch("daughters_vz",&_daughters_vz);
  }


  if(writeSoftLep)myTree->Branch("softLeptons",&_softLeptons);
  myTree->Branch("L1_tauEt",&_L1_tauEt);
  myTree->Branch("L1_tauEta",&_L1_tauEta);
  myTree->Branch("L1_tauPhi",&_L1_tauPhi);
  myTree->Branch("L1_tauIso",&_L1_tauIso);
  myTree->Branch("L1_jetEt",&_L1_jetEt);
  myTree->Branch("L1_jetEta",&_L1_jetEta);
  myTree->Branch("L1_jetPhi",&_L1_jetPhi);

  myTree->Branch("daughters_highestEt_L1IsoTauMatched", &_daughters_highestEt_L1IsoTauMatched);
  if( theisMC || IsEmbed){

    myTree->Branch("daughters_hasTES",&_daughters_hasTES);
    myTree->Branch("daughters_px_TauUp",&_daughters_px_TauUp);
    myTree->Branch("daughters_py_TauUp",&_daughters_py_TauUp);
    myTree->Branch("daughters_pz_TauUp",&_daughters_pz_TauUp);
    myTree->Branch("daughters_e_TauUp",&_daughters_e_TauUp);
    myTree->Branch("daughters_px_TauDown",&_daughters_px_TauDown);
    myTree->Branch("daughters_py_TauDown",&_daughters_py_TauDown);
    myTree->Branch("daughters_pz_TauDown",&_daughters_pz_TauDown);
    myTree->Branch("daughters_e_TauDown",&_daughters_e_TauDown);
    myTree->Branch("daughters_hasEES",&_daughters_hasEES);
    myTree->Branch("daughters_px_EleUp",&_daughters_px_EleUp);
    myTree->Branch("daughters_py_EleUp",&_daughters_py_EleUp);
    myTree->Branch("daughters_pz_EleUp",&_daughters_pz_EleUp);
    myTree->Branch("daughters_e_EleUp",&_daughters_e_EleUp);
    myTree->Branch("daughters_px_EleDown",&_daughters_px_EleDown);
    myTree->Branch("daughters_py_EleDown",&_daughters_py_EleDown);
    myTree->Branch("daughters_pz_EleDown",&_daughters_pz_EleDown);
    myTree->Branch("daughters_e_EleDown",&_daughters_e_EleDown);

    myTree->Branch("daughters_isTauMatched",&_daughters_isTauMatched);
    //myTree->Branch("genDaughters",&_genDaughters);
    //myTree->Branch("quarks_px",&_bquarks_px);
    //myTree->Branch("quarks_py",&_bquarks_py);
    //myTree->Branch("quarks_pz",&_bquarks_pz);
    //myTree->Branch("quarks_e",&_bquarks_e);
    //myTree->Branch("quarks_pdg",&_bquarks_pdg);
    //myTree->Branch("motmass",&_bmotmass);
    //myTree->Branch("genH_px",&_genH_px);
    //myTree->Branch("genH_py",&_genH_py);
    //myTree->Branch("genH_pz",&_genH_pz);
    //myTree->Branch("genH_e",&_genH_e);
    myTree->Branch("PUNumInteractions",&_PUNumInteractions,"PUNumInteractions/I");
    myTree->Branch("daughters_genindex",&_daughters_genindex);
    myTree->Branch("MC_weight",&_MC_weight);
    myTree->Branch("_nominal_wt",&_nominal_wt);
    myTree->Branch("TheoreticalPSUnc",&_TheoreticalPSUnc);
    myTree->Branch("TheoreticalScaleUnc1"  , &_TheoreticalScaleUncTab[0]);
    myTree->Branch("TheoreticalScaleUnc2"  , &_TheoreticalScaleUncTab[1]);
    myTree->Branch("TheoreticalScaleUnc3"  , &_TheoreticalScaleUncTab[2]);
    myTree->Branch("TheoreticalScaleUnc4"  , &_TheoreticalScaleUncTab[3]);
    myTree->Branch("TheoreticalScaleUnc5"  , &_TheoreticalScaleUncTab[4]);
    myTree->Branch("TheoreticalScaleUnc6"  , &_TheoreticalScaleUncTab[5]);
    myTree->Branch("TheoreticalScaleUnc7"  , &_TheoreticalScaleUncTab[6]);
    myTree->Branch("TheoreticalScaleUnc8"  , &_TheoreticalScaleUncTab[7]);
    myTree->Branch("TheoreticalScaleUnc9"  , &_TheoreticalScaleUncTab[8]);
    myTree->Branch("TheoreticalScaleUnc1001"  , &_TheoreticalScaleUncTab[9]);
    myTree->Branch("TheoreticalScaleUnc1002"  , &_TheoreticalScaleUncTab[10]);
    myTree->Branch("TheoreticalScaleUnc1003"  , &_TheoreticalScaleUncTab[11]);
    myTree->Branch("TheoreticalScaleUnc1004"  , &_TheoreticalScaleUncTab[12]);
    myTree->Branch("TheoreticalScaleUnc1005"  , &_TheoreticalScaleUncTab[13]);
    myTree->Branch("TheoreticalScaleUnc1006"  , &_TheoreticalScaleUncTab[14]);
    myTree->Branch("TheoreticalScaleUnc1007"  , &_TheoreticalScaleUncTab[15]);
    myTree->Branch("TheoreticalScaleUnc1008"  , &_TheoreticalScaleUncTab[16]);
    myTree->Branch("TheoreticalScaleUnc1009"  , &_TheoreticalScaleUncTab[17]);

    myTree->Branch("MC_weight_scale_muF0p5",&_MC_weight_scale_muF0p5);
    myTree->Branch("MC_weight_scale_muF2",&_MC_weight_scale_muF2);
    myTree->Branch("MC_weight_scale_muR0p5",&_MC_weight_scale_muR0p5);
    myTree->Branch("MC_weight_scale_muR2",&_MC_weight_scale_muR2);
    myTree->Branch("lheHt",&_lheHt,"lheHt/F");
    myTree->Branch("lheNOutPartons", &_lheNOutPartons, "lheNOutPartons/I");
    myTree->Branch("lheNOutB", &_lheNOutB, "lheNOutB/I");
    myTree->Branch("lheNOutC", &_lheNOutC, "lheNOutC/I");
    myTree->Branch("aMCatNLOweight",&_aMCatNLOweight);
    myTree->Branch("genpart_px", &_genpart_px);
    myTree->Branch("genpart_py", &_genpart_py);
    myTree->Branch("genpart_pz", &_genpart_pz);
    myTree->Branch("genpart_e", &_genpart_e);


    if(doCPVariables){
      myTree->Branch("genpart_pca_x",&_genpart_pca_x);
      myTree->Branch("genpart_pca_y",&_genpart_pca_y);
      myTree->Branch("genpart_pca_z",&_genpart_pca_z);
    }
    myTree->Branch("genpart_pdg", &_genpart_pdg);
    myTree->Branch("genpart_status", &_genpart_status);
    //myTree->Branch("genpart_mothInd", _genpart_mothInd);
    myTree->Branch("genpart_HMothInd", &_genpart_HMothInd);
    myTree->Branch("genpart_MSSMHMothInd", &_genpart_MSSMHMothInd);
    myTree->Branch("genpart_TopMothInd", &_genpart_TopMothInd);
    myTree->Branch("genpart_TauMothInd", &_genpart_TauMothInd);
    myTree->Branch("genpart_ZMothInd", &_genpart_ZMothInd);
    myTree->Branch("genpart_WMothInd", &_genpart_WMothInd);
    myTree->Branch("genpart_bMothInd", &_genpart_bMothInd);
    myTree->Branch("genpart_HZDecayMode", &_genpart_HZDecayMode);
    myTree->Branch("genpart_TopDecayMode", &_genpart_TopDecayMode);
    myTree->Branch("genpart_WDecayMode", &_genpart_WDecayMode);
    myTree->Branch("genpart_TauGenDecayMode", &_genpart_TauGenDecayMode);
    myTree->Branch("genpart_TauGenDetailedDecayMode", &_genpart_TauGenDetailedDecayMode);
    myTree->Branch("genpart_flags", &_genpart_flags);

    myTree->Branch("genjet_px", &_genjet_px);
    myTree->Branch("genjet_py", &_genjet_py);
    myTree->Branch("genjet_pz", &_genjet_pz);
    myTree->Branch("genjet_e" , &_genjet_e);
    myTree->Branch("genjet_partonFlavour" , &_genjet_partonFlavour);
    myTree->Branch("genjet_hadronFlavour" , &_genjet_hadronFlavour);

    //    myTree->Branch("NUP", &_nup,"NUP/I");
    // myTree->Branch("SVfit_fitMETPhiTauUp", &_SVMetPhiTauUp);
    // myTree->Branch("SVfit_fitMETPhiTauDown", &_SVMetPhiTauDown);
    // myTree->Branch("SVfit_fitMETRhoTauUp", &_SVMetRhoTauUp);
    // myTree->Branch("SVfit_fitMETRhoTauDown", &_SVMetRhoTauDown);
    // myTree->Branch("SVfit_phiUncTauUp", &_SVphiUncTauUp);
    // myTree->Branch("SVfit_phiUncTauDown", &_SVphiUncTauDown);
    // myTree->Branch("SVfit_phiTauUp", &_SVphiTauUp);
    // myTree->Branch("SVfit_phiTauDown", &_SVphiTauDown);
    // myTree->Branch("SVfit_etaUncTauUp", &_SVetaUncTauUp);
    // myTree->Branch("SVfit_etaUncTauDown", &_SVetaUncTauDown);
    // myTree->Branch("SVfit_etaTauUp", &_SVetaTauUp);
    // myTree->Branch("SVfit_etaTauDown", &_SVetaTauDown);
    // myTree->Branch("SVfit_ptUncTauUp", &_SVptUncTauUp);
    // myTree->Branch("SVfit_ptUncTauDown", &_SVptUncTauDown);
    // myTree->Branch("SVfit_ptTauUp", &_SVptTauUp);
    // myTree->Branch("SVfit_ptTauDown", &_SVptTauDown);
    // myTree->Branch("SVfitTransverseMassUncTauUp",&_SVmassTransverseUncTauUp);
    // myTree->Branch("SVfitTransverseMassUncTauDown",&_SVmassTransverseUncTauDown);
    // myTree->Branch("SVfitTransverseMassTauUp",&_SVmassTransverseTauUp);
    // myTree->Branch("SVfitTransverseMassTauDown",&_SVmassTransverseTauDown);
    // myTree->Branch("SVfitMassUncTauUp",&_SVmassUncTauUp);
    // myTree->Branch("SVfitMassUncTauDown",&_SVmassUncTauDown);
    // myTree->Branch("SVfitMassTauUp",&_SVmassTauUp);
    // myTree->Branch("SVfitMassTauDown",&_SVmassTauDown);

    // myTree->Branch("SVfit_fitMETPhiMETUp", &_SVMetPhiMETUp);
    // myTree->Branch("SVfit_fitMETPhiMETDown", &_SVMetPhiMETDown);
    // myTree->Branch("SVfit_fitMETRhoMETUp", &_SVMetRhoMETUp);
    // myTree->Branch("SVfit_fitMETRhoMETDown", &_SVMetRhoMETDown);
    // myTree->Branch("SVfit_phiUncMETUp", &_SVphiUncMETUp);
    // myTree->Branch("SVfit_phiUncMETDown", &_SVphiUncMETDown);
    // myTree->Branch("SVfit_phiMETUp", &_SVphiMETUp);
    // myTree->Branch("SVfit_phiMETDown", &_SVphiMETDown);
    // myTree->Branch("SVfit_etaUncMETUp", &_SVetaUncMETUp);
    // myTree->Branch("SVfit_etaUncMETDown", &_SVetaUncMETDown);
    // myTree->Branch("SVfit_etaMETUp", &_SVetaMETUp);
    // myTree->Branch("SVfit_etaMETDown", &_SVetaMETDown);
    // myTree->Branch("SVfit_ptUncMETUp", &_SVptUncMETUp);
    // myTree->Branch("SVfit_ptUncMETDown", &_SVptUncMETDown);
    // myTree->Branch("SVfit_ptMETUp", &_SVptMETUp);
    // myTree->Branch("SVfit_ptMETDown", &_SVptMETDown);
    // myTree->Branch("SVfitTransverseMassUncMETUp",&_SVmassTransverseUncMETUp);
    // myTree->Branch("SVfitTransverseMassUncMETDown",&_SVmassTransverseUncMETDown);
    // myTree->Branch("SVfitTransverseMassMETUp",&_SVmassTransverseMETUp);
    // myTree->Branch("SVfitTransverseMassMETDown",&_SVmassTransverseMETDown);
    // myTree->Branch("SVfitMassUncMETUp",&_SVmassUncMETUp);
    // myTree->Branch("SVfitMassUncMETDown",&_SVmassUncMETDown);
    // myTree->Branch("SVfitMassMETUp",&_SVmassMETUp);
    // myTree->Branch("SVfitMassMETDown",&_SVmassMETDown);

    // myTree->Branch("SVfit_fitMETPhiEleUp", &_SVMetPhiEleUp);
    // myTree->Branch("SVfit_fitMETPhiEleDown", &_SVMetPhiEleDown);
    // myTree->Branch("SVfit_fitMETRhoEleUp", &_SVMetRhoEleUp);
    // myTree->Branch("SVfit_fitMETRhoEleDown", &_SVMetRhoEleDown);
    // myTree->Branch("SVfit_phiUncEleUp", &_SVphiUncEleUp);
    // myTree->Branch("SVfit_phiUncEleDown", &_SVphiUncEleDown);
    // myTree->Branch("SVfit_phiEleUp", &_SVphiEleUp);
    // myTree->Branch("SVfit_phiEleDown", &_SVphiEleDown);
    // myTree->Branch("SVfit_etaUncEleUp", &_SVetaUncEleUp);
    // myTree->Branch("SVfit_etaUncEleDown", &_SVetaUncEleDown);
    // myTree->Branch("SVfit_etaEleUp", &_SVetaEleUp);
    // myTree->Branch("SVfit_etaEleDown", &_SVetaEleDown);
    // myTree->Branch("SVfit_ptUncEleUp", &_SVptUncEleUp);
    // myTree->Branch("SVfit_ptUncEleDown", &_SVptUncEleDown);
    // myTree->Branch("SVfit_ptEleUp", &_SVptEleUp);
    // myTree->Branch("SVfit_ptEleDown", &_SVptEleDown);
    // myTree->Branch("SVfitTransverseMassUncEleUp",&_SVmassTransverseUncEleUp);
    // myTree->Branch("SVfitTransverseMassUncEleDown",&_SVmassTransverseUncEleDown);
    // myTree->Branch("SVfitTransverseMassEleUp",&_SVmassTransverseEleUp);
    // myTree->Branch("SVfitTransverseMassEleDown",&_SVmassTransverseEleDown);
    // myTree->Branch("SVfitMassUncEleUp",&_SVmassUncEleUp);
    // myTree->Branch("SVfitMassUncEleDown",&_SVmassUncEleDown);
    // myTree->Branch("SVfitMassEleUp",&_SVmassEleUp);
    // myTree->Branch("SVfitMassEleDown",&_SVmassEleDown);
  }// end if isMC
  //myTree->Branch("daughters2",&_daughter2);
  if (do_MCSummary_ || IsEmbed) {
    myTree->Branch("MCSignalParticle_p4", &MCSignalParticle_p4);
    myTree->Branch("MCSignalParticle_pdgid", &MCSignalParticle_pdgid);
    myTree->Branch("MCSignalParticle_charge", &MCSignalParticle_charge);
    myTree->Branch("MCSignalParticle_Poca", &MCSignalParticle_Poca);
    myTree->Branch("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx);
    myTree->Branch("MCTauandProd_p4", &MCTauandProd_p4);
    myTree->Branch("MCTauandProd_Vertex", &MCTauandProd_Vertex);

    myTree->Branch("MCTauandProd_pdgid", &MCTauandProd_pdgid);
    myTree->Branch("MCTauandProd_midx", &MCTauandProd_midx);
    myTree->Branch("MCTauandProd_charge", &MCTauandProd_charge);
    myTree->Branch("MCTau_JAK", &MCTau_JAK);
    myTree->Branch("MCTau_DecayBitMask", &MCTau_DecayBitMask);
  }

  if (do_MCComplete_ || IsEmbed) {
    myTree->Branch("MC_p4", &MC_p4);
    myTree->Branch("MC_pdgid", &MC_pdgid);
    myTree->Branch("MC_charge", &MC_charge);
    myTree->Branch("MC_midx", &MC_midx);
    myTree->Branch("MC_childpdgid", &MC_childpdgid);
    myTree->Branch("MC_childidx", &MC_childidx);
    myTree->Branch("MC_status", &MC_status);
  }

  myTree->Branch("Event_isRealData", &Event_isRealData);
  // myTree->Branch("SVfitMass",&_SVmass);
  // myTree->Branch("SVfitMassUnc",&_SVmassUnc);
  // myTree->Branch("SVfitTransverseMass",&_SVmassTransverse);
  // myTree->Branch("SVfitTransverseMassUnc",&_SVmassTransverseUnc);
  // myTree->Branch("SVfit_pt", &_SVpt);
  // myTree->Branch("SVfit_ptUnc", &_SVptUnc);
  // myTree->Branch("SVfit_eta", &_SVeta);
  // myTree->Branch("SVfit_etaUnc", &_SVetaUnc);
  // myTree->Branch("SVfit_phi", &_SVphi);
  // myTree->Branch("SVfit_phiUnc", &_SVphiUnc);
  // myTree->Branch("SVfit_fitMETRho", &_SVMetRho);
  // myTree->Branch("SVfit_fitMETPhi", &_SVMetPhi);
  myTree->Branch("isOSCand",&_isOSCand);
  myTree->Branch("METx",&_metx);
  myTree->Branch("METy",&_mety);
  // myTree->Branch("METx_UP_JES",&_metx_up_jes);
  // myTree->Branch("METy_UP_JES",&_mety_up_jes);
  // myTree->Branch("METx_DOWN_JES",&_metx_down_jes);
  // myTree->Branch("METy_DOWN_JES",&_mety_down_jes);
  // myTree->Branch("METx_UP_TES",&_metx_up_tes);
  // myTree->Branch("METy_UP_TES",&_mety_up_tes);
  // myTree->Branch("METx_DOWN_TES",&_metx_down_tes);
  // myTree->Branch("METy_DOWN_TES",&_mety_down_tes);
  // myTree->Branch("METx_UP_EES",&_metx_up_ees);
  // myTree->Branch("METy_UP_EES",&_mety_up_ees);
  // myTree->Branch("METx_DOWN_EES",&_metx_down_ees);
  // myTree->Branch("METy_DOWN_EES",&_mety_down_ees);
  //myTree->Branch("uncorrMETx",&_uncorrmetx);
  //myTree->Branch("uncorrMETy",&_uncorrmety);
  // myTree->Branch("MET_cov00",&_metCov00);
  // myTree->Branch("MET_cov01",&_metCov01);
  // myTree->Branch("MET_cov10",&_metCov10);
  // myTree->Branch("MET_cov11",&_metCov11);
  //  myTree->Branch("MET_significance",&_metSignif);
  myTree->Branch("mT_Dau1",&_mTDau1);
  myTree->Branch("mT_Dau2",&_mTDau2);
  myTree->Branch("PDGIdDaughters",&_pdgdau);
  myTree->Branch("indexDau1",&_indexDau1);
  myTree->Branch("indexDau2",&_indexDau2);
  myTree->Branch("particleType",&_particleType);
  //myTree->Branch("discriminator",&_discriminator);
  myTree->Branch("daughters_muonID",&_daughters_muonID);
  myTree->Branch("daughters_typeOfMuon",&_daughters_typeOfMuon);
  myTree->Branch("dxy",&_dxy);
  myTree->Branch("dz",&_dz);
  myTree->Branch("dxy_innerTrack",&_dxy_innerTrack);
  myTree->Branch("dz_innerTrack",&_dz_innerTrack);
  myTree->Branch("daughters_rel_error_trackpt",&_daughters_rel_error_trackpt);
  myTree->Branch("SIP",&_SIP);
  //  myTree->Branch("daughters_iseleBDT",&_daughters_iseleBDT);
  myTree->Branch("daughters_iseleWPLoose",&_daughters_iseleWPLoose);
  myTree->Branch("daughters_iseleWP80",&_daughters_iseleWP80);
  myTree->Branch("daughters_iseleWP90",&_daughters_iseleWP90);
  myTree->Branch("daughters_iseleNoIsoWPLoose",&_daughters_iseleNoIsoWPLoose);
  myTree->Branch("daughters_iseleNoIsoWP80",&_daughters_iseleNoIsoWP80);
  myTree->Branch("daughters_iseleNoIsoWP90",&_daughters_iseleNoIsoWP90);
  myTree->Branch("daughters_iseleCutBased",&_daughters_iseleCutBased);
  myTree->Branch("daughters_eleMVAnt",&_daughters_eleMVAnt);
  //  myTree->Branch("daughters_eleMVAntNoIso",&_daughters_eleMVAntNoIso);
  myTree->Branch("daughters_eleMVA_HZZ",&_daughters_eleMVA_HZZ);
  myTree->Branch("daughters_passConversionVeto",&_daughters_passConversionVeto);
  myTree->Branch("daughters_eleMissingHits",&_daughters_eleMissingHits);
  //myTree->Branch("daughters_eleMissingLostHits",&_daughters_eleMissingLostHits);
  myTree->Branch("daughters_iseleChargeConsistent",&_daughters_iseleChargeConsistent);
  //myTree->Branch("daughters_eleCUTID",&_daughters_iseleCUT);
  myTree->Branch("decayMode",&_decayType);
  myTree->Branch("genmatch",&_genmatch);
  myTree->Branch("tauID",&_daughters_tauID);
  myTree->Branch("combreliso",& _combreliso);
  myTree->Branch("combreliso03",& _combreliso03);
  myTree->Branch("daughters_depositR03_tracker",&_daughters_depositR03_tracker);
  myTree->Branch("daughters_depositR03_ecal",&_daughters_depositR03_ecal);
  myTree->Branch("daughters_depositR03_hcal",&_daughters_depositR03_hcal);
  myTree->Branch("daughters_decayModeFindingOldDMs", &_daughters_decayModeFindingOldDMs);
  //  myTree->Branch("daughters_SCeta",&_daughters_SCeta);

  //myTree->Branch("MVADM2016v1",&_MVADM2016v1);
  myTree->Branch("MVADM2017v1",&_MVADM2017v1);
  //myTree->Branch("MVADM2017v1InTauFiller",&_MVADM2017v1InTauFiller);

  myTree->Branch("footprintCorrection",&_daughters_footprintCorrection);
  myTree->Branch("neutralIsoPtSumWeight",&_daughters_neutralIsoPtSumWeight);
  myTree->Branch("photonPtSumOutsideSignalCone",&_daughters_photonPtSumOutsideSignalCone);

  myTree->Branch("daughters_decayModeFindingNewDMs", &_daughters_decayModeFindingNewDMs);
  myTree->Branch("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  myTree->Branch("daughters_byIsolationMVArun2v1DBoldDMwLTraw",&_daughters_byIsolationMVArun2v1DBoldDMwLTraw);
  myTree->Branch("daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017",&_daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017); //FRA
  myTree->Branch("daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017",&_daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017); //FRA
  myTree->Branch("daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017",&_daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017); //FRA
  myTree->Branch("daughters_byIsolationMVArun2v1DBoldDMwLTraw",&_daughters_byIsolationMVArun2v1DBoldDMwLTraw);
  myTree->Branch("daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &_daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017); //FRA
  myTree->Branch("daughters_byDeepTau2017v2p1VSjetraw",&_daughters_byDeepTau2017v2p1VSjetraw);
  myTree->Branch("daughters_byDeepTau2017v2p1VSeraw",&_daughters_byDeepTau2017v2p1VSeraw);
  myTree->Branch("daughters_byDeepTau2017v2p1VSmuraw",&_daughters_byDeepTau2017v2p1VSmuraw);
  myTree->Branch("daughters_chargedIsoPtSum", &_daughters_chargedIsoPtSum);
  myTree->Branch("daughters_neutralIsoPtSum", &_daughters_neutralIsoPtSum);
  myTree->Branch("daughters_puCorrPtSum", &_daughters_puCorrPtSum);
  myTree->Branch("daughters_numChargedParticlesSignalCone", &_daughters_numChargedParticlesSignalCone);
  myTree->Branch("daughters_numNeutralHadronsSignalCone", &_daughters_numNeutralHadronsSignalCone);
  myTree->Branch("daughters_numPhotonsSignalCone", &_daughters_numPhotonsSignalCone);
  myTree->Branch("daughters_daughters_numParticlesSignalCone", &_daughters_numParticlesSignalCone);
  myTree->Branch("daughters_numChargedParticlesIsoCone", &_daughters_numChargedParticlesIsoCone);
  myTree->Branch("daughters_numNeutralHadronsIsoCone", &_daughters_numNeutralHadronsIsoCone);
  myTree->Branch("daughters_numPhotonsIsoCone", &_daughters_numPhotonsIsoCone);
  myTree->Branch("daughters_numParticlesIsoCone", &_daughters_numParticlesIsoCone);
  myTree->Branch("daughters_leadChargedParticlePt", &_daughters_leadChargedParticlePt);
  myTree->Branch("daughters_trackRefPt", &_daughters_trackRefPt);
    myTree->Branch("daughters_isLastTriggerObjectforPath", &_daughters_LFtrigger);
  myTree->Branch("daughters_trgMatched", &_daughters_trgMatched);
  myTree->Branch("daughters_isTriggerObjectforPath", &_daughters_L3trigger);
  myTree->Branch("daughters_FilterFired",&_daughters_FilterFired);
  myTree->Branch("daughters_isGoodTriggerType",&_daughters_isGoodTriggerType);
  myTree->Branch("daughters_L3FilterFired",&_daughters_L3FilterFired);
  myTree->Branch("daughters_L3FilterFiredLast",&_daughters_L3FilterFiredLast);
  myTree->Branch("daughters_HLTpt",&_daughters_HLTpt);
    myTree->Branch("daughters_isL1IsoTau28Matched", &_daughters_isL1IsoTau28Matched);


  myTree->Branch("Muon_trackCharge", &Muon_trackCharge);
  myTree->Branch("Muon_pdgid", &Muon_pdgid);
  myTree->Branch("Muon_B", &Muon_B);
  myTree->Branch("Muon_M", &Muon_M);
  myTree->Branch("Muon_par", &Muon_par);
  myTree->Branch("Muon_cov", &Muon_cov);


  myTree->Branch("PFTau_Track_par", &PFTau_Track_par);
  myTree->Branch("PFTau_Track_cov", &PFTau_Track_cov);
  myTree->Branch("PFTau_Track_charge", &PFTau_Track_charge);
  myTree->Branch("PFTau_Track_pdgid", &PFTau_Track_pdgid);
  myTree->Branch("PFTau_Track_B", &PFTau_Track_B);
  myTree->Branch("PFTau_Track_M", &PFTau_Track_M);


  myTree->Branch("PFTauSVPos", &_PFTauSVPos);
  myTree->Branch("PFTauSVCov", &_PFTauSVCov);


  myTree->Branch("PFTauPionsP4", &_PFTauPionsP4);
  myTree->Branch("PFTauRefitPionsP4", &_PFTauRefitPionsP4);
  myTree->Branch("PFTauPionsCharge", &_PFTauPionsCharge);
  myTree->Branch("PFTauRefitPionsCharge", &_PFTauRefitPionsCharge);
  myTree->Branch("PFTauSVChi2NDofMatchingQuality", &_PFTauSVChi2NDofMatchingQuality);

  myTree->Branch("TauFLSignificance", &TauFLSignificance);
  myTree->Branch("PFTauGEOMFlightLenght", &PFTauGEOMFlightLenght);
  myTree->Branch("PFTauGEOMFlightLenghtSignificance", &PFTauGEOMFlightLenghtSignificance);


    myTree->Branch("daughters_jetNDauChargedMVASel",&_daughters_jetNDauChargedMVASel);
  myTree->Branch("daughters_miniRelIsoCharged",&_daughters_miniRelIsoCharged);
  myTree->Branch("daughters_miniRelIsoNeutral",&_daughters_miniRelIsoNeutral);
  myTree->Branch("daughters_jetPtRel",&_daughters_jetPtRel);
  myTree->Branch("daughters_jetPtRatio",&_daughters_jetPtRatio);
  myTree->Branch("daughters_jetBTagCSV",&_daughters_jetBTagCSV);
  myTree->Branch("daughters_jetBTagDeepCSV",&_daughters_jetBTagDeepCSV);
  myTree->Branch("daughters_jetBTagDeepFlavor",&_daughters_jetBTagDeepFlavor);
  //  myTree->Branch("daughters_lepMVA_mvaId",&_daughters_lepMVA_mvaId);*/

  /*if(doCPVariables){
    myTree->Branch("daughters_pca_x",&_daughters_pca_x);
    myTree->Branch("daughters_pca_y",&_daughters_pca_y);
    myTree->Branch("daughters_pca_z",&_daughters_pca_z);
    myTree->Branch("daughters_pcaRefitPV_x",&_daughters_pcaRefitPV_x);
    myTree->Branch("daughters_pcaRefitPV_y",&_daughters_pcaRefitPV_y);
    myTree->Branch("daughters_pcaRefitPV_z",&_daughters_pcaRefitPV_z);
    myTree->Branch("daughters_pcaRefitPVBS_x",&_daughters_pcaRefitPVBS_x);
    myTree->Branch("daughters_pcaRefitPVBS_y",&_daughters_pcaRefitPVBS_y);
    myTree->Branch("daughters_pcaRefitPVBS_z",&_daughters_pcaRefitPVBS_z);
    myTree->Branch("daughters_pcaGenPV_x",&_daughters_pcaGenPV_x);
    myTree->Branch("daughters_pcaGenPV_y",&_daughters_pcaGenPV_y);
    myTree->Branch("daughters_pcaGenPV_z",&_daughters_pcaGenPV_z);
  }

  myTree->Branch("PFTau_a1_lvp", &PFTau_a1_lvp);
  myTree->Branch("PFTau_a1_cov", &PFTau_a1_cov);
  myTree->Branch("PFTau_a1_charge", &PFTau_a1_charge);
  myTree->Branch("PFTau_a1_pdgid", &PFTau_a1_pdgid);
  myTree->Branch("PFTau_a1_B", &PFTau_a1_B);
  myTree->Branch("PFTau_a1_M", &PFTau_a1_M);
  myTree->Branch("PFTauTrack_deltaR", &PFTauTrack_deltaR);
  myTree->Branch("PFTauLeadTrackLV", &PFTauLeadTrackLV);*/

  //myTree->Branch("JetsNumber",&_numberOfJets,"JetsNumber/I");
  //myTree->Branch("jets_VBFleadFilterMatch"   , &_jets_VBFleadFilterMatch   ); //FRA
  //myTree->Branch("jets_VBFsubleadFilterMatch", &_jets_VBFsubleadFilterMatch); //FRA
  //myTree->Branch("jets_px",&_jets_px);
  //myTree->Branch("jets_py",&_jets_py);
  //myTree->Branch("jets_pz",&_jets_pz);
  //myTree->Branch("jets_e",&_jets_e);
  //myTree->Branch("jets_rawPt", &_jets_rawPt);
  /*myTree->Branch("jets_area", &_jets_area);
  myTree->Branch("jets_mT", &_jets_mT);
  myTree->Branch("jets_Flavour",&_jets_Flavour);
  myTree->Branch("jets_HadronFlavour",&_jets_HadronFlavour);
  myTree->Branch("jets_genjetIndex", &_jets_genjetIndex);
  myTree->Branch("jets_vtxPt", &_jets_vtxPt);
  myTree->Branch("jets_vtxMass", &_jets_vtxMass);
  myTree->Branch("jets_vtx3dL", &_jets_vtx3dL);
  myTree->Branch("jets_vtxNtrk", &_jets_vtxNtrk);
  myTree->Branch("jets_vtx3deL", &_jets_vtx3deL);
  myTree->Branch("jets_leadTrackPt", &_jets_leadTrackPt);
  myTree->Branch("jets_leptonPtRel", &_jets_leptonPtRel);
  myTree->Branch("jets_leptonPt", &_jets_leptonPt);
  myTree->Branch("jets_leptonDeltaR", &_jets_leptonDeltaR);
  myTree->Branch("jets_chEmEF" , &_jets_chEmEF);
  myTree->Branch("jets_chHEF"  , &_jets_chHEF);
  myTree->Branch("jets_nEmEF"  , &_jets_nEmEF);
  myTree->Branch("jets_nHEF"   , &_jets_nHEF);
  myTree->Branch("jets_MUF"    , &_jets_MUF);
  myTree->Branch("jets_neMult" , &_jets_neMult);
  myTree->Branch("jets_chMult" , &_jets_chMult);
  myTree->Branch("jets_jecUnc" , &_jets_jecUnc);
  myTree->Branch("pileupjetidMVA" , &_pileupjetidMVA);

  // JEC Regrouped uncertainty sources
  myTree->Branch("jets_jetUncRegrouped_FlavorQCD_up"  , &_SourceUncValRegrouped_up["FlavorQCD"]); // up variations
  myTree->Branch("jets_jetUncRegrouped_RelativeBal_up", &_SourceUncValRegrouped_up["RelativeBal"]);
  myTree->Branch("jets_jetUncRegrouped_HF_up"         , &_SourceUncValRegrouped_up["HF"]);
  myTree->Branch("jets_jetUncRegrouped_BBEC1_up"      , &_SourceUncValRegrouped_up["BBEC1"]);
  myTree->Branch("jets_jetUncRegrouped_EC2_up"        , &_SourceUncValRegrouped_up["EC2"]);
  myTree->Branch("jets_jetUncRegrouped_Absolute_up"   , &_SourceUncValRegrouped_up["Absolute"]);
  myTree->Branch("jets_jetUncRegrouped_Total_up"      , &_SourceUncValRegrouped_up["Total"]);
  myTree->Branch("jets_jetUncRegrouped_FlavorQCD_dw"  , &_SourceUncValRegrouped_dw["FlavorQCD"]); // down variations
  myTree->Branch("jets_jetUncRegrouped_RelativeBal_dw", &_SourceUncValRegrouped_dw["RelativeBal"]);
  myTree->Branch("jets_jetUncRegrouped_HF_dw"         , &_SourceUncValRegrouped_dw["HF"]);
  myTree->Branch("jets_jetUncRegrouped_BBEC1_dw"      , &_SourceUncValRegrouped_dw["BBEC1"]);
  myTree->Branch("jets_jetUncRegrouped_EC2_dw"        , &_SourceUncValRegrouped_dw["EC2"]);
  myTree->Branch("jets_jetUncRegrouped_Absolute_dw"   , &_SourceUncValRegrouped_dw["Absolute"]);
  myTree->Branch("jets_jetUncRegrouped_Total_dw"      , &_SourceUncValRegrouped_dw["Total"]);
  if (theYear == 2016)
    {
      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_up"         , &_SourceUncValRegrouped_up["BBEC1_2016"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_up"           , &_SourceUncValRegrouped_up["EC2_2016"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_up"      , &_SourceUncValRegrouped_up["Absolute_2016"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_up"            , &_SourceUncValRegrouped_up["HF_2016"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_up", &_SourceUncValRegrouped_up["RelativeSample_2016"]);

      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_dw"         , &_SourceUncValRegrouped_dw["BBEC1_2016"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_dw"           , &_SourceUncValRegrouped_dw["EC2_2016"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_dw"      , &_SourceUncValRegrouped_dw["Absolute_2016"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_dw"            , &_SourceUncValRegrouped_dw["HF_2016"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_dw", &_SourceUncValRegrouped_dw["RelativeSample_2016"]);
    }
  if (theYear == 2017)
    {
      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_up"         , &_SourceUncValRegrouped_up["BBEC1_2017"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_up"           , &_SourceUncValRegrouped_up["EC2_2017"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_up"      , &_SourceUncValRegrouped_up["Absolute_2017"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_up"            , &_SourceUncValRegrouped_up["HF_2017"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_up", &_SourceUncValRegrouped_up["RelativeSample_2017"]);

      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_dw"         , &_SourceUncValRegrouped_dw["BBEC1_2017"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_dw"           , &_SourceUncValRegrouped_dw["EC2_2017"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_dw"      , &_SourceUncValRegrouped_dw["Absolute_2017"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_dw"            , &_SourceUncValRegrouped_dw["HF_2017"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_dw", &_SourceUncValRegrouped_dw["RelativeSample_2017"]);
    }
  if (theYear == 2018)
    {
      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_up"         , &_SourceUncValRegrouped_up["BBEC1_2018"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_up"           , &_SourceUncValRegrouped_up["EC2_2018"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_up"      , &_SourceUncValRegrouped_up["Absolute_2018"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_up"            , &_SourceUncValRegrouped_up["HF_2018"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_up", &_SourceUncValRegrouped_up["RelativeSample_2018"]);

      myTree->Branch("jets_jetUncRegrouped_BBEC1_YEAR_dw"         , &_SourceUncValRegrouped_dw["BBEC1_2018"]); // up variations
      myTree->Branch("jets_jetUncRegrouped_EC2_YEAR_dw"           , &_SourceUncValRegrouped_dw["EC2_2018"]);
      myTree->Branch("jets_jetUncRegrouped_Absolute_YEAR_dw"      , &_SourceUncValRegrouped_dw["Absolute_2018"]);
      myTree->Branch("jets_jetUncRegrouped_HF_YEAR_dw"            , &_SourceUncValRegrouped_dw["HF_2018"]);
      myTree->Branch("jets_jetUncRegrouped_RelativeSample_YEAR_dw", &_SourceUncValRegrouped_dw["RelativeSample_2018"]);
    }

  myTree->Branch("jetsDown_px",&_jetsDown_px);
  myTree->Branch("jetsDown_py",&_jetsDown_py);
  myTree->Branch("jetsDown_pz",&_jetsDown_pz);
  myTree->Branch("jetsDown_e",&_jetsDown_e);
  myTree->Branch("jetsDown_area", &_jetsDown_area);
  myTree->Branch("jetsDown_mT", &_jetsDown_mT);
  myTree->Branch("jetsDown_Flavour",&_jetsDown_Flavour);
  myTree->Branch("jetsDown_HadronFlavour",&_jetsDown_HadronFlavour);
  myTree->Branch("jetsDown_genjetIndex", &_jetsDown_genjetIndex);
  myTree->Branch("jetsDown_leadTrackPt", &_jetsDown_leadTrackPt);
  myTree->Branch("jetsDown_leptonPtRel", &_jetsDown_leptonPtRel);
  myTree->Branch("jetsDown_leptonPt", &_jetsDown_leptonPt);
  myTree->Branch("jetsDown_leptonDeltaR", &_jetsDown_leptonDeltaR);
  myTree->Branch("jetsDown_chEmEF" , &_jetsDown_chEmEF);
  myTree->Branch("jetsDown_chHEF"  , &_jetsDown_chHEF);
  myTree->Branch("jetsDown_nEmEF"  , &_jetsDown_nEmEF);
  myTree->Branch("jetsDown_nHEF"   , &_jetsDown_nHEF);
  myTree->Branch("jetsDown_MUF"    , &_jetsDown_MUF);
  myTree->Branch("jetsDown_neMult" , &_jetsDown_neMult);
  myTree->Branch("jetsDown_chMult" , &_jetsDown_chMult);
  myTree->Branch("pileupjetidMVADown" , &_pileupjetidMVADown);

  myTree->Branch("jetsUp_px",&_jetsUp_px);
  myTree->Branch("jetsUp_py",&_jetsUp_py);
  myTree->Branch("jetsUp_pz",&_jetsUp_pz);
  myTree->Branch("jetsUp_e",&_jetsUp_e);
  myTree->Branch("jetsUp_area", &_jetsUp_area);
  myTree->Branch("jetsUp_mT", &_jetsUp_mT);
  myTree->Branch("jetsUp_Flavour",&_jetsUp_Flavour);
  myTree->Branch("jetsUp_HadronFlavour",&_jetsUp_HadronFlavour);
  myTree->Branch("jetsUp_genjetIndex", &_jetsUp_genjetIndex);
  myTree->Branch("jetsUp_leadTrackPt", &_jetsUp_leadTrackPt);
  myTree->Branch("jetsUp_leptonPtRel", &_jetsUp_leptonPtRel);
  myTree->Branch("jetsUp_leptonPt", &_jetsUp_leptonPt);
  myTree->Branch("jetsUp_leptonDeltaR", &_jetsUp_leptonDeltaR);
  myTree->Branch("jetsUp_chEmEF" , &_jetsUp_chEmEF);
  myTree->Branch("jetsUp_chHEF"  , &_jetsUp_chHEF);
  myTree->Branch("jetsUp_nEmEF"  , &_jetsUp_nEmEF);
  myTree->Branch("jetsUp_nHEF"   , &_jetsUp_nHEF);
  myTree->Branch("jetsUp_MUF"    , &_jetsUp_MUF);
  myTree->Branch("jetsUp_neMult" , &_jetsUp_neMult);
  myTree->Branch("jetsUp_chMult" , &_jetsUp_chMult);
  myTree->Branch("pileupjetidMVAUp" , &_pileupjetidMVAUp);*/

  // JEC Uncertainty sources
  /*myTree->Branch("jets_jetUnc_AbsoluteFlavMap_up"  , &_SourceUncVal_up["AbsoluteFlavMap"]); // up variations
  myTree->Branch("jets_jetUnc_AbsoluteMPFBias_up"  , &_SourceUncVal_up["AbsoluteMPFBias"]);
  myTree->Branch("jets_jetUnc_AbsoluteSample_up"   , &_SourceUncVal_up["AbsoluteSample"]);
  myTree->Branch("jets_jetUnc_AbsoluteScale_up"    , &_SourceUncVal_up["AbsoluteScale"]);
  myTree->Branch("jets_jetUnc_AbsoluteStat_up"     , &_SourceUncVal_up["AbsoluteStat"]);
  myTree->Branch("jets_jetUnc_FlavorQCD_up"        , &_SourceUncVal_up["FlavorQCD"]);
  myTree->Branch("jets_jetUnc_Fragmentation_up"    , &_SourceUncVal_up["Fragmentation"]);
  myTree->Branch("jets_jetUnc_PileUpDataMC_up"     , &_SourceUncVal_up["PileUpDataMC"]);
  myTree->Branch("jets_jetUnc_PileUpPtBB_up"       , &_SourceUncVal_up["PileUpPtBB"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC1_up"      , &_SourceUncVal_up["PileUpPtEC1"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC2_up"      , &_SourceUncVal_up["PileUpPtEC2"]);
  myTree->Branch("jets_jetUnc_PileUpPtHF_up"       , &_SourceUncVal_up["PileUpPtHF"]);
  myTree->Branch("jets_jetUnc_PileUpPtRef_up"      , &_SourceUncVal_up["PileUpPtRef"]);
  myTree->Branch("jets_jetUnc_RelativeBal_up"      , &_SourceUncVal_up["RelativeBal"]);
  myTree->Branch("jets_jetUnc_RelativeFSR_up"      , &_SourceUncVal_up["RelativeFSR"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC1_up"   , &_SourceUncVal_up["RelativeJEREC1"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC2_up"   , &_SourceUncVal_up["RelativeJEREC2"]);
  myTree->Branch("jets_jetUnc_RelativeJERHF_up"    , &_SourceUncVal_up["RelativeJERHF"]);
  myTree->Branch("jets_jetUnc_RelativePtBB_up"     , &_SourceUncVal_up["RelativePtBB"]);
  myTree->Branch("jets_jetUnc_RelativePtEC1_up"    , &_SourceUncVal_up["RelativePtEC1"]);
  myTree->Branch("jets_jetUnc_RelativePtEC2_up"    , &_SourceUncVal_up["RelativePtEC2"]);
  myTree->Branch("jets_jetUnc_RelativePtHF_up"     , &_SourceUncVal_up["RelativePtHF"]);
  myTree->Branch("jets_jetUnc_RelativeSample_up"   , &_SourceUncVal_up["RelativeSample"]);
  myTree->Branch("jets_jetUnc_RelativeStatEC_up"   , &_SourceUncVal_up["RelativeStatEC"]);
  myTree->Branch("jets_jetUnc_RelativeStatFSR_up"  , &_SourceUncVal_up["RelativeStatFSR"]);
  myTree->Branch("jets_jetUnc_RelativeStatHF_up"   , &_SourceUncVal_up["RelativeStatHF"]);
  myTree->Branch("jets_jetUnc_SinglePionECAL_up"   , &_SourceUncVal_up["SinglePionECAL"]);
  myTree->Branch("jets_jetUnc_SinglePionHCAL_up"   , &_SourceUncVal_up["SinglePionHCAL"]);
  myTree->Branch("jets_jetUnc_TimePtEta_up"        , &_SourceUncVal_up["TimePtEta"]);
  myTree->Branch("jets_jetUnc_AbsoluteFlavMap_dw"  , &_SourceUncVal_dw["AbsoluteFlavMap"]); // down variations
  myTree->Branch("jets_jetUnc_AbsoluteMPFBias_dw"  , &_SourceUncVal_dw["AbsoluteMPFBias"]);
  myTree->Branch("jets_jetUnc_AbsoluteSample_dw"   , &_SourceUncVal_dw["AbsoluteSample"]);
  myTree->Branch("jets_jetUnc_AbsoluteScale_dw"    , &_SourceUncVal_dw["AbsoluteScale"]);
  myTree->Branch("jets_jetUnc_AbsoluteStat_dw"     , &_SourceUncVal_dw["AbsoluteStat"]);
  myTree->Branch("jets_jetUnc_FlavorQCD_dw"        , &_SourceUncVal_dw["FlavorQCD"]);
  myTree->Branch("jets_jetUnc_Fragmentation_dw"    , &_SourceUncVal_dw["Fragmentation"]);
  myTree->Branch("jets_jetUnc_PileUpDataMC_dw"     , &_SourceUncVal_dw["PileUpDataMC"]);
  myTree->Branch("jets_jetUnc_PileUpPtBB_dw"       , &_SourceUncVal_dw["PileUpPtBB"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC1_dw"      , &_SourceUncVal_dw["PileUpPtEC1"]);
  myTree->Branch("jets_jetUnc_PileUpPtEC2_dw"      , &_SourceUncVal_dw["PileUpPtEC2"]);
  myTree->Branch("jets_jetUnc_PileUpPtHF_dw"       , &_SourceUncVal_dw["PileUpPtHF"]);
  myTree->Branch("jets_jetUnc_PileUpPtRef_dw"      , &_SourceUncVal_dw["PileUpPtRef"]);
  myTree->Branch("jets_jetUnc_RelativeBal_dw"      , &_SourceUncVal_dw["RelativeBal"]);
  myTree->Branch("jets_jetUnc_RelativeFSR_dw"      , &_SourceUncVal_dw["RelativeFSR"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC1_dw"   , &_SourceUncVal_dw["RelativeJEREC1"]);
  myTree->Branch("jets_jetUnc_RelativeJEREC2_dw"   , &_SourceUncVal_dw["RelativeJEREC2"]);
  myTree->Branch("jets_jetUnc_RelativeJERHF_dw"    , &_SourceUncVal_dw["RelativeJERHF"]);
  myTree->Branch("jets_jetUnc_RelativePtBB_dw"     , &_SourceUncVal_dw["RelativePtBB"]);
  myTree->Branch("jets_jetUnc_RelativePtEC1_dw"    , &_SourceUncVal_dw["RelativePtEC1"]);
  myTree->Branch("jets_jetUnc_RelativePtEC2_dw"    , &_SourceUncVal_dw["RelativePtEC2"]);
  myTree->Branch("jets_jetUnc_RelativePtHF_dw"     , &_SourceUncVal_dw["RelativePtHF"]);
  myTree->Branch("jets_jetUnc_RelativeSample_dw"   , &_SourceUncVal_dw["RelativeSample"]);
  myTree->Branch("jets_jetUnc_RelativeStatEC_dw"   , &_SourceUncVal_dw["RelativeStatEC"]);
  myTree->Branch("jets_jetUnc_RelativeStatFSR_dw"  , &_SourceUncVal_dw["RelativeStatFSR"]);
  myTree->Branch("jets_jetUnc_RelativeStatHF_dw"   , &_SourceUncVal_dw["RelativeStatHF"]);
  myTree->Branch("jets_jetUnc_SinglePionECAL_dw"   , &_SourceUncVal_dw["SinglePionECAL"]);
  myTree->Branch("jets_jetUnc_SinglePionHCAL_dw"   , &_SourceUncVal_dw["SinglePionHCAL"]);
  myTree->Branch("jets_jetUnc_TimePtEta_dw"        , &_SourceUncVal_dw["TimePtEta"]);

  myTree->Branch("bDiscriminator",&_bdiscr);
  myTree->Branch("bCSVscore",&_bdiscr2);
  myTree->Branch("pfCombinedMVAV2BJetTags",&_bdiscr3);
  myTree->Branch("bDeepCSV_probb",&_bdiscr4);
  myTree->Branch("bDeepCSV_probbb",&_bdiscr5);
  myTree->Branch("bDeepCSV_probudsg",&_bdiscr6);
  myTree->Branch("bDeepCSV_probc",&_bdiscr7);
  myTree->Branch("bDeepCSV_probcc",&_bdiscr8);
  myTree->Branch("bDeepFlavor_probb",&_bdiscr9);
  myTree->Branch("bDeepFlavor_probbb",&_bdiscr10);
  myTree->Branch("bDeepFlavor_problepb",&_bdiscr11);
  myTree->Branch("bDeepFlavor_probc",&_bdiscr12);
  myTree->Branch("bDeepFlavor_probuds",&_bdiscr13);
  myTree->Branch("bDeepFlavor_probg",&_bdiscr14);*/
  //myTree->Branch("PFjetID",&_jetID);
  /*myTree->Branch("looseJetID",&_looseJetID);
  myTree->Branch("tightJetID",&_tightJetID);
  myTree->Branch("tightLepVetoJetID",&_tightLepVetoJetID);

  myTree->Branch("jetRawf",&_jetrawf);
  myTree->Branch("jets_JER",&_jets_JER);

  myTree->Branch("bDiscriminatorDown",&_bdiscrDown);
  myTree->Branch("bCSVscoreDown",&_bdiscr2Down);
  myTree->Branch("pfCombinedMVAV2BJetTagsDown",&_bdiscr3Down);
  myTree->Branch("bDeepCSV_probbDown",&_bdiscr4Down);
  myTree->Branch("bDeepCSV_probbbDown",&_bdiscr5Down);
  myTree->Branch("bDeepCSV_probudsgDown",&_bdiscr6Down);
  myTree->Branch("bDeepCSV_probcDown",&_bdiscr7Down);
  myTree->Branch("bDeepCSV_probccDown",&_bdiscr8Down);
  myTree->Branch("bDeepFlavor_probbDown",&_bdiscr9Down);
  myTree->Branch("bDeepFlavor_probbbDown",&_bdiscr10Down);
  myTree->Branch("bDeepFlavor_problepbDown",&_bdiscr11Down);
  myTree->Branch("bDeepFlavor_probcDown",&_bdiscr12Down);
  myTree->Branch("bDeepFlavor_probudsDown",&_bdiscr13Down);
  myTree->Branch("bDeepFlavor_probgDown",&_bdiscr14Down);
  myTree->Branch("looseJetIDDown",&_looseJetIDDown);
  myTree->Branch("tightJetIDDown",&_tightJetIDDown);
  myTree->Branch("tightLepVetoJetIDDown",&_tightLepVetoJetIDDown);
  myTree->Branch("jetRawfDown",&_jetrawfDown);

  myTree->Branch("bDiscriminatorUp",&_bdiscrUp);
  myTree->Branch("bCSVscoreUp",&_bdiscr2Up);
  myTree->Branch("pfCombinedMVAV2BJetTagsUp",&_bdiscr3Up);
  myTree->Branch("bDeepCSV_probbUp",&_bdiscr4Up);
  myTree->Branch("bDeepCSV_probbbUp",&_bdiscr5Up);
  myTree->Branch("bDeepCSV_probudsgUp",&_bdiscr6Up);
  myTree->Branch("bDeepCSV_probcUp",&_bdiscr7Up);
  myTree->Branch("bDeepCSV_probccUp",&_bdiscr8Up);
  myTree->Branch("bDeepFlavor_probbUp",&_bdiscr9Up);
  myTree->Branch("bDeepFlavor_probbbUp",&_bdiscr10Up);
  myTree->Branch("bDeepFlavor_problepbUp",&_bdiscr11Up);
  myTree->Branch("bDeepFlavor_probcUp",&_bdiscr12Up);
  myTree->Branch("bDeepFlavor_probudsUp",&_bdiscr13Up);
  myTree->Branch("bDeepFlavor_probgUp",&_bdiscr14Up);
  myTree->Branch("looseJetIDUp",&_looseJetIDUp);
  myTree->Branch("tightJetIDUp",&_tightJetIDUp);
  myTree->Branch("tightLepVetoJetIDUp",&_tightLepVetoJetIDUp);
  myTree->Branch("jetRawfUp",&_jetrawfUp);

  myTree->Branch("susyModel",&_susyModel);
  if (computeQGVar) myTree->Branch("jets_QGdiscr" , &_jets_QGdiscr);
  myTree->Branch("ak8jets_px", &_ak8jets_px);
  myTree->Branch("ak8jets_py", &_ak8jets_py);
  myTree->Branch("ak8jets_pz", &_ak8jets_pz);
  myTree->Branch("ak8jets_e", &_ak8jets_e);
  myTree->Branch("ak8jets_SoftDropMass", &_ak8jets_SoftDropMass);
  myTree->Branch("ak8jets_PrunedMass", &_ak8jets_PrunedMass);
  myTree->Branch("ak8jets_TrimmedMass", &_ak8jets_TrimmedMass);
  myTree->Branch("ak8jets_FilteredMass", &_ak8jets_FilteredMass);
  myTree->Branch("ak8jets_tau1", &_ak8jets_tau1);
  myTree->Branch("ak8jets_tau2", &_ak8jets_tau2);
  myTree->Branch("ak8jets_tau3", &_ak8jets_tau3);
  myTree->Branch("ak8jets_tau4", &_ak8jets_tau4);
  myTree->Branch("ak8jets_CSV", &_ak8jets_CSV);
  myTree->Branch("ak8jets_deepCSV_probb", &_ak8jets_deepCSV_probb);
  myTree->Branch("ak8jets_deepCSV_probbb", &_ak8jets_deepCSV_probbb);
  myTree->Branch("ak8jets_deepFlavor_probb", &_ak8jets_deepFlavor_probb);
  myTree->Branch("ak8jets_deepFlavor_probbb", &_ak8jets_deepFlavor_probbb);
  myTree->Branch("ak8jets_deepFlavor_problepb", &_ak8jets_deepFlavor_problepb);
  myTree->Branch("ak8jets_nsubjets", &_ak8jets_nsubjets);

  myTree->Branch("subjets_px", &_subjets_px);
  myTree->Branch("subjets_py", &_subjets_py);
  myTree->Branch("subjets_pz", &_subjets_pz);
  myTree->Branch("subjets_e", &_subjets_e);
  myTree->Branch("subjets_CSV", &_subjets_CSV);
  myTree->Branch("subjets_deepCSV_probb", &_subjets_deepCSV_probb);
  myTree->Branch("subjets_deepCSV_probbb", &_subjets_deepCSV_probbb);
  myTree->Branch("subjets_deepFlavor_probb", &_subjets_deepFlavor_probb);
  myTree->Branch("subjets_deepFlavor_probbb", &_subjets_deepFlavor_probbb);
  myTree->Branch("subjets_deepFlavor_problepb", &_subjets_deepFlavor_problepb);
  myTree->Branch("subjets_ak8MotherIdx", &_subjets_ak8MotherIdx);


  myTree->Branch("pv_x", &_pv_x);
  myTree->Branch("pv_y", &_pv_y);
  myTree->Branch("pv_z", &_pv_z);
  myTree->Branch("pv_cov", &_pvRefit_cov);

  myTree->Branch("pvRefit_x", &_pvRefit_x);
  myTree->Branch("pvRefit_y", &_pvRefit_y);
  myTree->Branch("pvRefit_z", &_pvRefit_z);

  myTree->Branch("RefitPVBS_x", &_RefitPVBS_x);
  myTree->Branch("RefitPVBS_y", &_RefitPVBS_y);
  myTree->Branch("RefitPVBS_z", &_RefitPVBS_z);

  myTree->Branch("RefitPVNoBS_x", &_RefitPVNoBS_x);
  myTree->Branch("RefitPVNoBS_y", &_RefitPVNoBS_y);
  myTree->Branch("RefitPVNoBS_z", &_RefitPVNoBS_z);

  myTree->Branch("RefitPVBS_xError", &_RefitPVBS_xError);
  myTree->Branch("RefitPVBS_yError", &_RefitPVBS_yError);
  myTree->Branch("RefitPVBS_zError", &_RefitPVBS_zError);

  myTree->Branch("RefitPVBS_Cov", &_RefitPVBS_Cov);

  myTree->Branch("RefitPVNoBS_xError", &_RefitPVNoBS_xError);
  myTree->Branch("RefitPVNoBS_yError", &_RefitPVNoBS_yError);
  myTree->Branch("RefitPVNoBS_zError", &_RefitPVNoBS_zError);

  myTree->Branch("VertexHashBS1", &_VertexHashBS1);
  myTree->Branch("VertexHashBS2", &_VertexHashBS2);
  myTree->Branch("VertexHashNoBS1", &_VertexHashNoBS1);
  myTree->Branch("VertexHashNoBS2", &_VertexHashNoBS2);

  myTree->Branch("VertexHashBSTracksRemovedOld1", &_VertexHashBSTracksRemovedOld1);
  myTree->Branch("VertexHashBSTracksRemovedOld2", &_VertexHashBSTracksRemovedOld2);
  myTree->Branch("VertexHashNoBSTracksRemovedOld1", &_VertexHashNoBSTracksRemovedOld1);
  myTree->Branch("VertexHashNoBSTracksRemovedOld2", &_VertexHashNoBSTracksRemovedOld2);

  myTree->Branch("LeptonHash", &_LeptonHash);

  myTree->Branch("pvGen_x", &_pvGen_x);
  myTree->Branch("pvGen_y", &_pvGen_y);
  myTree->Branch("pvGen_z", &_pvGen_z);

  myTree->Branch("isRefitPV", &_isRefitPV);*/

  //Synchro

  //myTree->Branch("run",&_runNumber,"RunNumber/I");
  //myTree->Branch("lumi",&_lumi,"lumi/I");
  //myTree->Branch("evt",&_indexevents,"EventNumber/l");
  //myTree->Branch("SelectedPairs",&SelectedPairs,"SelectedPairs/I");
  //myTree->Branch("dilepton_veto", &dilepton_veto);
  /*myTree->Branch("extraelec_veto", &eleveto);
  myTree->Branch("extramuon_veto", &muonveto);

  myTree->Branch("trg_singleelectron", &trg_singleelectron);
  myTree->Branch("trg_singlemuon", &trg_singlemuon);
  myTree->Branch("trg_singletau", &trg_singletau);
  myTree->Branch("trg_muonelectron", &trg_muonelectron);
  myTree->Branch("trg_mutaucross", &trg_mutaucross);
  myTree->Branch("trg_doubletau", &trg_doubletau);

  myTree->Branch("pt_1", &pt_1);
  myTree->Branch("phi_1", &phi_1);
  myTree->Branch("eta_1", &eta_1);
  myTree->Branch("m_1", &m_1);
  myTree->Branch("q_1", &q_1);
  myTree->Branch("d0_1", &d0_1);
  myTree->Branch("dz_1", &dz_1);
  myTree->Branch("mt_1", &mt_1);*/
  //myTree->Branch("pfmt_1", &pfmt_1);
  /*myTree->Branch("puppimt_1", &puppimt_1);
  myTree->Branch("iso_1", &iso_1);
  myTree->Branch("gen_match_1", &gen_match_1);
  myTree->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1);
  myTree->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1);
  myTree->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
  myTree->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
  myTree->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1);
  myTree->Branch("againstMuonLoose3_1", &againstMuonLoose3_1);
  myTree->Branch("againstMuonTight3_1", &againstMuonTight3_1);
  myTree->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1);
  myTree->Branch("idisoweight_1", &idisoweight_1);
  myTree->Branch("antieweight_1", &antieweight_1);
  myTree->Branch("antimuweight_1", &antimuweight_1);

  myTree->Branch("pt_2", &pt_2);
  myTree->Branch("phi_2", &phi_2);
  myTree->Branch("eta_2", &eta_2);
  myTree->Branch("m_2", &m_2);
  myTree->Branch("q_2", &q_2);
  myTree->Branch("d0_2 ", &d0_2);
  myTree->Branch("dz_2", &dz_2);
  myTree->Branch("mt_2", &mt_2);*/
  //myTree->Branch("pfmt_2", &pfmt_2);
  /*myTree->Branch("puppimt_2", &puppimt_2);
  myTree->Branch("iso_2", &iso_2);
  myTree->Branch("gen_match_2", &gen_match_2);
  myTree->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2);
  myTree->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2);
  myTree->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
  myTree->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
  myTree->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2);
  myTree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
  myTree->Branch("againstMuonTight3_2", &againstMuonTight3_2);
  myTree->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
  myTree->Branch("idisoweight_2", &idisoweight_2);
  myTree->Branch("antieweight_2", &antieweight_2);
  myTree->Branch("antimuweight_2", &antimuweight_2);

  myTree->Branch("pt_tt", &pt_tt);
  myTree->Branch("pt_vis", &pt_vis);
  myTree->Branch("mt_tot", &mt_tot);
  myTree->Branch("m_vis", &m_vis);*/
  //myTree->Branch("m_sv", &m_sv);
  //myTree->Branch("mt_sv", &mt_sv);

  //myTree->Branch("met",&_met,"met/F");
  //myTree->Branch("metphi",&_metphi,"metphi/F");
  //myTree->Branch("PUPPImet",&_PUPPImet,"PUPPImet/F");
  //myTree->Branch("PUPPImetphi",&_PUPPImetphi,"PUPPImetphi/F");
  //myTree->Branch("pzetavis", &pzetavis);
  //myTree->Branch("pzetamiss", &pzetamiss);
  //myTree->Branch("pfpzetamiss", &pfpzetamiss);
 // myTree->Branch("puppipzetamiss", &puppipzetamiss);
  // myTree->Branch("metcov00",&_PFMETCov00,"metcov00/F");
  // myTree->Branch("metcov01",&_PFMETCov01,"metcov01/F");
  // myTree->Branch("metcov10",&_PFMETCov10,"metcov10/F");
  // myTree->Branch("metcov11",&_PFMETCov11,"metcov11/F");
 /* myTree->Branch("puppimetcov00",&_PUPPIMETCov00,"puppimetcov00/F");
  myTree->Branch("puppimetcov01",&_PUPPIMETCov01,"puppimetcov01/F");
  myTree->Branch("puppimetcov10",&_PUPPIMETCov10,"puppimetcov10/F");
  myTree->Branch("puppimetcov11",&_PUPPIMETCov11,"puppimetcov11/F");
  myTree->Branch("puppimet_ex_JetEnUp",&_puppimet_ex_JetEnUp,"puppimet_ex_JetEnUp/F");
  myTree->Branch("puppimet_ey_JetEnUp",&_puppimet_ey_JetEnUp,"puppimet_ey_JetEnUp/F");
  myTree->Branch("puppimet_ex_JetEnDown",&_puppimet_ex_JetEnDown,"puppimet_ex_JetEnDown/F");
  myTree->Branch("puppimet_ey_JetEnDown",&_puppimet_ey_JetEnDown,"puppimet_ey_JetEnDown/F");
  myTree->Branch("puppimet_ex_UnclusteredEnUp",&_puppimet_ex_UnclusteredEnUp,"puppimet_ex_UnclusteredEnUp/F");
  myTree->Branch("puppimet_ey_UnclusteredEnUp",&_puppimet_ey_UnclusteredEnUp,"puppimet_ey_UnclusteredEnUp/F");
  myTree->Branch("puppimet_ex_UnclusteredEnDown",&_puppimet_ex_UnclusteredEnDown,"puppimet_ex_UnclusteredEnDown/F");
  myTree->Branch("puppimet_ey_UnclusteredEnDown",&_puppimet_ey_UnclusteredEnDown,"puppimet_ey_UnclusteredEnDown/F");
  myTree->Branch("puppimet_ex_JetResUp",&_puppimet_ex_JetResUp,"puppimet_ex_JetResUp/F");
  myTree->Branch("puppimet_ey_JetResUp",&_puppimet_ey_JetResUp,"puppimet_ey_JetResUp/F");
  myTree->Branch("puppimet_ex_JetResDown",&_puppimet_ex_JetResDown,"puppimet_ex_JetResDown/F");
  myTree->Branch("puppimet_ey_JetResDown",&_puppimet_ey_JetResDown,"puppimet_ey_JetResDown/F");

  myTree->Branch("mjj", &mjj);
  myTree->Branch("jdeta", &jdeta);
  myTree->Branch("njetingap", &njetingap);
  myTree->Branch("njetingap20", &njetingap20);
  myTree->Branch("jdphi", &jdphi);
  myTree->Branch("dijetpt", &dijetpt);
  myTree->Branch("dijetphi", &dijetphi);
  myTree->Branch("ptvis", &ptvis);

  myTree->Branch("nbtag",&nbtag);
  myTree->Branch("njets", &njets);
  myTree->Branch("njetspt20", &njetspt20);
  myTree->Branch("jpt_1", &jpt_1);
  myTree->Branch("jeta_1", &jeta_1);
  myTree->Branch("jphi_1", &jphi_1);
  myTree->Branch("jcsv_1", &jcsv_1);
  myTree->Branch("jpt_2", &jpt_2);
  myTree->Branch("jeta_2", &jeta_2);
  myTree->Branch("jphi_2", &jphi_2);
  myTree->Branch("jcsv_2", &jcsv_2);
  myTree->Branch("bpt_1", &bpt_1);
  myTree->Branch("beta_1", &beta_1);
  myTree->Branch("bphi_1", &bphi_1);
  myTree->Branch("bcsv_1", &bcsv_1);
  myTree->Branch("bpt_2", &bpt_2);
  myTree->Branch("beta_2", &beta_2);
  myTree->Branch("bphi_2", &bphi_2);
  myTree->Branch("bcsv_2 ", &bcsv_2);

  myTree->Branch("mjjDown", &mjjDown);
  myTree->Branch("jdetaDown", &jdetaDown);
  myTree->Branch("njetingapDown", &njetingapDown);
  myTree->Branch("njetingap20Down", &njetingap20Down);
  myTree->Branch("jdphiDown", &jdphiDown);
  myTree->Branch("dijetptDown", &dijetptDown);
  myTree->Branch("dijetphiDown", &dijetphiDown);
  myTree->Branch("ptvisDown", &ptvisDown);

  myTree->Branch("nbtagDown",&nbtagDown);
  myTree->Branch("njetsDown", &njetsDown);
  myTree->Branch("njetspt20Down", &njetspt20Down);
  myTree->Branch("jptDown_1", &jptDown_1);
  myTree->Branch("jetaDown_1", &jetaDown_1);
  myTree->Branch("jphiDown_1", &jphiDown_1);
  myTree->Branch("jcsvDown_1", &jcsvDown_1);
  myTree->Branch("jptDown_2", &jptDown_2);
  myTree->Branch("jetaDown_2", &jetaDown_2);
  myTree->Branch("jphiDown_2", &jphiDown_2);
  myTree->Branch("jcsvDown_2", &jcsvDown_2);
  myTree->Branch("bptDown_1", &bptDown_1);
  myTree->Branch("betaDown_1", &betaDown_1);
  myTree->Branch("bphiDown_1", &bphiDown_1);
  myTree->Branch("bcsvDown_1", &bcsvDown_1);
  myTree->Branch("bptDown_2", &bptDown_2);
  myTree->Branch("betaDown_2", &betaDown_2);
  myTree->Branch("bphiDown_2", &bphiDown_2);
  myTree->Branch("bcsvDown_2 ", &bcsvDown_2);

  myTree->Branch("mjjUp", &mjjUp);
  myTree->Branch("jdetaUp", &jdetaUp);
  myTree->Branch("njetingapUp", &njetingapUp);
  myTree->Branch("njetingap20Up", &njetingap20Up);
  myTree->Branch("jdphiUp", &jdphiUp);
  myTree->Branch("dijetptUp", &dijetptUp);
  myTree->Branch("dijetphiUp", &dijetphiUp);
  myTree->Branch("ptvisUp", &ptvisUp);

  myTree->Branch("nbtagUp",&nbtagUp);
  myTree->Branch("njetsUp", &njetsUp);
  myTree->Branch("njetspt20Up", &njetspt20Up);
  myTree->Branch("jptUp_1", &jptUp_1);
  myTree->Branch("jetaUp_1", &jetaUp_1);
  myTree->Branch("jphiUp_1", &jphiUp_1);
  myTree->Branch("jcsvUp_1", &jcsvUp_1);
  myTree->Branch("jptUp_2", &jptUp_2);
  myTree->Branch("jetaUp_2", &jetaUp_2);
  myTree->Branch("jphiUp_2", &jphiUp_2);
  myTree->Branch("jcsvUp_2", &jcsvUp_2);
  myTree->Branch("bptUp_1", &bptUp_1);
  myTree->Branch("betaUp_1", &betaUp_1);
  myTree->Branch("bphiUp_1", &bphiUp_1);
  myTree->Branch("bcsvUp_1", &bcsvUp_1);
  myTree->Branch("bptUp_2", &bptUp_2);
  myTree->Branch("betaUp_2", &betaUp_2);
  myTree->Branch("bphiUp_2", &bphiUp_2);
  myTree->Branch("bcsvUp_2 ", &bcsvUp_2);

  myTree->Branch("puweight", &puweight);

  myTree->Branch("NUP", &_nup,"NUP/I");

  myTree->Branch("weight", &weight);

  myTree->Branch("jpfid_1", &jpfid_1);
  myTree->Branch("jpuid_1", &jpuid_1);
  myTree->Branch("jpfid_2", &jpfid_2);
  myTree->Branch("jpuid_2", &jpuid_2);
  myTree->Branch("bpfid_1", &bpfid_1);
  myTree->Branch("bpuid_1", &bpuid_1);
  myTree->Branch("bpfid_2", &bpfid_2);
  myTree->Branch("bpuid_2", &bpuid_2);

  myTree->Branch("npv",&_npv,"npv/I");
  myTree->Branch("npu",&_npu,"npu/F");
  myTree->Branch("rho",&_rho,"rho/F");

  myTree->Branch("pt_sv", &pt_sv);
  myTree->Branch("eta_sv", &eta_sv);
  myTree->Branch("phi_sv", &phi_sv);
  myTree->Branch("met_sv", &met_sv);


  myTree->Branch("deepTauVsJetRaw_1", &byDeepTau2017v2p1VSjetraw_1);
  myTree->Branch("byVVVLooseDeepTau2017v2p1VSjet_1", &byVVVLooseDeepTau2017v2p1VSjet_1);
  myTree->Branch("byVVLooseDeepTau2017v2p1VSjet_1", &byVVLooseDeepTau2017v2p1VSjet_1);
  myTree->Branch("byVLooseDeepTau2017v2p1VSjet_1", &byVLooseDeepTau2017v2p1VSjet_1);
  myTree->Branch("byLooseDeepTau2017v2p1VSjet_1", &byLooseDeepTau2017v2p1VSjet_1);
  myTree->Branch("byMediumDeepTau2017v2p1VSjet_1", &byMediumDeepTau2017v2p1VSjet_1);
  myTree->Branch("byTightDeepTau2017v2p1VSjet_1", &byTightDeepTau2017v2p1VSjet_1);
  myTree->Branch("byVTightDeepTau2017v2p1VSjet_1", &byVTightDeepTau2017v2p1VSjet_1);
  myTree->Branch("byVVTightDeepTau2017v2p1VSjet_1", &byVVTightDeepTau2017v2p1VSjet_1);

  myTree->Branch("deepTauVsEleRaw_1", &byDeepTau2017v2p1VSeleraw_1);
  myTree->Branch("byVVVLooseDeepTau2017v2p1VSe_1",&byVVVLooseDeepTau2017v2p1VSe_1);
  myTree->Branch("byVVLooseDeepTau2017v2p1VSe_1",&byVVLooseDeepTau2017v2p1VSe_1);
  myTree->Branch("byVLooseDeepTau2017v2p1VSe_1",&byVLooseDeepTau2017v2p1VSe_1);
  myTree->Branch("byLooseDeepTau2017v2p1VSe_1",&byLooseDeepTau2017v2p1VSe_1);
  myTree->Branch("byMediumDeepTau2017v2p1VSe_1",&byMediumDeepTau2017v2p1VSe_1);
  myTree->Branch("byTightDeepTau2017v2p1VSe_1",&byTightDeepTau2017v2p1VSe_1);
  myTree->Branch("byVTightDeepTau2017v2p1VSe_1",&byVTightDeepTau2017v2p1VSe_1);
  myTree->Branch("byVVTightDeepTau2017v2p1VSe_1",&byVVTightDeepTau2017v2p1VSe_1);

  myTree->Branch("deepTauVsMuRaw_1", &byDeepTau2017v2p1VSmuraw_1);
  myTree->Branch("byDeepTau2017v2p1VSmuraw_1",&byDeepTau2017v2p1VSmuraw_1);
  myTree->Branch("byVLooseDeepTau2017v2p1VSmu_1",&byVLooseDeepTau2017v2p1VSmu_1);
  myTree->Branch("byLooseDeepTau2017v2p1VSmu_1",&byLooseDeepTau2017v2p1VSmu_1);
  myTree->Branch("byMediumDeepTau2017v2p1VSmu_1",&byMediumDeepTau2017v2p1VSmu_1);
  myTree->Branch("byTightDeepTau2017v2p1VSmu_1",&byTightDeepTau2017v2p1VSmu_1);

  myTree->Branch("deepTauVsJetRaw_2", &byDeepTau2017v2p1VSjetraw_2);
  myTree->Branch("byVVVLooseDeepTau2017v2p1VSjet_2", &byVVVLooseDeepTau2017v2p1VSjet_2);
  myTree->Branch("byVVLooseDeepTau2017v2p1VSjet_2", &byVVLooseDeepTau2017v2p1VSjet_2);
  myTree->Branch("byVLooseDeepTau2017v2p1VSjet_2", &byVLooseDeepTau2017v2p1VSjet_2);
  myTree->Branch("byLooseDeepTau2017v2p1VSjet_2", &byLooseDeepTau2017v2p1VSjet_2);
  myTree->Branch("byMediumDeepTau2017v2p1VSjet_2", &byMediumDeepTau2017v2p1VSjet_2);
  myTree->Branch("byTightDeepTau2017v2p1VSjet_2", &byTightDeepTau2017v2p1VSjet_2);
  myTree->Branch("byVTightDeepTau2017v2p1VSjet_2", &byVTightDeepTau2017v2p1VSjet_2);
  myTree->Branch("byVVTightDeepTau2017v2p1VSjet_2", &byVVTightDeepTau2017v2p1VSjet_2);

  myTree->Branch("deepTauVsEleRaw_2", &byDeepTau2017v2p1VSeleraw_2);
  myTree->Branch("byVVVLooseDeepTau2017v2p1VSe_2",&byVVVLooseDeepTau2017v2p1VSe_2);
  myTree->Branch("byVVLooseDeepTau2017v2p1VSe_2",&byVVLooseDeepTau2017v2p1VSe_2);
  myTree->Branch("byVLooseDeepTau2017v2p1VSe_2",&byVLooseDeepTau2017v2p1VSe_2);
  myTree->Branch("byLooseDeepTau2017v2p1VSe_2",&byLooseDeepTau2017v2p1VSe_2);
  myTree->Branch("byMediumDeepTau2017v2p1VSe_2",&byMediumDeepTau2017v2p1VSe_2);
  myTree->Branch("byTightDeepTau2017v2p1VSe_2",&byTightDeepTau2017v2p1VSe_2);
  myTree->Branch("byVTightDeepTau2017v2p1VSe_2",&byVTightDeepTau2017v2p1VSe_2);
  myTree->Branch("byVVTightDeepTau2017v2p1VSe_2",&byVVTightDeepTau2017v2p1VSe_2);

  myTree->Branch("deepTauVsMuRaw_2", &byDeepTau2017v2p1VSmuraw_2);
  myTree->Branch("byVLooseDeepTau2017v2p1VSmu_2",&byVLooseDeepTau2017v2p1VSmu_2);
  myTree->Branch("byLooseDeepTau2017v2p1VSmu_2",&byLooseDeepTau2017v2p1VSmu_2);
  myTree->Branch("byMediumDeepTau2017v2p1VSmu_2",&byMediumDeepTau2017v2p1VSmu_2);
  myTree->Branch("byTightDeepTau2017v2p1VSmu_2",&byTightDeepTau2017v2p1VSmu_2);

  myTree->Branch("pvx", &pvx);
  myTree->Branch("pvy ", &pvy);
  myTree->Branch("pvz", &pvz);
  myTree->Branch("pvcov", &pvcov);
  myTree->Branch("pcax", &pcax);
  myTree->Branch("pcay", &pcay);
  myTree->Branch("pcaz", &pcaz);

  myTree->Branch("dm_1", &dm_1);
  myTree->Branch("dmMVA_1", &dmMVA_1);

  myTree->Branch("dm_2", &dm_2);
  myTree->Branch("dmMVA_2", &dmMVA_2);

  myTree->Branch("chpt_1", &_chpt_1);
  myTree->Branch("chphi_1", &_chphi_1);
  myTree->Branch("cheta_1", &_cheta_1);
  myTree->Branch("chm_1", &_chm_1);

  myTree->Branch("chpt_2", &_chpt_2);
  myTree->Branch("chphi_2", &_chphi_2);
  myTree->Branch("cheta_2", &_cheta_2);
  myTree->Branch("chm_2", &_chm_2);

  myTree->Branch("npt_1", &_npt_1);
  myTree->Branch("nphi_1", &_nphi_1);
  myTree->Branch("neta_1", &_neta_1);
  myTree->Branch("nm_1", &_nm_1);

  myTree->Branch("npt_2", &_npt_2);
  myTree->Branch("nphi_2", &_nphi_2);
  myTree->Branch("neta_2", &_neta_2);
  myTree->Branch("nm_2", &_nm_2);

  myTree->Branch("svx_1", &svx_1);
  myTree->Branch("svy_1", &svy_1);
  myTree->Branch("svz_1", &svz_1);

  myTree->Branch("svx_2", &svx_2);
  myTree->Branch("svy_2", &svy_2);
  myTree->Branch("svz_2", &svz_2);

  myTree->Branch("tau1IndexVect",&tau1IndexVect);
  myTree->Branch("tau2IndexVect",&tau2IndexVect);*/

  //myTree->Branch("tauspinnerH", &tauspinnerH);
  //myTree->Branch("tauspinnerA", &tauspinnerA);
  //myTree->Branch("tauspinnerMaxMix", &tauspinnerMaxMix);
}

Int_t HTauTauNtuplizer::FindCandIndex(const reco::Candidate& cand,Int_t iCand=0){
  const reco::Candidate *daughter = cand.daughter(iCand);
  for(UInt_t iLeptons=0;iLeptons<_softLeptons.size();iLeptons++){
    if(daughter->masterClone().get()==_softLeptons.at(iLeptons)){
      return iLeptons;
    }
  }
  return -1;
}
// ----Analyzer (main) ----
// ------------ method called for each event  ------------
void HTauTauNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{

  //  if(DEBUG)printf("\n\n\n===New Event===   ");
  //std::cout<<" =========  new event =======" << std::endl;
  Initialize();
  findPrimaryVertices(event, eSetup);
  Handle<vector<reco::Vertex> >  vertexs;
  //event.getByLabel("offlineSlimmedPrimaryVertices",vertex);
  event.getByToken(theVtxTag,vertexs);

  // JEC Reduced fill map with "sourceName, var"
  // Up variations
  _SourceUncValRegrouped_up.emplace("FlavorQCD", _jets_jetUncRegrouped_FlavorQCD_up);
  _SourceUncValRegrouped_up.emplace("RelativeBal", _jets_jetUncRegrouped_RelativeBal_up);
  _SourceUncValRegrouped_up.emplace("HF", _jets_jetUncRegrouped_HF_up);
  _SourceUncValRegrouped_up.emplace("BBEC1", _jets_jetUncRegrouped_BBEC1_up);
  _SourceUncValRegrouped_up.emplace("EC2", _jets_jetUncRegrouped_EC2_up);
  _SourceUncValRegrouped_up.emplace("Absolute", _jets_jetUncRegrouped_Absolute_up);
  _SourceUncValRegrouped_up.emplace("Total", _jets_jetUncRegrouped_Total_up);
  // Down variations
  _SourceUncValRegrouped_dw.emplace("FlavorQCD", _jets_jetUncRegrouped_FlavorQCD_dw);
  _SourceUncValRegrouped_dw.emplace("RelativeBal", _jets_jetUncRegrouped_RelativeBal_dw);
  _SourceUncValRegrouped_dw.emplace("HF", _jets_jetUncRegrouped_HF_dw);
  _SourceUncValRegrouped_dw.emplace("BBEC1", _jets_jetUncRegrouped_BBEC1_dw);
  _SourceUncValRegrouped_dw.emplace("EC2", _jets_jetUncRegrouped_EC2_dw);
  _SourceUncValRegrouped_dw.emplace("Absolute", _jets_jetUncRegrouped_Absolute_dw);
  _SourceUncValRegrouped_dw.emplace("Total", _jets_jetUncRegrouped_Total_dw);

  if (theYear == 2016)
    {
      _SourceUncValRegrouped_up.emplace("BBEC1_2016"         , _jets_jetUncRegrouped_BBEC1_YEAR_up);
      _SourceUncValRegrouped_up.emplace("EC2_2016"           , _jets_jetUncRegrouped_EC2_YEAR_up);
      _SourceUncValRegrouped_up.emplace("Absolute_2016"      , _jets_jetUncRegrouped_Absolute_YEAR_up);
      _SourceUncValRegrouped_up.emplace("HF_2016"            , _jets_jetUncRegrouped_HF_YEAR_up);
      _SourceUncValRegrouped_up.emplace("RelativeSample_2016", _jets_jetUncRegrouped_RelativeSample_YEAR_up);

      _SourceUncValRegrouped_dw.emplace("BBEC1_2016"         , _jets_jetUncRegrouped_BBEC1_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("EC2_2016"           , _jets_jetUncRegrouped_EC2_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("Absolute_2016"      , _jets_jetUncRegrouped_Absolute_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("HF_2016"            , _jets_jetUncRegrouped_HF_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("RelativeSample_2016", _jets_jetUncRegrouped_RelativeSample_YEAR_dw);
    }
  else if (theYear == 2017)
    {
      _SourceUncValRegrouped_up.emplace("BBEC1_2017"         , _jets_jetUncRegrouped_BBEC1_YEAR_up);
      _SourceUncValRegrouped_up.emplace("EC2_2017"           , _jets_jetUncRegrouped_EC2_YEAR_up);
      _SourceUncValRegrouped_up.emplace("Absolute_2017"      , _jets_jetUncRegrouped_Absolute_YEAR_up);
      _SourceUncValRegrouped_up.emplace("HF_2017"            , _jets_jetUncRegrouped_HF_YEAR_up);
      _SourceUncValRegrouped_up.emplace("RelativeSample_2017", _jets_jetUncRegrouped_RelativeSample_YEAR_up);

      _SourceUncValRegrouped_dw.emplace("BBEC1_2017"         , _jets_jetUncRegrouped_BBEC1_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("EC2_2017"           , _jets_jetUncRegrouped_EC2_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("Absolute_2017"      , _jets_jetUncRegrouped_Absolute_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("HF_2017"            , _jets_jetUncRegrouped_HF_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("RelativeSample_2017", _jets_jetUncRegrouped_RelativeSample_YEAR_dw);
    }
  else if (theYear == 2018)
    {
      _SourceUncValRegrouped_up.emplace("BBEC1_2018"         , _jets_jetUncRegrouped_BBEC1_YEAR_up);
      _SourceUncValRegrouped_up.emplace("EC2_2018"           , _jets_jetUncRegrouped_EC2_YEAR_up);
      _SourceUncValRegrouped_up.emplace("Absolute_2018"      , _jets_jetUncRegrouped_Absolute_YEAR_up);
      _SourceUncValRegrouped_up.emplace("HF_2018"            , _jets_jetUncRegrouped_HF_YEAR_up);
      _SourceUncValRegrouped_up.emplace("RelativeSample_2018", _jets_jetUncRegrouped_RelativeSample_YEAR_up);

      _SourceUncValRegrouped_dw.emplace("BBEC1_2018"         , _jets_jetUncRegrouped_BBEC1_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("EC2_2018"           , _jets_jetUncRegrouped_EC2_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("Absolute_2018"      , _jets_jetUncRegrouped_Absolute_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("HF_2018"            , _jets_jetUncRegrouped_HF_YEAR_dw);
      _SourceUncValRegrouped_dw.emplace("RelativeSample_2018", _jets_jetUncRegrouped_RelativeSample_YEAR_dw);
    }

  // JEC fill map with "sourceName, var"
  // Up variations
  _SourceUncVal_up.emplace("AbsoluteFlavMap" ,_jets_jetUnc_AbsoluteFlavMap_up);
  _SourceUncVal_up.emplace("AbsoluteMPFBias" ,_jets_jetUnc_AbsoluteMPFBias_up);
  _SourceUncVal_up.emplace("AbsoluteSample"  ,_jets_jetUnc_AbsoluteSample_up);
  _SourceUncVal_up.emplace("AbsoluteScale"   ,_jets_jetUnc_AbsoluteScale_up);
  _SourceUncVal_up.emplace("AbsoluteStat"    ,_jets_jetUnc_AbsoluteStat_up);
  _SourceUncVal_up.emplace("FlavorQCD"       ,_jets_jetUnc_FlavorQCD_up);
  _SourceUncVal_up.emplace("Fragmentation"   ,_jets_jetUnc_Fragmentation_up);
  _SourceUncVal_up.emplace("PileUpDataMC"    ,_jets_jetUnc_PileUpDataMC_up);
  _SourceUncVal_up.emplace("PileUpPtBB"      ,_jets_jetUnc_PileUpPtBB_up);
  _SourceUncVal_up.emplace("PileUpPtEC1"     ,_jets_jetUnc_PileUpPtEC1_up);
  _SourceUncVal_up.emplace("PileUpPtEC2"     ,_jets_jetUnc_PileUpPtEC2_up);
  _SourceUncVal_up.emplace("PileUpPtHF"      ,_jets_jetUnc_PileUpPtHF_up);
  _SourceUncVal_up.emplace("PileUpPtRef"     ,_jets_jetUnc_PileUpPtRef_up);
  _SourceUncVal_up.emplace("RelativeBal"     ,_jets_jetUnc_RelativeBal_up);
  _SourceUncVal_up.emplace("RelativeFSR"     ,_jets_jetUnc_RelativeFSR_up);
  _SourceUncVal_up.emplace("RelativeJEREC1"  ,_jets_jetUnc_RelativeJEREC1_up);
  _SourceUncVal_up.emplace("RelativeJEREC2"  ,_jets_jetUnc_RelativeJEREC2_up);
  _SourceUncVal_up.emplace("RelativeJERHF"   ,_jets_jetUnc_RelativeJERHF_up);
  _SourceUncVal_up.emplace("RelativePtBB"    ,_jets_jetUnc_RelativePtBB_up);
  _SourceUncVal_up.emplace("RelativePtEC1"   ,_jets_jetUnc_RelativePtEC1_up);
  _SourceUncVal_up.emplace("RelativePtEC2"   ,_jets_jetUnc_RelativePtEC2_up);
  _SourceUncVal_up.emplace("RelativePtHF"    ,_jets_jetUnc_RelativePtHF_up);
  _SourceUncVal_up.emplace("RelativeSample"  ,_jets_jetUnc_RelativeSample_up);
  _SourceUncVal_up.emplace("RelativeStatEC"  ,_jets_jetUnc_RelativeStatEC_up);
  _SourceUncVal_up.emplace("RelativeStatFSR" ,_jets_jetUnc_RelativeStatFSR_up);
  _SourceUncVal_up.emplace("RelativeStatHF"  ,_jets_jetUnc_RelativeStatHF_up);
  _SourceUncVal_up.emplace("SinglePionECAL"  ,_jets_jetUnc_SinglePionECAL_up);
  _SourceUncVal_up.emplace("SinglePionHCAL"  ,_jets_jetUnc_SinglePionHCAL_up);
  _SourceUncVal_up.emplace("TimePtEta"       ,_jets_jetUnc_TimePtEta_up);
  // Down variations
  _SourceUncVal_dw.emplace("AbsoluteFlavMap" ,_jets_jetUnc_AbsoluteFlavMap_dw);
  _SourceUncVal_dw.emplace("AbsoluteMPFBias" ,_jets_jetUnc_AbsoluteMPFBias_dw);
  _SourceUncVal_dw.emplace("AbsoluteSample"  ,_jets_jetUnc_AbsoluteSample_dw);
  _SourceUncVal_dw.emplace("AbsoluteScale"   ,_jets_jetUnc_AbsoluteScale_dw);
  _SourceUncVal_dw.emplace("AbsoluteStat"    ,_jets_jetUnc_AbsoluteStat_dw);
  _SourceUncVal_dw.emplace("FlavorQCD"       ,_jets_jetUnc_FlavorQCD_dw);
  _SourceUncVal_dw.emplace("Fragmentation"   ,_jets_jetUnc_Fragmentation_dw);
  _SourceUncVal_dw.emplace("PileUpDataMC"    ,_jets_jetUnc_PileUpDataMC_dw);
  _SourceUncVal_dw.emplace("PileUpPtBB"      ,_jets_jetUnc_PileUpPtBB_dw);
  _SourceUncVal_dw.emplace("PileUpPtEC1"     ,_jets_jetUnc_PileUpPtEC1_dw);
  _SourceUncVal_dw.emplace("PileUpPtEC2"     ,_jets_jetUnc_PileUpPtEC2_dw);
  _SourceUncVal_dw.emplace("PileUpPtHF"      ,_jets_jetUnc_PileUpPtHF_dw);
  _SourceUncVal_dw.emplace("PileUpPtRef"     ,_jets_jetUnc_PileUpPtRef_dw);
  _SourceUncVal_dw.emplace("RelativeBal"     ,_jets_jetUnc_RelativeBal_dw);
  _SourceUncVal_dw.emplace("RelativeFSR"     ,_jets_jetUnc_RelativeFSR_dw);
  _SourceUncVal_dw.emplace("RelativeJEREC1"  ,_jets_jetUnc_RelativeJEREC1_dw);
  _SourceUncVal_dw.emplace("RelativeJEREC2"  ,_jets_jetUnc_RelativeJEREC2_dw);
  _SourceUncVal_dw.emplace("RelativeJERHF"   ,_jets_jetUnc_RelativeJERHF_dw);
  _SourceUncVal_dw.emplace("RelativePtBB"    ,_jets_jetUnc_RelativePtBB_dw);
  _SourceUncVal_dw.emplace("RelativePtEC1"   ,_jets_jetUnc_RelativePtEC1_dw);
  _SourceUncVal_dw.emplace("RelativePtEC2"   ,_jets_jetUnc_RelativePtEC2_dw);
  _SourceUncVal_dw.emplace("RelativePtHF"    ,_jets_jetUnc_RelativePtHF_dw);
  _SourceUncVal_dw.emplace("RelativeSample"  ,_jets_jetUnc_RelativeSample_dw);
  _SourceUncVal_dw.emplace("RelativeStatEC"  ,_jets_jetUnc_RelativeStatEC_dw);
  _SourceUncVal_dw.emplace("RelativeStatFSR" ,_jets_jetUnc_RelativeStatFSR_dw);
  _SourceUncVal_dw.emplace("RelativeStatHF"  ,_jets_jetUnc_RelativeStatHF_dw);
  _SourceUncVal_dw.emplace("SinglePionECAL"  ,_jets_jetUnc_SinglePionECAL_dw);
  _SourceUncVal_dw.emplace("SinglePionHCAL"  ,_jets_jetUnc_SinglePionHCAL_dw);
  _SourceUncVal_dw.emplace("TimePtEta"       ,_jets_jetUnc_TimePtEta_dw);

  //----------------------------------------------------------------------
  // Analyze MC history. THIS HAS TO BE DONE BEFORE ANY RETURN STATEMENT
  // (eg skim or trigger), in order to update the gen counters correctly!!!

  if (event.isRealData() && !IsEmbed) {
    _DataMC_Type = DataMCType::Data;
  }
  Event_isRealData = event.isRealData();
  _npv = vertexs->size();
  if (theisMC) {
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByToken(thePUTag, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
        _PUNumInteractions  = PVI->getPU_NumInteractions();
        float nTrueInt = PVI->getTrueNumInteractions();
        _npu = nTrueInt;
        break;
      }
    }
  }


  // pile up information -- rho
  edm::Handle<double> rhoHandle;
  event.getByToken(theRhoTag, rhoHandle);
  _rho = *rhoHandle;

  edm::Handle<edm::TriggerResults> triggerBits;
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  if(DEBUG){
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      _trigger_name.push_back( (string) names.triggerName(i));
      _trigger_accept.push_back( (int) triggerBits->accept(i));
    }
  }

  if(theisMC) {
    Handle<LHEEventProduct> lheEventProduct;
    try {event.getByToken(theLHEPTag, lheEventProduct);} catch (...) {;}
    if (lheEventProduct.isValid())
      {
	const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
	std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
	double lheHt = 0.;
	int lheNOutPartons = 0;
	int lheNOutB = 0;
	int lheNOutC = 0;
	size_t numParticles = lheParticles.size();
	for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
	  int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
	  int status = lheEvent.ISTUP[idxParticle];
	  if ( status == 1 && ((absPdgId >= 1 &&  absPdgId<= 6) ||  absPdgId== 21) ) { // quarks and gluons
            //cout << "DEBUG: APDGID: " << absPdgId << endl;
            lheHt += TMath::Sqrt((lheParticles[idxParticle][0])*(lheParticles[idxParticle][0]) + (lheParticles[idxParticle][1])*(lheParticles[idxParticle][1])); // first entry is px, second py
            ++lheNOutPartons;
            if (absPdgId == 5) ++lheNOutB ;
            if (absPdgId == 4) ++lheNOutC ;
	  }
	}
	_lheHt = lheHt;
	_lheNOutPartons = lheNOutPartons;
	_lheNOutB = lheNOutB;
	_lheNOutC = lheNOutC;
      }
  }

  _triggerbit = myTriggerHelper->FindTriggerBit(event,foundPaths,indexOfPath,triggerBits);
  _metfilterbit = myTriggerHelper->FindMETBit(event, metFilterBits_);

  Long64_t tbit = _triggerbit;
  for(int itr=0;itr<myTriggerHelper->GetNTriggers();itr++) {
    if(myTriggerHelper->IsTriggerFired(tbit,itr)) hCounter->Fill(itr+3);
  }

  hYear->Fill(theYear);

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate>>candHandle;
  edm::Handle<edm::View<reco::Candidate>>dauHandle;
  edm::Handle<edm::View<pat::Jet>>SmearedjetHandle;
  edm::Handle<edm::View<pat::Jet>>SmearedjetDownHandle;
  edm::Handle<edm::View<pat::Jet>>SmearedjetUpHandle;
  edm::Handle<edm::View<pat::Jet>>fatjetHandle;
  edm::Handle<edm::ValueMap<float>>qgTaggerHandle;
  edm::Handle<BXVector<l1t::Tau>>L1TauHandle;
  edm::Handle<BXVector<l1t::Jet>>L1JetHandle;
  edm::Handle<pat::METCollection> PUPPImetHandle;
  edm::Handle<math::Error<2>::type> covHandle;
  edm::Handle<math::Error<2>::type> PUPPIcovHandle;
  edm::Handle<double> PUPPIMETsignificanceHandle;
  edm::Handle<GenFilterInfo> embeddingWeightHandle;
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle< bool > passecalBadCalibFilterUpdate ;
  edm::Handle< double > theprefweight;
  edm::Handle< double > theprefweightup;
  edm::Handle< double > theprefweightdown;

  // protect in case of events where trigger hasn't fired --> no collection created
  event.getByToken(theCandTag,candHandle);
  if (!candHandle.isValid()) return;
  event.getByToken(theCandTag,candHandle);
  event.getByToken(theSmearedJetTag,SmearedjetHandle);
  event.getByToken(theSmearedJetDownTag,SmearedjetDownHandle);
  event.getByToken(theSmearedJetUpTag,SmearedjetUpHandle);
  if(computeQGVar)
    event.getByToken(theQGTaggerTag,qgTaggerHandle);
  event.getByToken(theL1TauTag,L1TauHandle);
  event.getByToken(theL1JetTag,L1JetHandle);
  event.getByToken(theFatJetTag,fatjetHandle);
  event.getByToken(theLepTag,dauHandle);
  event.getByToken(thePUPPIMetTag,PUPPImetHandle);
  event.getByToken(thePUPPIMETCovTag,PUPPIcovHandle);
  event.getByToken(thePUPPIMETSignifTag,PUPPIMETsignificanceHandle);
  event.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
  event.getByToken(prefweight_token, theprefweight);
  event.getByToken(prefweightup_token, theprefweightup);
  event.getByToken(prefweightdown_token, theprefweightdown);

  if(theisMC || IsEmbed){
    edm::Handle<LHEEventProduct> lheeventinfo;
    event.getByToken(theLHEPTag,lheeventinfo);

    edm::Handle<GenEventInfoProduct> genEvt;
    event.getByToken(theGenTag,genEvt);
    _aMCatNLOweight=genEvt->weight();
    _MC_weight = _aMCatNLOweight; // duplicated

    if (lheeventinfo.isValid()) {
      _nup=lheeventinfo->hepeup().NUP;
      _nominal_wt = lheeventinfo->hepeup().XWGTUP;
      if (lheeventinfo->weights().size() > 6) // access weights only if weights() is filled
	{
	  _MC_weight_scale_muF0p5 = _aMCatNLOweight*(lheeventinfo->weights()[2].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 0.5 | muR = 1
	  _MC_weight_scale_muF2 = _aMCatNLOweight*(lheeventinfo->weights()[1].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 2 | muR = 1
	  _MC_weight_scale_muR0p5 = _aMCatNLOweight*(lheeventinfo->weights()[6].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 1 | muR = 0.5
	  _MC_weight_scale_muR2 = _aMCatNLOweight*(lheeventinfo->weights()[3].wgt)/(lheeventinfo->originalXWGTUP()); // muF = 1 | muR = 2
	}
      std::string TheoId="";
      for(unsigned int iweight=0;iweight<lheeventinfo->weights().size();iweight++){
	TheoId=lheeventinfo->weights()[iweight].id;
      	if(TheoId=="1")_TheoreticalScaleUncTab[0]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="2")_TheoreticalScaleUncTab[1]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="3")_TheoreticalScaleUncTab[2]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="4")_TheoreticalScaleUncTab[3]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="5")_TheoreticalScaleUncTab[4]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="6")_TheoreticalScaleUncTab[5]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="7")_TheoreticalScaleUncTab[6]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="8")_TheoreticalScaleUncTab[7]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="9")_TheoreticalScaleUncTab[8]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1001")_TheoreticalScaleUncTab[9]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1002")_TheoreticalScaleUncTab[10]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1003")_TheoreticalScaleUncTab[11]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1004")_TheoreticalScaleUncTab[12]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1005")_TheoreticalScaleUncTab[13]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1006")_TheoreticalScaleUncTab[14]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1007")_TheoreticalScaleUncTab[15]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1008")_TheoreticalScaleUncTab[16]=lheeventinfo->weights()[iweight].wgt;
	else if(TheoId=="1009")_TheoreticalScaleUncTab[17]=lheeventinfo->weights()[iweight].wgt;
      }
    }
    if (genEvt.isValid()) {
      for(unsigned int iweight=0;iweight<genEvt->weights().size();iweight++){
	_TheoreticalPSUnc.push_back(genEvt->weights()[iweight]);
      }
    }
  }

  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  const edm::View<reco::Candidate>* daus = dauHandle.product();
  const edm::View<pat::Jet>* Smearedjets = SmearedjetHandle.product();
  const edm::View<pat::Jet>* SmearedjetsDown = SmearedjetDownHandle.product();
  const edm::View<pat::Jet>* SmearedjetsUp = SmearedjetUpHandle.product();
  //const edm::View<pat::Jet>* fatjets = fatjetHandle.product();
  const pat::MET &PUPPImet = PUPPImetHandle->front();
  const BXVector<l1t::Tau>* L1Tau = L1TauHandle.product();
  //const BXVector<l1t::Jet>* L1Jet = L1JetHandle.product();

  _indexevents = event.id().event();

  //std::cout<<"--------------------------------------------------  " << event.id().event() <<" ------------"<<std::endl;
  _runNumber = event.id().run();
  //cout<<"_runNumber: "<<_runNumber<<endl;
  _lumi=event.luminosityBlock();
  _year=theYear;
  _PUPPImet = PUPPImet.pt();
  _PUPPImetphi = PUPPImet.phi();
  _PUPPIMETCov00 = (*PUPPIcovHandle)(0,0);
  _PUPPIMETCov10 = (*PUPPIcovHandle)(1,0);
  _PUPPIMETCov01 = _PUPPIMETCov10; // (1,0) is the only one saved
  _PUPPIMETCov11 = (*PUPPIcovHandle)(1,1);
  _PUPPIMETsignif = (*PUPPIMETsignificanceHandle);

  _puppimet_ex_JetEnUp = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_JetEnUp = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::JetEnUp,pat::MET::METCorrectionLevel::Type1);

  _puppimet_ex_JetEnDown = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_JetEnDown = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::JetEnDown,pat::MET::METCorrectionLevel::Type1);

  _puppimet_ex_UnclusteredEnUp = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_UnclusteredEnUp = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,pat::MET::METCorrectionLevel::Type1);

  _puppimet_ex_UnclusteredEnDown = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_UnclusteredEnDown = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,pat::MET::METCorrectionLevel::Type1);

  _puppimet_ex_JetResUp = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_JetResUp = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::JetResUp,pat::MET::METCorrectionLevel::Type1);

  _puppimet_ex_JetResDown = (*PUPPImetHandle)[0].shiftedPx(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);
  _puppimet_ey_JetResDown = (*PUPPImetHandle)[0].shiftedPy(pat::MET::METUncertainty::JetResDown,pat::MET::METCorrectionLevel::Type1);

  _passecalBadCalibFilterUpdate =  (*passecalBadCalibFilterUpdate );
  if (theisMC && (theYear==2016 || theYear==2017) && !IsEmbed)
    {
      _prefiringweight = (*theprefweight);
      _prefiringweightup =(*theprefweightup);
      _prefiringweightdown =(*theprefweightdown);
    }

  //Do all the stuff here
  //Compute the variables needed for the output and store them in the ntuple
  //if(DEBUG)printf("===New Event===\n");

  //Loop over generated b quarks
  //if(theisMC)FillbQuarks(event);
  if(theisMC || IsEmbed)
    {
      FillGenInfo(event); // gen particles
      if(theisMC && !IsEmbed)FillGenJetInfo(event); // gen jets
      fillMCTruth(event,eSetup);
    }
  //Loop of softleptons and fill them
  FillSoftLeptons(daus,event,eSetup,theFSR,Smearedjets,L1Tau);
  edm::Handle<pat::TauCollection> TauHandle;
  event.getByToken(theNewTauTag,TauHandle);

  //Loop on Jets
  // Accessing the JEC uncertainties
  //ak4
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  eSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc (JetCorPar);

  bool Run2016B = (_runNumber >= 272007 && _runNumber <= 275376);
  bool Run2016C = (_runNumber >= 275657 && _runNumber <= 276283);
  bool Run2016D = (_runNumber >= 276315 && _runNumber <= 276811);
  bool Run2016E = (_runNumber >= 276831 && _runNumber <= 277420);
  bool Run2016F = (_runNumber >= 277772 && _runNumber <= 278808);
  bool Run2016G = (_runNumber >= 278820 && _runNumber <= 280385);
  bool Run2016H = (_runNumber >= 280919 && _runNumber <= 284044);
  bool Run2017B = (_runNumber >= 297046 && _runNumber <= 299329);
  bool Run2017C = (_runNumber >= 299368 && _runNumber <= 302029);
  bool Run2017D = (_runNumber >= 302030 && _runNumber <= 303434);
  bool Run2017E = (_runNumber >= 303824 && _runNumber <= 304797);
  bool Run2017F = (_runNumber >= 305040 && _runNumber <= 306462);
  bool Run2018A = (_runNumber >= 315252 && _runNumber <= 316995);
  bool Run2018B = (_runNumber >= 316998 && _runNumber <= 319312);
  bool Run2018C = (_runNumber >= 319313 && _runNumber <= 320393);
  bool Run2018D = (_runNumber >= 320394 && _runNumber <= 325273);

  // JEC Regrouped uncertainty sources
  if (theYear == 2016)
    {
      for (const auto& source: m_jec_sources_regrouped_2016) {
	JetCorrectorParameters source_parameters_reduced("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Regrouped_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt", source);
	std::unique_ptr<JetCorrectionUncertainty> source_uncertainty_reduced(new JetCorrectionUncertainty(source_parameters_reduced));
	jecSourceUncRegroupedProviders.emplace(source, std::move(source_uncertainty_reduced));
      }
    }
  else if (theYear == 2017)
    {
      for (const auto& source: m_jec_sources_regrouped_2017) {
	JetCorrectorParameters source_parameters_reduced("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Regrouped_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt", source);
	std::unique_ptr<JetCorrectionUncertainty> source_uncertainty_reduced(new JetCorrectionUncertainty(source_parameters_reduced));
	jecSourceUncRegroupedProviders.emplace(source, std::move(source_uncertainty_reduced));
      }
    }
  else if (theYear == 2018)
    {
      for (const auto& source: m_jec_sources_regrouped_2018) {
	JetCorrectorParameters source_parameters_reduced("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt", source);
        //JetCorrectorParameters source_parameters_reduced("/opt/sbg/cms/safe1/cms/msessini/MuTauProducer/CMSSW_10_2_23/src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt", source);
	std::unique_ptr<JetCorrectionUncertainty> source_uncertainty_reduced(new JetCorrectionUncertainty(source_parameters_reduced));
	jecSourceUncRegroupedProviders.emplace(source, std::move(source_uncertainty_reduced));
      }
    }

  // Full JEC uncertainty sources

  if(theisMC)
    {
      if (theYear == 2016)
	{
	  for (const auto& source: m_jec_sources_2016) {
	    JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt", source);
	    std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	    jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	  }
	}
      else if (theYear == 2017)
	{
	  for (const auto& source: m_jec_sources_2017) {
	    JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt", source);
	    std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	    jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	  }
	}
      else if (theYear == 2018)
	{
	  for (const auto& source: m_jec_sources_2018) {
	    JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt", source);
            //JetCorrectorParameters source_parameters("/opt/sbg/cms/safe1/cms/msessini/MuTauProducer/CMSSW_10_2_23/src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt", source);
	    std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	    jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	  }
	}
    }
  else
    {
      if(Run2016B || Run2016C || Run2016D){
	for (const auto& source: m_jec_sources_2016) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Summer16_07Aug2017BCD_V11_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if (Run2016E || Run2016F){
	for (const auto& source: m_jec_sources_2016) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Summer16_07Aug2017EF_V11_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if (Run2016G || Run2016H){
	for (const auto& source: m_jec_sources_2016) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Summer16_07Aug2017GH_V11_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2017B){
	for (const auto& source: m_jec_sources_2017) {
	  //JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2017C){
	for (const auto& source: m_jec_sources_2017) {
	  //JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017C_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Fall17_17Nov2017C_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2017D || Run2017E){
	for (const auto& source: m_jec_sources_2017) {
	  //JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017DE_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Fall17_17Nov2017DE_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2017F){
	for (const auto& source: m_jec_sources_2017) {
	  //JetCorrectorParameters source_parameters("JECUncertaintySources/Fall17_17Nov2017F_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Fall17_17Nov2017F_V32_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2018A){
	for (const auto& source: m_jec_sources_2018) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_RunA_V19_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2018B){
	for (const auto& source: m_jec_sources_2018) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_RunB_V19_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2018C){
	for (const auto& source: m_jec_sources_2018) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_RunC_V19_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }else if(Run2018D){
	for (const auto& source: m_jec_sources_2018) {
	  JetCorrectorParameters source_parameters("src/LLRHiggsTauTau/NtupleProducer/data/JECUncertaintySources/Autumn18_RunD_V19_DATA_UncertaintySources_AK4PFchs.txt", source);
	  std::unique_ptr<JetCorrectionUncertainty> source_uncertainty(new JetCorrectionUncertainty(source_parameters));
	  jecSourceUncProviders.emplace(source, std::move(source_uncertainty));
	}
      }
    }
  //Fill Triggers
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  if(theYear==2016){
    //single
    MatchedTriggerNames.push_back("HLT_IsoMu22_v");
    MatchedTriggerNames.push_back("HLT_IsoMu22_eta2p1_v");
    MatchedTriggerNames.push_back("HLT_IsoTkMu22_v");
    MatchedTriggerNames.push_back("HLT_IsoTkMu22_eta2p1_v");
    //cross
    MatchedTriggerNames.push_back("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v");
    MatchedTriggerNames.push_back("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v");
  }
  else if(theYear==2017){
    //single
    MatchedTriggerNames.push_back("HLT_IsoMu24_v");
    MatchedTriggerNames.push_back("HLT_IsoMu27_v");
    //cross
    MatchedTriggerNames.push_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v");
  }
  else if(theYear==2018){
    //single
    MatchedTriggerNames.push_back("HLT_IsoMu24_v");
    MatchedTriggerNames.push_back("HLT_IsoMu27_v");
    //cross
    if(theisMC || (!theisMC && _runNumber>317509) || IsEmbed) MatchedTriggerNames.push_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v");
    else MatchedTriggerNames.push_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v");
    //if(!theisMC){
    //if(315974<=_runNumber<=317509 || IsEmbed) MatchedTriggerNames.push_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v");
    //if(_runNumber>317509 || IsEmbed) MatchedTriggerNames.push_back("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v");
    //}
  }

  std::vector<int> out;
  for(unsigned int i=0; i<_trigger_accept.size();i++){
    TString name=_trigger_name.at(i);
    for(unsigned int j=0; j<MatchedTriggerNames.size(); j++){
      if(name.Contains(MatchedTriggerNames.at(j))) out.push_back(i) ;
    }
  }


  bool trig=false;
  TriggerIndexVector=out;
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(_trigger_accept.at(TriggerIndexVector.at(itrig))){
      trig=true;
    }
  }

  //MET filters
  bool METFilters=true;

  if(theYear==2017 || theYear==2018){
    if(!_passecalBadCalibFilterUpdate){
      METFilters=false;
    }
  }
  if(_metfilterbit!=127) METFilters=false;
  if (!METFilters) return;

  //TauSpinner weights
  Handle <double> tauSpinnerWTEvenHandle;
  event.getByToken(theTauSpinnerWTEventag,tauSpinnerWTEvenHandle);
  const double tauSpinnerWTEven = (*tauSpinnerWTEvenHandle.product());
  //
  Handle <double> tauSpinnerWTOddHandle;
  event.getByToken(theTauSpinnerWTOddtag,tauSpinnerWTOddHandle);
  const double tauSpinnerWTOdd = (*tauSpinnerWTOddHandle.product());
  //
  Handle <double> tauSpinnerWTMMHandle;
  event.getByToken(theTauSpinnerWTMMtag,tauSpinnerWTMMHandle);
  const double tauSpinnerWTMM = (*tauSpinnerWTMMHandle.product());
  
  //Theoretical Uncertainties
  std::map<std::string, double> TheoreticalUncmap;
  TheoreticalUncmap["wMC"]        =  _MC_weight;
  TheoreticalUncmap["wPSISRUp"]   =  _TheoreticalPSUnc[6]/_MC_weight;
  TheoreticalUncmap["wPSISRDown"] =  _TheoreticalPSUnc[8]/_MC_weight;
  TheoreticalUncmap["wPSFSRUp"]   =  _TheoreticalPSUnc[7]/_MC_weight;
  TheoreticalUncmap["wPSFSRDown"] =  _TheoreticalPSUnc[9]/_MC_weight;
  TheoreticalUncmap["wScaleUp"]   =  _TheoreticalScaleUncTab[13]/_nominal_wt;
  TheoreticalUncmap["wScaleDown"] =  _TheoreticalScaleUncTab[17]/_nominal_wt;

  mySysHelper->GetCollections(cands, daus, Smearedjets, SmearedjetsUp, SmearedjetsDown);
  if(!Event_isRealData) {
    mySysHelper->GetGenInfo(theGenericTag, TheoreticalUncmap, _DataMC_Type);
    mySysHelper->FillGenTaus(MCSignalParticle_Tauidx, MCTauandProd_p4, MCTauandProd_charge, MCTauandProd_Vertex, MCTauandProd_pdgid, MCTau_JAK);
  }
  mySysHelper->GetEventInfo(IsEmbed, Event_isRealData, theisMC, _runNumber, _npu, _indexevents, _lumi);
  mySysHelper->GetDecayProducts(A1LVP, MuonTrack, _PFTauRefitPionsP4, _PFTauRefitPionsCharge);
  mySysHelper->GetPV(_RefitPVBS_x, _RefitPVBS_y, _RefitPVBS_z, _RefitPVBS_Cov, _VertexHashBS1, _VertexHashBS2, _LeptonHash);
  mySysHelper->GetMETCov(_PUPPIMETCov00, _PUPPIMETCov10, _PUPPIMETCov11);
  mySysHelper->GetJECUnc(_SourceUncValRegrouped_up, _SourceUncValRegrouped_dw, &jecSourceUncRegroupedProviders, &jecUnc);
  mySysHelper->GetTauSpinnerWeights(tauSpinnerWTEven, tauSpinnerWTOdd, tauSpinnerWTMM);
  /////
  if(Event_isRealData) {
    mySysHelper->FillTree(Nominal, "Nominal", "Nom", event, trig, _daughters_trgMatched, LeptonP4, PUPPImet, _lheNOutPartons);
  }
  else {
    for(std::map<std::string, TTree*>::iterator iMap = SystematicsMap.begin(); iMap != SystematicsMap.end(); ++iMap) {
      mySysHelper->ResetVariables();
      if(iMap->first == "Nominal") mySysHelper->FillTree(iMap->second, "Nominal", "Nom", event, trig, _daughters_trgMatched, LeptonP4, PUPPImet, _lheNOutPartons);
      else {
        auto foundUp = (iMap->first).find("Up");
        auto foundDown = (iMap->first).find("Down");
        std::string Sys = ""; //Name of systematic
        std::string Unc = ""; //Up or Down
        if(foundUp != std::string::npos) {
          Sys = iMap->first.substr(0, foundUp);
          Unc = iMap->first.substr(foundUp);
        }
        else {
	  Sys = iMap->first.substr(0, foundDown);
          Unc = iMap->first.substr(foundDown);
        }
        //
        mySysHelper->FillTree(iMap->second, Sys, Unc, event, trig, _daughters_trgMatched, LeptonP4, PUPPImet, _lheNOutPartons);
      }
    }
  }
}


//Fill jets quantities
int HTauTauNtuplizer::FillJet(const edm::View<pat::Jet> *jets, const edm::View<pat::Jet> *jetsDown,const edm::View<pat::Jet> *jetsUp,const edm::Event& event, edm::EventSetup const& iSetup, JetCorrectionUncertainty* jecUnc, myJECMap* jecSourceUncProviders, myJECMap* jecSourceUncRegroupedProviders,const reco::Candidate* cand1,const reco::Candidate* cand2, int ipair) {

  // TriggerBits and TriggerObjets (for VBF trigger matching)
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);

  // Save a subsample of TOs with (type==85) && (hasLeadFilter || hasSubLeadFilter)
  std::vector<pat::TriggerObjectStandAlone> goodTOs;
  for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto)
    {
      // Get the TO and unpack the filter labels
      pat::TriggerObjectStandAlone TO = triggerObjects->at(idxto);
      TO.unpackFilterLabels(event,*triggerBits);

      if ( (TO.type(85)==true) && ( TO.hasFilterLabel("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20") || TO.hasFilterLabel("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20") || TO.hasFilterLabel("hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTauHPS20") || TO.hasFilterLabel("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTauHPS20") ) )
	goodTOs.push_back(TO);
    }

  // Getting the primary Vertex (FRA 2017)
  Handle<vector<reco::Vertex> >  vertexs;
  event.getByToken(theVtxTag,vertexs);
  const auto & pv = (*vertexs)[0]; // get the firstPV


  // Getting the secondaryVertexCollection (FRA 2017)
  Handle<edm::View<reco::VertexCompositePtrCandidate>> secVtxsHandle;
  event.getByToken(theSecVtxTag,secVtxsHandle);
  if (!secVtxsHandle.isValid())
    {
      throw cms::Exception("ProductNotValid") << "slimmedSecondaryVertices product not valid";
    }
  const edm::View<reco::VertexCompositePtrCandidate> * secVtxs = secVtxsHandle.product();

  // Getting the Rho for Jet Energy Resolution
  edm::Handle<double> rhoJERHandle;
  event.getByToken(theRhoForJERTag, rhoJERHandle);

  // Initialize Jet Energy Resolution (FRA 2017)
  JME::JetResolution JERresolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");

  vector <pair<float, int>> softLeptInJet; // pt, idx
  vector<pat::Jet> jetsVect;
  vector<pat::Jet> njetspt20Vect;
  vector<pat::Jet> btagjetsVect;
  vector<pat::Jet> btagjetsCSVVect;
  vector<bool> njetspt20GapVect;

  for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet) {

    float jecFactor = ijet->jecFactor("Uncorrected");
    float jetRawPt = jecFactor * ijet->pt();

    float pileupjetidMVA = ijet->userFloat("pileupJetId:fullDiscriminant");

    //PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float MUF = ijet->muonEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    int NumNeutralParticles =ijet->neutralMultiplicity();
    int chargedMult = ijet->chargedMultiplicity();
    int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());

    _jets_chEmEF .push_back(CEMF);
    _jets_chHEF  .push_back(CHF);
    _jets_nEmEF  .push_back(NEMF);
    _jets_nHEF   .push_back(NHF);
    _jets_chMult .push_back(chargedMult);
    _jets_neMult .push_back(NumNeutralParticles);
    _jets_MUF    .push_back(MUF);

    // JetID
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data

    // 2016 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    // 2017 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    // 2018 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018

    //int jetid=0;


    bool looseJetID = false;
    bool tightJetID = false;
    bool tightLepVetoJetID = false;

    //2016 data
    if (theYear == 2016)
      {
        if (absjeta <= 2.7)
          {
            looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
            tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
            tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4) );
          }
        else if (absjeta <= 3.0)
          {
            looseJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
            tightJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
          }
        else
          {
            looseJetID = (NEMF<0.90 && NumNeutralParticles>10 );
            tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
          }

        if(looseJetID) _looseJetID.push_back(true);
        else _looseJetID.push_back(false);
        if(tightJetID) _tightJetID.push_back(true);
        else _tightJetID.push_back(false);
        if(tightLepVetoJetID) _tightLepVetoJetID.push_back(true);
        else _tightLepVetoJetID.push_back(false);
      }
   // 2017 data
    else if (theYear == 2017)
      {
        if (absjeta <= 2.7)
          {
            tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0) || absjeta>2.4) );
            tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.80) || absjeta>2.4) );
          }
        else if (absjeta <= 3.0)
          {
            tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
          }
        else
          {
            tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 );
          }
        if(tightJetID) _tightJetID.push_back(true);
        else _tightJetID.push_back(false);
        if(tightLepVetoJetID) _tightLepVetoJetID.push_back(true);
        else _tightLepVetoJetID.push_back(false);
      }
    // 2018 data
    else if (theYear == 2018)
      {
        if (absjeta <= 2.6)
          {
            tightJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 );
            tightLepVetoJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80 );
          }
        else if (absjeta>2.6 && absjeta <= 2.7)
          {
            tightJetID = ( NHF<0.90 && NEMF<0.99 && CHM>0 );
            tightLepVetoJetID = ( NHF<0.90 && NEMF<0.99 && MUF<0.8 && CHM>0 && CEMF<0.80 );
          }
        else if (absjeta>2.7 && absjeta <= 3.0)
          {
            tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
          }
        else
          {
            tightJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );
          }
        if(tightJetID) _tightJetID.push_back(true);
        else _tightJetID.push_back(false);
        if(tightLepVetoJetID) _tightLepVetoJetID.push_back(true);
        else _tightLepVetoJetID.push_back(false);
      }

    if(ijet->pt()<20 || abs(ijet->eta())>4.7 || deltaR(ijet->p4(),cand1->p4())<0.5 || deltaR(ijet->p4(),cand2->p4())<0.5 || !tightJetID)continue;

    if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;

    if(theYear==2016){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.6321 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsVect.push_back(*ijet);}} //Medium
    if(theYear==2017){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4941 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsVect.push_back(*ijet);}} //Medium
    if(theYear==2018){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4184 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsVect.push_back(*ijet);}} //Medium

    if(ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsCSVVect.push_back(*ijet);} //Medium
    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){njetspt20Vect.push_back(*ijet);}

    _jets_px.push_back( (float) ijet->px());
    _jets_py.push_back( (float) ijet->py());
    _jets_pz.push_back( (float) ijet->pz());
    _jets_e.push_back( (float) ijet->energy());
    _jets_mT.push_back( (float) ijet->mt());
    _jets_Flavour.push_back(ijet->partonFlavour());
    _jets_HadronFlavour.push_back(ijet->hadronFlavour());

    //FRA: new way to access these variables (in 2017)
    Float_t vtxPt   = 0.0;
    Float_t vtxMass = 0.0;
    Float_t vtx3dL  = 0.0;
    Float_t vtxNtrk = 0.0;
    Float_t vtx3deL = 0.0;

    VertexDistance3D vdist;
    float maxFoundSignificance=0;
    for(edm::View<reco::VertexCompositePtrCandidate>::const_iterator isv = secVtxs->begin(); isv!=secVtxs->end(); ++isv)
      {
	GlobalVector flightDir(isv->vertex().x() - pv.x(), isv->vertex().y() - pv.y(), isv->vertex().z() - pv.z());
	GlobalVector jetDir(ijet->px(),ijet->py(),ijet->pz());
	if( deltaR2( flightDir, jetDir ) < 0.09 )
	  {
	    Measurement1D dl= vdist.distance(pv,VertexState(RecoVertex::convertPos(isv->position()),RecoVertex::convertError(isv->error())));
	    if(dl.significance() > maxFoundSignificance)
	      {
		maxFoundSignificance=dl.significance();
		vtxPt   = isv->pt();
		vtxMass = isv->p4().M();
		vtx3dL  = dl.value();
		vtx3deL = dl.error();
		vtxNtrk = isv->numberOfSourceCandidatePtrs();
	      }
	  }
      }


    _jets_vtxPt.  push_back(vtxPt);
    _jets_vtxMass.push_back(vtxMass);
    _jets_vtx3dL. push_back(vtx3dL);
    _jets_vtxNtrk.push_back(vtxNtrk);
    _jets_vtx3deL.push_back(vtx3deL);

    _bdiscr.push_back(ijet->bDiscriminator("pfJetProbabilityBJetTags"));
    _bdiscr2.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    _bdiscr3.push_back(ijet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    // DeepCSV
    _bdiscr4.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
    _bdiscr5.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    _bdiscr6.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
    _bdiscr7.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probc"));
    _bdiscr8.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probcc"));

    // DeepFlavor
    _bdiscr9.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probb"));
    _bdiscr10.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    _bdiscr11.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    _bdiscr12.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probc"));
    _bdiscr13.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    _bdiscr14.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probg"));

    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && tightJetID){njetspt20GapVect.push_back(true);}
    else if(ijet->pt()>20 && abs(ijet->eta())<4.7) {njetspt20GapVect.push_back(false);}

    _jets_area.push_back (ijet->jetArea());
    _jetrawf.push_back(jecFactor);

    _pileupjetidMVA.push_back(pileupjetidMVA);

    // loop on jet contituents to retrieve info for b jet regression
    int nDau = ijet -> numberOfDaughters();
    //cout << "JET: " << (ijet - jets->begin()) << " N daught: " << nDau << endl;

    float leadTrackPt = 0.;
    softLeptInJet.clear();
    for (int iDau = 0; iDau < nDau; ++iDau)
      {
	// pdg id for packed pf candidates meaning is:
	// the particle charge and pdgId: 11, 13, 22 for ele/mu/gamma, 211 for charged hadrons, 130 for neutral hadrons, 1 and 2 for hadronic and em particles in HF.
	const Candidate * dau = ijet->daughter(iDau);
	if (abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13)
	  {
	    softLeptInJet.push_back( make_pair(dau->pt(), iDau) );
	  }

	if (dau->charge() != 0 ) // tracks -> charged
	  {
	    float ptBuf = dau->pt();
	    if (ptBuf > leadTrackPt) leadTrackPt = ptBuf;
	  }
      }

    _jets_leadTrackPt.push_back(leadTrackPt);
    float leptonPtRel = -1.;
    float leptonPt = -1.;
    float leptonDeltaR = -1.;
    int softLeptIdx = -1;
    if (softLeptInJet.size() > 0)
      {
	sort(softLeptInJet.begin(), softLeptInJet.end());
	softLeptIdx = softLeptInJet.back().second;
      }
    if (softLeptIdx >= 0)
      {
	const Candidate * dau = ijet->daughter(softLeptIdx);
	leptonPtRel = dau->pt() / ijet->pt() ;
	leptonPt = dau->pt();
	leptonDeltaR = deltaR(*dau, *ijet) ;
      }
    _jets_leptonPtRel .push_back (leptonPtRel);
    _jets_leptonPt    .push_back (leptonPt);
    _jets_leptonDeltaR.push_back (leptonDeltaR);


    jecUnc->setJetEta(ijet->eta());
    jecUnc->setJetPt(ijet->pt()); // here you must use the CORRECTED jet pt
    _jets_jecUnc.push_back(jecUnc->getUncertainty(true));

    // JEC Regrouped uncertainty sources
    for (myJECMap::iterator it=jecSourceUncRegroupedProviders->begin(); it!=jecSourceUncRegroupedProviders->end(); ++it)
      {
	// up variations
	it->second->setJetEta(ijet->eta());
	it->second->setJetPt(ijet->pt());
	float uncertainty_up = it->second->getUncertainty(true);
	_SourceUncValRegrouped_up[it->first].push_back(uncertainty_up);

	// down variations
	it->second->setJetEta(ijet->eta());
	it->second->setJetPt(ijet->pt());
	float uncertainty_dw = it->second->getUncertainty(false);
	_SourceUncValRegrouped_dw[it->first].push_back(uncertainty_dw);
      }

    // JEC uncertainties sources
    for (myJECMap::iterator it=jecSourceUncProviders->begin(); it!=jecSourceUncProviders->end(); ++it)
      {
	// up variations
	it->second->setJetEta(ijet->eta());
	it->second->setJetPt(ijet->pt());
	float uncertainty_up = it->second->getUncertainty(true);
	_SourceUncVal_up[it->first].push_back(uncertainty_up);

	// down variations
	it->second->setJetEta(ijet->eta());
	it->second->setJetPt(ijet->pt());
	float uncertainty_dw = it->second->getUncertainty(false);
	_SourceUncVal_dw[it->first].push_back(uncertainty_dw);
      }

    // Jet Energy Resolution
    JME::JetParameters JER_parameters;
    JER_parameters.setJetPt(ijet->pt());
    JER_parameters.setJetEta(ijet->eta());
    JER_parameters.setRho(*rhoJERHandle);
    _jets_JER.push_back( JERresolution.getResolution(JER_parameters) * ijet->energy() ); // JER*energy because cmssw gives the % of JER, while KinFit wants resolution in GeV

    // VBF trigger matching
    Long64_t VBFleadFilterMatch    = 0;
    Long64_t VBFsubleadFilterMatch = 0;
    // loop on the TOs that are already of type==85 and passed one of the VBF jet filters
    for (size_t idxto = 0; idxto < goodTOs.size(); ++idxto)
      {
	// Get the TO
	pat::TriggerObjectStandAlone obj = goodTOs.at(idxto);

	// bools to save if it has the filter or not
	bool isVBFsubleadFilterMatched = false;
	bool isVBFleadFilterMatched = false;

	// DeltaR2 match between jet and trigger object
	if(deltaR2(obj,*ijet)<0.25)
	  {
	    // Unpacking Filter Labels and Path Names
	    obj.unpackFilterLabels(event,*triggerBits);
	    obj.unpackPathNames(names);

	    //Get HLT path names
	    std::vector<std::string> pathNamesAll  = obj.pathNames(false);

	    // Loop on the HLT path names in the Trigger Object
	    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
	      {
		int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
		if (triggerbit < 0) continue ; // not a path I want to save

		triggerMapper trgmap = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);

		// TO filter labels
		const std::vector<std::string>& vLabels = obj.filterLabels();

		// Loop on Leg3 filter (2 jets with pt40) ---> subleadFilter
		if (trgmap.GetNfiltersleg3()>0)
		  {
		    for(int ifilt=0;ifilt<trgmap.GetNfiltersleg3();ifilt++)
		      {
			string label = trgmap.GetfilterVBF(true,ifilt);
			if (label.empty()) continue;
			if (find(vLabels.begin(), vLabels.end(), label) != vLabels.end()) isVBFsubleadFilterMatched=true;
		      }
		  }
		else isVBFsubleadFilterMatched = false;

		// Loop on Leg4 filter (1 jet with pt115)  ---> leadFilter
		if (trgmap.GetNfiltersleg4()>0)
		  {
		    for(int ifilt=0;ifilt<trgmap.GetNfiltersleg4();ifilt++) //change to leg4
		      {
			string label = trgmap.GetfilterVBF(false,ifilt);
			if (label.empty()) continue;
			if (find(vLabels.begin(), vLabels.end(), label) != vLabels.end()) isVBFleadFilterMatched=true;
		      }
		  }
		else isVBFleadFilterMatched = false;

	      } // end loop on HLT paths

	  } // end if dR2<0.25

	if(isVBFsubleadFilterMatched) VBFsubleadFilterMatch |= (Long64_t(1) <<idxto);
	if(isVBFleadFilterMatched)    VBFleadFilterMatch    |= (Long64_t(1) <<idxto);

      } //end loop on goodTOs

    // Fill branches with result of VBF trig matching
    _jets_VBFleadFilterMatch.push_back(VBFleadFilterMatch);
    _jets_VBFsubleadFilterMatch.push_back(VBFsubleadFilterMatch);

  }// end loop on jets

  vector <pair<float, int>> softLeptInJetDown; // pt, idx
  vector<pat::Jet> jetsDownVect;
  vector<pat::Jet> njetspt20DownVect;
  vector<pat::Jet> btagjetsDownVect;
  vector<pat::Jet> btagjetsCSVDownVect;
  vector<bool> njetspt20GapDownVect;

  for(edm::View<pat::Jet>::const_iterator ijet = jetsDown->begin(); ijet!=jetsDown->end();++ijet) {

    float jecFactor = ijet->jecFactor("Uncorrected");
    float jetRawPt = jecFactor * ijet->pt();

    float pileupjetidMVA = ijet->userFloat("pileupJetId:fullDiscriminant");

    //PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float MUF = ijet->muonEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    int NumNeutralParticles =ijet->neutralMultiplicity();
    int chargedMult = ijet->chargedMultiplicity();
    int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());

    _jetsDown_chEmEF .push_back(CEMF);
    _jetsDown_chHEF  .push_back(CHF);
    _jetsDown_nEmEF  .push_back(NEMF);
    _jetsDown_nHEF   .push_back(NHF);
    _jetsDown_chMult .push_back(chargedMult);
    _jetsDown_neMult .push_back(NumNeutralParticles);
    _jetsDown_MUF    .push_back(MUF);


    // JetID
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data

    // 2016 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    // 2017 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    // 2018 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018

    //int jetid=0;

    bool looseJetID = false;
    bool tightJetID = false;
    bool tightLepVetoJetID = false;

    //2016 data
    if (theYear == 2016)
      {
  if (absjeta <= 2.7)
    {
      looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
      tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4) );
    }
  else if (absjeta <= 3.0)
    {
      looseJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
      tightJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
    }
  else
    {
      looseJetID = (NEMF<0.90 && NumNeutralParticles>10 );
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    }

  if(looseJetID) _looseJetIDDown.push_back(true);
  else _looseJetIDDown.push_back(false);
  if(tightJetID) _tightJetIDDown.push_back(true);
  else _tightJetIDDown.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDDown.push_back(true);
  else _tightLepVetoJetIDDown.push_back(false);
      }
    // 2017 data
    else if (theYear == 2017)
      {
  if (absjeta <= 2.7)
    {
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0) || absjeta>2.4) );
      tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.80) || absjeta>2.4) );
    }
  else if (absjeta <= 3.0)
    {
      tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    }
  else
    {
      tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 );
    }
  if(tightJetID) _tightJetIDDown.push_back(true);
  else _tightJetIDDown.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDDown.push_back(true);
  else _tightLepVetoJetIDDown.push_back(false);
      }
    // 2018 data
    else if (theYear == 2018)
      {
  if (absjeta <= 2.6)
    {
      tightJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 );
      tightLepVetoJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80 );
    }
  else if (absjeta>2.6 && absjeta <= 2.7)
    {
      tightJetID = ( NHF<0.90 && NEMF<0.99 && CHM>0 );
      tightLepVetoJetID = ( NHF<0.90 && NEMF<0.99 && MUF<0.8 && CHM>0 && CEMF<0.80 );
    }
  else if (absjeta>2.7 && absjeta <= 3.0)
    {
      tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    }
  else
    {
      tightJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );
    }
  if(tightJetID) _tightJetIDDown.push_back(true);
  else _tightJetIDDown.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDDown.push_back(true);
  else _tightLepVetoJetIDDown.push_back(false);
      }

    if(ijet->pt()<20 || abs(ijet->eta())>4.7 || deltaR(ijet->p4(),cand1->p4())<0.5 || deltaR(ijet->p4(),cand2->p4())<0.5 || !tightJetID)continue;

    if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;

    if(theYear==2016){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.6321 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsDownVect.push_back(*ijet);}} //Medium
    if(theYear==2017){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4941 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsDownVect.push_back(*ijet);}} //Medium
    if(theYear==2018){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4184 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsDownVect.push_back(*ijet);}} //Medium

    if(ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsCSVDownVect.push_back(*ijet);} //Medium

    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){njetspt20DownVect.push_back(*ijet);}

    _jetsDown_px.push_back( (float) ijet->px());
    _jetsDown_py.push_back( (float) ijet->py());
    _jetsDown_pz.push_back( (float) ijet->pz());
    _jetsDown_e.push_back( (float) ijet->energy());
    _jetsDown_mT.push_back( (float) ijet->mt());
    _jetsDown_Flavour.push_back(ijet->partonFlavour());
    _jetsDown_HadronFlavour.push_back(ijet->hadronFlavour());


    _bdiscrDown.push_back(ijet->bDiscriminator("pfJetProbabilityBJetTags"));
    _bdiscr2Down.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    _bdiscr3Down.push_back(ijet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    // DeepCSV
    _bdiscr4Down.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
    _bdiscr5Down.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    _bdiscr6Down.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
    _bdiscr7Down.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probc"));
    _bdiscr8Down.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probcc"));

    // DeepFlavor
    _bdiscr9Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probb"));
    _bdiscr10Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    _bdiscr11Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    _bdiscr12Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probc"));
    _bdiscr13Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    _bdiscr14Down.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probg"));

    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && tightJetID){njetspt20GapDownVect.push_back(true);}
    else if(ijet->pt()>20 && abs(ijet->eta())<4.7) {njetspt20GapDownVect.push_back(false);}

    _jetsDown_area.push_back (ijet->jetArea());
    _jetrawfDown.push_back(jecFactor);

    _pileupjetidMVADown.push_back(pileupjetidMVA);

    // loop on jet contituents to retrieve info for b jet regression
    int nDau = ijet -> numberOfDaughters();

    float leadTrackPt = 0.;
    softLeptInJetDown.clear();
    for (int iDau = 0; iDau < nDau; ++iDau)
      {
	// pdg id for packed pf candidates meaning is:
	// the particle charge and pdgId: 11, 13, 22 for ele/mu/gamma, 211 for charged hadrons, 130 for neutral hadrons, 1 and 2 for hadronic and em particles in HF.
	const Candidate * dau = ijet->daughter(iDau);
	if (abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13)
	  {
	    softLeptInJetDown.push_back( make_pair(dau->pt(), iDau) );
	  }

	if (dau->charge() != 0 ) // tracks -> charged
	  {
	    float ptBuf = dau->pt();
	    if (ptBuf > leadTrackPt) leadTrackPt = ptBuf;
	  }

      }

    _jetsDown_leadTrackPt.push_back(leadTrackPt);
    float leptonPtRel = -1.;
    float leptonPt = -1.;
    float leptonDeltaR = -1.;
    int softLeptIdx = -1;
    if (softLeptInJetDown.size() > 0)
      {
	sort(softLeptInJetDown.begin(), softLeptInJetDown.end());
	softLeptIdx = softLeptInJetDown.back().second;
      }
    if (softLeptIdx >= 0)
      {
	const Candidate * dau = ijet->daughter(softLeptIdx);
	leptonPtRel = dau->pt() / ijet->pt() ;
	leptonPt = dau->pt();
	leptonDeltaR = deltaR(*dau, *ijet) ;
      }
    _jetsDown_leptonPtRel .push_back (leptonPtRel);
    _jetsDown_leptonPt    .push_back (leptonPt);
    _jetsDown_leptonDeltaR.push_back (leptonDeltaR);

  }// end loop on jets


  vector <pair<float, int>> softLeptInJetUp; // pt, idx
  vector<pat::Jet> jetsUpVect;
  vector<pat::Jet> njetspt20UpVect;
  vector<pat::Jet> btagjetsUpVect;
  vector<pat::Jet> btagjetsCSVUpVect;
  vector<bool> njetspt20GapUpVect;

  for(edm::View<pat::Jet>::const_iterator ijet = jetsUp->begin(); ijet!=jetsUp->end();++ijet) {

    float jecFactor = ijet->jecFactor("Uncorrected");
    float jetRawPt = jecFactor * ijet->pt();

    float pileupjetidMVA = ijet->userFloat("pileupJetId:fullDiscriminant");

    //PF jet ID
    float NHF = ijet->neutralHadronEnergyFraction();
    float NEMF = ijet->neutralEmEnergyFraction();
    float CHF = ijet->chargedHadronEnergyFraction();
    float MUF = ijet->muonEnergyFraction();
    float CEMF = ijet->chargedEmEnergyFraction();
    int NumNeutralParticles =ijet->neutralMultiplicity();
    int chargedMult = ijet->chargedMultiplicity();
    int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
    float CHM = ijet->chargedMultiplicity();
    float absjeta = fabs(ijet->eta());

    _jetsUp_chEmEF .push_back(CEMF);
    _jetsUp_chHEF  .push_back(CHF);
    _jetsUp_nEmEF  .push_back(NEMF);
    _jetsUp_nHEF   .push_back(NHF);
    _jetsUp_chMult .push_back(chargedMult);
    _jetsUp_neMult .push_back(NumNeutralParticles);
    _jetsUp_MUF    .push_back(MUF);


    // JetID
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data

    // 2016 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    // 2017 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    // 2018 Data: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018

    //int jetid=0;

    bool looseJetID = false;
    bool tightJetID = false;
    bool tightLepVetoJetID = false;

    //2016 data
    if (theYear == 2016)
      {
  if (absjeta <= 2.7)
    {
      looseJetID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) );
      tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4) );
    }
  else if (absjeta <= 3.0)
    {
      looseJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
      tightJetID = (NEMF>0.01 && NHF<0.98 && NumNeutralParticles>2 );
    }
  else
    {
      looseJetID = (NEMF<0.90 && NumNeutralParticles>10 );
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10 );
    }

  if(looseJetID) _looseJetIDUp.push_back(true);
  else _looseJetIDUp.push_back(false);
  if(tightJetID) _tightJetIDUp.push_back(true);
  else _tightJetIDUp.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDUp.push_back(true);
  else _tightLepVetoJetIDUp.push_back(false);
      }
    // 2017 data
    else if (theYear == 2017)
      {
  if (absjeta <= 2.7)
    {
      tightJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0) || absjeta>2.4) );
      tightLepVetoJetID = ( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.80) || absjeta>2.4) );
    }
  else if (absjeta <= 3.0)
    {
      tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    }
  else
    {
      tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 );
    }
  if(tightJetID) _tightJetIDUp.push_back(true);
  else _tightJetIDUp.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDUp.push_back(true);
  else _tightLepVetoJetIDUp.push_back(false);
      }
    // 2018 data
    else if (theYear == 2018)
      {
  if (absjeta <= 2.6)
    {
      tightJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 );
      tightLepVetoJetID = ( NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80 );
    }
  else if (absjeta>2.6 && absjeta <= 2.7)
    {
      tightJetID = ( NHF<0.90 && NEMF<0.99 && CHM>0 );
      tightLepVetoJetID = ( NHF<0.90 && NEMF<0.99 && MUF<0.8 && CHM>0 && CEMF<0.80 );
    }
  else if (absjeta>2.7 && absjeta <= 3.0)
    {
      tightJetID = ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    }
  else
    {
      tightJetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );
    }
  if(tightJetID) _tightJetIDUp.push_back(true);
  else _tightJetIDUp.push_back(false);
  if(tightLepVetoJetID) _tightLepVetoJetIDUp.push_back(true);
  else _tightLepVetoJetIDUp.push_back(false);
      }

    if(ijet->pt()<20 || abs(ijet->eta())>4.7 || deltaR(ijet->p4(),cand1->p4())<0.5 || deltaR(ijet->p4(),cand2->p4())<0.5 || !tightJetID)continue;

    if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;

    if(theYear==2016){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.6321 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsUpVect.push_back(*ijet);}} //Medium
    if(theYear==2017){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4941 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsUpVect.push_back(*ijet);}} //Medium
    if(theYear==2018){if((ijet->bDiscriminator("pfDeepCSVJetTags:probb") + ijet->bDiscriminator("pfDeepCSVJetTags:probbb"))>0.4184 && ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsUpVect.push_back(*ijet);}} //Medium

    if(ijet->pt()>20 && abs(ijet->eta())<2.4 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){btagjetsCSVUpVect.push_back(*ijet);} //Medium

    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && deltaR(ijet->p4(),cand1->p4())>0.5 && deltaR(ijet->p4(),cand2->p4())>0.5){njetspt20UpVect.push_back(*ijet);}

    _jetsUp_px.push_back( (float) ijet->px());
    _jetsUp_py.push_back( (float) ijet->py());
    _jetsUp_pz.push_back( (float) ijet->pz());
    _jetsUp_e.push_back( (float) ijet->energy());
    _jetsUp_mT.push_back( (float) ijet->mt());
    _jetsUp_Flavour.push_back(ijet->partonFlavour());
    _jetsUp_HadronFlavour.push_back(ijet->hadronFlavour());


    _bdiscrUp.push_back(ijet->bDiscriminator("pfJetProbabilityBJetTags"));
    _bdiscr2Up.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    _bdiscr3Up.push_back(ijet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    // DeepCSV
    _bdiscr4Up.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
    _bdiscr5Up.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    _bdiscr6Up.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
    _bdiscr7Up.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probc"));
    _bdiscr8Up.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probcc"));

    // DeepFlavor
    _bdiscr9Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probb"));
    _bdiscr10Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    _bdiscr11Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    _bdiscr12Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probc"));
    _bdiscr13Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    _bdiscr14Up.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probg"));

    if(ijet->pt()>20 && abs(ijet->eta())<4.7 && tightJetID){njetspt20GapUpVect.push_back(true);}
    else if(ijet->pt()>20 && abs(ijet->eta())<4.7) {njetspt20GapUpVect.push_back(false);}

    _jetsUp_area.push_back (ijet->jetArea());
    _jetrawfUp.push_back(jecFactor);

    _pileupjetidMVAUp.push_back(pileupjetidMVA);

    // loop on jet contituents to retrieve info for b jet regression
    int nDau = ijet -> numberOfDaughters();

    float leadTrackPt = 0.;
    softLeptInJetUp.clear();
    for (int iDau = 0; iDau < nDau; ++iDau)
      {
	const Candidate * dau = ijet->daughter(iDau);
	if (abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13)
	  {
	    softLeptInJetUp.push_back( make_pair(dau->pt(), iDau) );
	  }

	if (dau->charge() != 0 ) // tracks -> charged
	  {
	    float ptBuf = dau->pt();
	    if (ptBuf > leadTrackPt) leadTrackPt = ptBuf;
	  }

      }

    _jetsUp_leadTrackPt.push_back(leadTrackPt);
    float leptonPtRel = -1.;
    float leptonPt = -1.;
    float leptonDeltaR = -1.;
    int softLeptIdx = -1;
    if (softLeptInJetUp.size() > 0)
      {
	sort(softLeptInJetUp.begin(), softLeptInJetUp.end());
	softLeptIdx = softLeptInJetUp.back().second;
      }
    if (softLeptIdx >= 0)
      {
	const Candidate * dau = ijet->daughter(softLeptIdx);
	leptonPtRel = dau->pt() / ijet->pt() ;
	leptonPt = dau->pt();
	leptonDeltaR = deltaR(*dau, *ijet) ;
      }
    _jetsUp_leptonPtRel .push_back (leptonPtRel);
    _jetsUp_leptonPt    .push_back (leptonPt);
    _jetsUp_leptonDeltaR.push_back (leptonDeltaR);

  }// end loop on jets

  if(njetspt20Vect.size()>1)std::sort(njetspt20Vect.begin(),njetspt20Vect.end(),ComparePairsbyPt);

  if(btagjetsVect.size()>1)std::sort(btagjetsVect.begin(),btagjetsVect.end(),ComparePairsbyPt);
  if(btagjetsCSVVect.size()>1)std::sort(btagjetsCSVVect.begin(),btagjetsCSVVect.end(),ComparePairsbyPt);

  if(njetspt20DownVect.size()>1)std::sort(njetspt20DownVect.begin(),njetspt20DownVect.end(),ComparePairsbyPt);

  if(btagjetsDownVect.size()>1)std::sort(btagjetsDownVect.begin(),btagjetsDownVect.end(),ComparePairsbyPt);
  if(btagjetsCSVDownVect.size()>1)std::sort(btagjetsCSVDownVect.begin(),btagjetsCSVDownVect.end(),ComparePairsbyPt);

  if(njetspt20UpVect.size()>1)std::sort(njetspt20UpVect.begin(),njetspt20UpVect.end(),ComparePairsbyPt);

  if(btagjetsUpVect.size()>1)std::sort(btagjetsUpVect.begin(),btagjetsUpVect.end(),ComparePairsbyPt);
  if(btagjetsCSVUpVect.size()>1)std::sort(btagjetsCSVUpVect.begin(),btagjetsCSVUpVect.end(),ComparePairsbyPt);


  int _njetingap=0;
  int _njetingap20=0;
  if(njetspt20Vect.size()>=2)
    {

      mjj.push_back((njetspt20Vect.at(0).p4()+njetspt20Vect.at(1).p4()).M());
      jdeta.push_back(abs(njetspt20Vect.at(0).eta()-njetspt20Vect.at(1).eta()));

      for(unsigned int i=0;i<njetspt20Vect.size();i++)
	{

	  if(deltaR(njetspt20Vect.at(i).p4(),njetspt20Vect.at(0).p4())>0.00000000001 && deltaR(njetspt20Vect.at(i).p4(),njetspt20Vect.at(1).p4())>0.00000000001 && njetspt20Vect.at(i).eta()>njetspt20Vect.at(1).eta() &&  njetspt20Vect.at(i).eta() < njetspt20Vect.at(0).eta()) {

	    if(njetspt20Vect.at(i).pt()>30 && njetspt20GapVect.at(i)==1)_njetingap++;

	    if(njetspt20GapVect.at(i)==1)_njetingap20++;

	  }
	}
      jdphi.push_back(deltaPhi(njetspt20Vect.at(0).phi(),njetspt20Vect.at(1).phi()));
      dijetpt.push_back((njetspt20Vect.at(0).p4()+njetspt20Vect.at(1).p4()).pt());
      dijetphi.push_back((njetspt20Vect.at(0).p4()+njetspt20Vect.at(1).p4()).phi());
      //ptvis

    }
  else
    {
      mjj.push_back(-99);
      jdeta.push_back(-99);
      jdphi.push_back(-99);
      dijetpt.push_back(-99);
      dijetphi.push_back(-99);
    }

  int _njetingapDown=0;
  int _njetingap20Down=0;
  if(njetspt20DownVect.size()>=2)
    {

      mjjDown.push_back((njetspt20DownVect.at(0).p4()+njetspt20DownVect.at(1).p4()).M());
      jdetaDown.push_back(abs(njetspt20DownVect.at(0).eta()-njetspt20DownVect.at(1).eta()));

      for(unsigned int i=0;i<njetspt20DownVect.size();i++)
	{

	  if(deltaR(njetspt20DownVect.at(i).p4(),njetspt20DownVect.at(0).p4())>0.00000000001 && deltaR(njetspt20DownVect.at(i).p4(),njetspt20DownVect.at(1).p4())>0.00000000001 && njetspt20DownVect.at(i).eta()>njetspt20DownVect.at(1).eta() &&  njetspt20DownVect.at(i).eta() < njetspt20DownVect.at(0).eta()) {

	    if(njetspt20DownVect.at(i).pt()>30 && njetspt20GapDownVect.at(i)==1)_njetingapDown++;

	    if(njetspt20GapDownVect.at(i)==1)_njetingap20Down++;

	  }
	}
      jdphiDown.push_back(deltaPhi(njetspt20DownVect.at(0).phi(),njetspt20DownVect.at(1).phi()));
      dijetptDown.push_back((njetspt20DownVect.at(0).p4()+njetspt20DownVect.at(1).p4()).pt());
      dijetphiDown.push_back((njetspt20DownVect.at(0).p4()+njetspt20DownVect.at(1).p4()).phi());
    }
  else
    {
      mjjDown.push_back(-99);
      jdetaDown.push_back(-99);
      jdphiDown.push_back(-99);
      dijetptDown.push_back(-99);
      dijetphiDown.push_back(-99);
    }


  int _njetingapUp=0;
  int _njetingap20Up=0;
  if(njetspt20UpVect.size()>=2)
    {

      mjjUp.push_back((njetspt20UpVect.at(0).p4()+njetspt20UpVect.at(1).p4()).M());
      jdetaUp.push_back(abs(njetspt20UpVect.at(0).eta()-njetspt20UpVect.at(1).eta()));

      for(unsigned int i=0;i<njetspt20UpVect.size();i++)
	{

	  if(deltaR(njetspt20UpVect.at(i).p4(),njetspt20UpVect.at(0).p4())>0.00000000001 && deltaR(njetspt20UpVect.at(i).p4(),njetspt20UpVect.at(1).p4())>0.00000000001 && njetspt20UpVect.at(i).eta()>njetspt20UpVect.at(1).eta() &&  njetspt20UpVect.at(i).eta() < njetspt20UpVect.at(0).eta()) {

	    if(njetspt20UpVect.at(i).pt()>30 && njetspt20GapUpVect.at(i)==1)_njetingapUp++;

	    if(njetspt20GapUpVect.at(i)==1)_njetingap20Up++;

	  }
	}
      jdphiUp.push_back(deltaPhi(njetspt20UpVect.at(0).phi(),njetspt20UpVect.at(1).phi()));
      dijetptUp.push_back((njetspt20UpVect.at(0).p4()+njetspt20UpVect.at(1).p4()).pt());
      dijetphiUp.push_back((njetspt20UpVect.at(0).p4()+njetspt20UpVect.at(1).p4()).phi());
    }
  else
    {
      mjjUp.push_back(-99);
      jdetaUp.push_back(-99);
      jdphiUp.push_back(-99);
      dijetptUp.push_back(-99);
      dijetphiUp.push_back(-99);
    }


  njetingap.push_back(_njetingap);
  njetingap20.push_back(_njetingap20);

  njetspt20.push_back(njetspt20Vect.size());
  int _njets=0;
  for(int i=0;i<njetspt20.at(ipair);i++){
    if(njetspt20Vect.at(i).pt()>30)_njets++;
  }

  njetingapDown.push_back(_njetingapDown);
  njetingap20Down.push_back(_njetingap20Down);

  njetspt20Down.push_back(njetspt20DownVect.size());
  int _njetsDown=0;
  for(int i=0;i<njetspt20Down.at(ipair);i++){
    if(njetspt20DownVect.at(i).pt()>30)_njetsDown++;
  }

  njetingapUp.push_back(_njetingapUp);
  njetingap20Up.push_back(_njetingap20Up);

  njetspt20Up.push_back(njetspt20UpVect.size());
  int _njetsUp=0;
  for(int i=0;i<njetspt20Up.at(ipair);i++){
    if(njetspt20UpVect.at(i).pt()>30)_njetsUp++;
  }


  if(njetspt20Vect.size()>0)
    {
      jpt_1.push_back(njetspt20Vect.at(0).pt());
      jeta_1.push_back(njetspt20Vect.at(0).eta());
      jphi_1.push_back(njetspt20Vect.at(0).phi());
      //jcsv_1
    }

  if(njetspt20Vect.size()>1)
    {
      jpt_2.push_back(njetspt20Vect.at(1).pt());
      jeta_2.push_back(njetspt20Vect.at(1).eta());
      jphi_2.push_back(njetspt20Vect.at(1).phi());
      //jcsv_2
    }

  if(njetspt20DownVect.size()>0)
    {
      jptDown_1.push_back(njetspt20DownVect.at(0).pt());
      jetaDown_1.push_back(njetspt20DownVect.at(0).eta());
      jphiDown_1.push_back(njetspt20DownVect.at(0).phi());
    }

  if(njetspt20DownVect.size()>1)
    {
      jptDown_2.push_back(njetspt20DownVect.at(1).pt());
      jetaDown_2.push_back(njetspt20DownVect.at(1).eta());
      jphiDown_2.push_back(njetspt20DownVect.at(1).phi());
    }

  if(njetspt20UpVect.size()>0)
    {
      jptUp_1.push_back(njetspt20UpVect.at(0).pt());
      jetaUp_1.push_back(njetspt20UpVect.at(0).eta());
      jphiUp_1.push_back(njetspt20UpVect.at(0).phi());
    }

  if(njetspt20UpVect.size()>1)
    {
      jptUp_2.push_back(njetspt20UpVect.at(1).pt());
      jetaUp_2.push_back(njetspt20UpVect.at(1).eta());
      jphiUp_2.push_back(njetspt20UpVect.at(1).phi());
    }

  nbtag.push_back(btagjetsVect.size());
  if(btagjetsVect.size()>0)
    {
      bpt_1.push_back(btagjetsVect.at(0).pt());
      beta_1.push_back(btagjetsVect.at(0).eta());
      bphi_1.push_back(btagjetsVect.at(0).phi());
    }
  if(btagjetsCSVVect.size()>0)
    {
      bcsv_1.push_back(btagjetsCSVVect.at(0).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVVect.at(0).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }
  if(btagjetsVect.size()>1)
    {
      bpt_2.push_back(btagjetsVect.at(1).pt());
      beta_2.push_back(btagjetsVect.at(1).eta());
      bphi_2.push_back(btagjetsVect.at(1).phi());
    }
  if(btagjetsCSVVect.size()>1)
    {
      bcsv_2.push_back(btagjetsCSVVect.at(1).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVVect.at(1).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }

  nbtagDown.push_back(btagjetsDownVect.size());
  if(btagjetsDownVect.size()>0)
    {
      bptDown_1.push_back(btagjetsDownVect.at(0).pt());
      betaDown_1.push_back(btagjetsDownVect.at(0).eta());
      bphiDown_1.push_back(btagjetsDownVect.at(0).phi());
    }
  if(btagjetsCSVDownVect.size()>0)
    {
      bcsvDown_1.push_back(btagjetsCSVDownVect.at(0).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVDownVect.at(0).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }
  if(btagjetsDownVect.size()>1)
    {
      bptDown_2.push_back(btagjetsDownVect.at(1).pt());
      betaDown_2.push_back(btagjetsDownVect.at(1).eta());
      bphiDown_2.push_back(btagjetsDownVect.at(1).phi());
    }
  if(btagjetsCSVDownVect.size()>1)
    {
      bcsvDown_2.push_back(btagjetsCSVDownVect.at(1).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVDownVect.at(1).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }


  nbtagUp.push_back(btagjetsUpVect.size());
  if(btagjetsUpVect.size()>0)
    {
      bptUp_1.push_back(btagjetsUpVect.at(0).pt());
      betaUp_1.push_back(btagjetsUpVect.at(0).eta());
      bphiUp_1.push_back(btagjetsUpVect.at(0).phi());
    }
  if(btagjetsCSVUpVect.size()>0)
    {
      bcsvUp_1.push_back(btagjetsCSVUpVect.at(0).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVUpVect.at(0).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }
  if(btagjetsUpVect.size()>1)
    {
      bptUp_2.push_back(btagjetsUpVect.at(1).pt());
      betaUp_2.push_back(btagjetsUpVect.at(1).eta());
      bphiUp_2.push_back(btagjetsUpVect.at(1).phi());
    }
  if(btagjetsCSVUpVect.size()>1)
    {
      bcsvUp_2.push_back(btagjetsCSVUpVect.at(1).bDiscriminator("pfDeepCSVJetTags:probb") + btagjetsCSVUpVect.at(1).bDiscriminator("pfDeepCSVJetTags:probbb"));
    }

  if(theisMC && !IsEmbed)
    {
      edm::Handle<edm::View<reco::GenJet>> genJetHandle;
      //event.getByLabel ("slimmedGenJets", genJetHandle);
      event.getByToken (theGenJetTag, genJetHandle);
      std::vector <const reco::GenJet*> ptrGenJets;
      unsigned int genJetSize = genJetHandle->size();
      for (unsigned int igj = 0; igj < genJetSize; igj++)
	{
	  ptrGenJets.push_back ( &((*genJetHandle)[igj]) );
	}

      // now gather associated gen jets and save their address
      // relies on the fact that gen jets are filled in the same order as they appear in the collection
      std::vector <const reco::GenJet*>::iterator it;
      for(edm::View<pat::Jet>::const_iterator ijet = jets->begin(); ijet!=jets->end();++ijet)
	{
	  float jecFactor = ijet->jecFactor("Uncorrected");
	  float jetRawPt = jecFactor * ijet->pt();
	  if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;
	  const reco::GenJet * thisGenJet =  ijet->genJet ();
	  int genindex = -1;
	  it = std::find(ptrGenJets.begin(), ptrGenJets.end(), thisGenJet);
	  if (it != ptrGenJets.end())
	    genindex = std::distance (ptrGenJets.begin(), it);
	  _jets_genjetIndex.push_back(genindex);
	}
      for(edm::View<pat::Jet>::const_iterator ijet = jetsDown->begin(); ijet!=jetsDown->end();++ijet)
	{
	  float jecFactor = ijet->jecFactor("Uncorrected");
	  float jetRawPt = jecFactor * ijet->pt();
	  if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;
	  const reco::GenJet * thisGenJet =  ijet->genJet ();
	  int genindex = -1;
	  it = std::find(ptrGenJets.begin(), ptrGenJets.end(), thisGenJet);
	  if (it != ptrGenJets.end())
	    genindex = std::distance (ptrGenJets.begin(), it);
	  _jetsDown_genjetIndex.push_back(genindex);
	}

      for(edm::View<pat::Jet>::const_iterator ijet = jetsUp->begin(); ijet!=jetsUp->end();++ijet)
	{
	  float jecFactor = ijet->jecFactor("Uncorrected");
	  float jetRawPt = jecFactor * ijet->pt();
	  if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;
	  const reco::GenJet * thisGenJet =  ijet->genJet ();
	  int genindex = -1;
	  it = std::find(ptrGenJets.begin(), ptrGenJets.end(), thisGenJet);
	  if (it != ptrGenJets.end())
	    genindex = std::distance (ptrGenJets.begin(), it);
	  _jetsUp_genjetIndex.push_back(genindex);
	}
    }
  njetsDown.push_back(_njetsDown);
  njetsUp.push_back(_njetsUp);
  return _njets;
}

void HTauTauNtuplizer::FillFatJet(const edm::View<pat::Jet>* fatjets, const edm::Event&)
{
  for(edm::View<pat::Jet>::const_iterator ijet = fatjets->begin(); ijet!=fatjets->end();++ijet)
    {
      float jecFactor = ijet->jecFactor("Uncorrected");
      float jetRawPt = jecFactor * ijet->pt();

      if(theYear==2017 && jetRawPt<50 &&  fabs(ijet->eta()) < 3.139 && fabs(ijet->eta()) > 2.65) continue;

      _ak8jets_px.push_back( (float) ijet->px());
      _ak8jets_py.push_back( (float) ijet->py());
      _ak8jets_pz.push_back( (float) ijet->pz());
      _ak8jets_e.push_back( (float) ijet->energy());
      _ak8jets_SoftDropMass.push_back (ijet->hasUserFloat("ak8PFJetsPuppiSoftDropMass") ? ijet->userFloat("ak8PFJetsPuppiSoftDropMass") : -999 );
      _ak8jets_PrunedMass.push_back   (ijet->hasUserFloat("ak8PFJetsPuppiPrunedMass")   ? ijet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   : -999 );
      _ak8jets_TrimmedMass.push_back  (ijet->hasUserFloat("ak8PFJetsPuppiTrimmedMass")  ? ijet->userFloat("ak8PFJetsPuppiTrimmedMass")  : -999 );
      _ak8jets_FilteredMass.push_back (ijet->hasUserFloat("ak8PFJetsPuppiFilteredMass") ? ijet->userFloat("ak8PFJetsPuppiFilteredMass") : -999 );
      _ak8jets_tau1.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau1")    ? ijet->userFloat("NjettinessAK8Puppi:tau1")    : -999 );
      _ak8jets_tau2.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau2")    ? ijet->userFloat("NjettinessAK8Puppi:tau2")    : -999 );
      _ak8jets_tau3.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau3")    ? ijet->userFloat("NjettinessAK8Puppi:tau3")    : -999 );
      _ak8jets_tau4.push_back         (ijet->hasUserFloat("NjettinessAK8Puppi:tau4")    ? ijet->userFloat("NjettinessAK8Puppi:tau4")    : -999 );
      _ak8jets_CSV.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      _ak8jets_deepCSV_probb.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probb"));
      _ak8jets_deepCSV_probbb.push_back(ijet->bDiscriminator("pfDeepCSVJetTags:probbb"));
      _ak8jets_deepFlavor_probb.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probb"));
      _ak8jets_deepFlavor_probbb.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
      _ak8jets_deepFlavor_problepb.push_back(ijet->bDiscriminator("pfDeepFlavourJetTags:problepb"));

      // store subjets for soft drop
      int nsubj = 0;
      //if (ijet->hasSubjets("SoftDrop"))
      if (ijet->hasSubjets("SoftDropPuppi"))
	{
	  //pat::JetPtrCollection const & subj = ijet->subjets("SoftDrop");
	  pat::JetPtrCollection const & subj = ijet->subjets("SoftDropPuppi");
	  // cout << "============= IJET " << ijet - fatjets->begin() << " ============= " << endl;
	  for (auto isubj = subj.begin(); isubj!=subj.end(); isubj++)
	    {
	      nsubj += 1;
	      _subjets_px.push_back((*isubj)->px());
	      _subjets_py.push_back((*isubj)->py());
	      _subjets_pz.push_back((*isubj)->pz());
	      _subjets_e.push_back((*isubj)->energy());
	      _subjets_CSV.push_back((*isubj)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	      _subjets_deepCSV_probb.push_back((*isubj)->bDiscriminator("pfDeepCSVJetTags:probb"));
	      _subjets_deepCSV_probbb.push_back((*isubj)->bDiscriminator("pfDeepCSVJetTags:probbb"));
	      _subjets_deepFlavor_probb.push_back((*isubj)->bDiscriminator("pfDeepFlavourJetTags:probb"));
	      _subjets_deepFlavor_probbb.push_back((*isubj)->bDiscriminator("pfDeepFlavourJetTags:probbb"));
	      _subjets_deepFlavor_problepb.push_back((*isubj)->bDiscriminator("pfDeepFlavourJetTags:problepb"));
	      _subjets_ak8MotherIdx.push_back(ijet - fatjets->begin()); // idx of fatjet in fatjet vector
	      // cout << " * " << (*isubj)->pt() << " " << isubj - subj.begin() << endl;
	    }
	}
      _ak8jets_nsubjets.push_back(nsubj);
    }
}


//Fill all leptons (we keep them all for veto purposes
void HTauTauNtuplizer::FillSoftLeptons(const edm::View<reco::Candidate> *daus,
				       const edm::Event& event, const edm::EventSetup& setup,
				       bool theFSR, const edm::View<pat::Jet> *jets,const BXVector<l1t::Tau>* l1taus ) {

  //cout<<"-----------------------------"<<endl;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  event.getByToken(triggerObjects_, triggerObjects);
  event.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = event.triggerNames(*triggerBits);

  edm::Handle<edm::View<pat::PackedCandidate> >pfCandHandle;
  event.getByToken(thePFCandTag,pfCandHandle);
  const edm::View<pat::PackedCandidate>* pfCands = pfCandHandle.product();

  std::vector<const pat::PackedCandidate *> pfCands_charged;
  std::vector<const pat::PackedCandidate *> pfCands_neutral;
  LeptonIsoHelper::PFIso_particles(pfCands, pfCands_charged, pfCands_neutral);

  edm::Handle<double> rhoHandle_miniRelIso;
  event.getByToken(theRhoMiniRelIsoTag, rhoHandle_miniRelIso);
  float rho_miniRelIso = *rhoHandle_miniRelIso;
  unsigned int ntaus(0);

  int compteurCand=-1;
  for(edm::View<reco::Candidate>::const_iterator daui = daus->begin(); daui!=daus->end();++daui) {
    compteurCand++;
    const reco::Candidate* cand = &(*daui);
    math::XYZTLorentzVector pfour = cand->p4();
    LeptonP4.push_back(pfour);
    math::XYZTLorentzVector chargedP4 = cand->p4();
    math::XYZTLorentzVector neutralP4;

    TLorentzVector pfourTauUp;
    TLorentzVector pfourTauDown;
    TLorentzVector pfourEleUp;
    TLorentzVector pfourEleDown;

    bool existTESshift = userdatahelpers::hasUserInt(cand,"isTESShifted"); // simply to check if the userfloat exists
    bool existEESshift = userdatahelpers::hasUserInt(cand,"isEESShifted"); // simply to check if the userfloat exists
    bool isTauMatched = userdatahelpers::hasUserInt(cand,"isTauMatched"); // check if it is a tau

    int (hasTES) = ( existTESshift ? userdatahelpers::getUserInt(cand,"isTESShifted") : false) ;   // actual check of the value of the userfloat
    int (hasEES) = ( existEESshift ? userdatahelpers::getUserInt(cand,"isEESShifted") : false) ;   // actual check of the value of the userfloat

    if(hasTES)
      {
	pfourTauUp.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_TauUp"),userdatahelpers::getUserFloat(cand,"py_TauUp"),userdatahelpers::getUserFloat(cand,"pz_TauUp"),userdatahelpers::getUserFloat(cand,"e_TauUp"));
	pfourTauDown.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_TauDown"),userdatahelpers::getUserFloat(cand,"py_TauDown"),userdatahelpers::getUserFloat(cand,"pz_TauDown"),userdatahelpers::getUserFloat(cand,"e_TauDown"));
      }
    else
      {
	pfourTauUp.SetPxPyPzE(-999.,-999.,-999.,-999.);
	pfourTauDown.SetPxPyPzE(-999.,-999.,-999.,-999.);
      }
    if(hasEES)
      {
	pfourEleUp.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_EleUp"),userdatahelpers::getUserFloat(cand,"py_EleUp"),userdatahelpers::getUserFloat(cand,"pz_EleUp"),userdatahelpers::getUserFloat(cand,"e_EleUp"));
	pfourEleDown.SetPxPyPzE(userdatahelpers::getUserFloat(cand,"px_EleDown"),userdatahelpers::getUserFloat(cand,"py_EleDown"),userdatahelpers::getUserFloat(cand,"pz_EleDown"),userdatahelpers::getUserFloat(cand,"e_EleDown"));
      }
    else
      {
	pfourEleUp.SetPxPyPzE(-999.,-999.,-999.,-999.);
	pfourEleDown.SetPxPyPzE(-999.,-999.,-999.,-999.);
      }

    if(theFSR){
      const pat::PFParticle* fsr=0;
      double maxPT=-1;
      const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(cand);
      if (gammas!=0) {
	for (PhotonPtrVector::const_iterator g = gammas->begin();g!= gammas->end(); ++g) {
	  //const pat::Photon* gamma = g->get();
	  const pat::PFParticle* gamma = g->get();
	  double pt = gamma->pt();
	  if (pt>maxPT) {
	    maxPT  = pt;
	    fsr = gamma;
	  }
	}
      }
      if (fsr!=0) {
	pfour+=fsr->p4();
      }
    }
    //    std::cout<<"  Px   " <<(float) pfour.X() <<std::endl;
    _daughters_px.push_back( (float) pfour.X());
    _daughters_py.push_back( (float) pfour.Y());
    _daughters_pz.push_back( (float) pfour.Z());
    _daughters_e.push_back( (float) pfour.T());

    _daughters_hasTES.push_back( (hasTES ? 1 : 0) );
    _daughters_px_TauUp.push_back((float)pfourTauUp.Px());
    _daughters_py_TauUp.push_back((float)pfourTauUp.Py());
    _daughters_pz_TauUp.push_back((float)pfourTauUp.Pz());
    _daughters_e_TauUp.push_back((float)pfourTauUp.E());
    _daughters_px_TauDown.push_back((float)pfourTauDown.Px());
    _daughters_py_TauDown.push_back((float)pfourTauDown.Py());
    _daughters_pz_TauDown.push_back((float)pfourTauDown.Pz());
    _daughters_e_TauDown.push_back((float)pfourTauDown.E());

    _daughters_hasEES.push_back( (hasEES ? 1 : 0) );
    _daughters_px_EleUp.push_back((float)pfourEleUp.Px());
    _daughters_py_EleUp.push_back((float)pfourEleUp.Py());
    _daughters_pz_EleUp.push_back((float)pfourEleUp.Pz());
    _daughters_e_EleUp.push_back((float)pfourEleUp.E());
    _daughters_px_EleDown.push_back((float)pfourEleDown.Px());
    _daughters_py_EleDown.push_back((float)pfourEleDown.Py());
    _daughters_pz_EleDown.push_back((float)pfourEleDown.Pz());
    _daughters_e_EleDown.push_back((float)pfourEleDown.E());
    _daughters_isTauMatched.push_back( (isTauMatched ? 1 : 0) );


    // gen info

    if (theisMC || IsEmbed)
      {
	int iMatched = GetMatchedGen (cand, event);
	//cout<<"iMatched:"<<iMatched<<endl;
	_daughters_genindex.push_back(iMatched);
      }
    //math::XYZTLorentzVector pfour(userdatahelpers::getUserFloat(cand,"genPx"),userdatahelpers::getUserFloat(cand,"genPy"),userdatahelpers::getUserFloat(cand,"genPz"),userdatahelpers::getUserFloat(cand,"genE"));
    //if(theisMC)_genDaughters.push_back(userdatahelpers::getUserFloat(cand,"fromH"));

    _softLeptons.push_back(cand);//This is needed also for FindCandIndex

    _daughters_vx.push_back(cand->vx());
    _daughters_vy.push_back(cand->vy());
    _daughters_vz.push_back(cand->vz());

    _pdgdau.push_back(cand->pdgId());
    _combreliso.push_back(userdatahelpers::getUserFloat(cand,"combRelIsoPF"));
    _combreliso03.push_back( userdatahelpers::hasUserFloat(cand,"combRelIsoPF03") ? userdatahelpers::getUserFloat(cand,"combRelIsoPF03") : -1 );
    _dxy.push_back(userdatahelpers::getUserFloat(cand,"dxy"));
    _dz.push_back(userdatahelpers::getUserFloat(cand,"dz"));

    int type = ParticleType::TAU;
    if(cand->isMuon()) type = ParticleType::MUON;
    else if(cand->isElectron()) type = ParticleType::ELECTRON;
    _particleType.push_back(type);
    boost::hash<const reco::Candidate*> hasher;
    _LeptonHash.push_back(hasher(cand));

    //Find closest jet for lepton MVA
    float dRmin_cand_jet = 0.4;
    pat::Jet closest_jet;

    for(edm::View<pat::Jet>::const_iterator jeti = jets->begin(); jeti!=jets->end();++jeti){

      float jecFactor = jeti->jecFactor("Uncorrected");
      float jetRawPt = jecFactor * jeti->pt();

      if(theYear==2017 && jetRawPt<50 &&  fabs(jeti->eta()) < 3.139 && fabs(jeti->eta()) > 2.65) continue;

      float dR_cand_jet = deltaR(*cand,*jeti);
      if(dR_cand_jet<dRmin_cand_jet){
	closest_jet = (*jeti);
	dRmin_cand_jet = dR_cand_jet;
      }
    }
    int muIDflag = 0;
    bool iseleLoose = false;
    bool isele80=false;
    bool isele90=false;
    bool iselenoisoLoose = false;
    bool iselenoiso80=false;
    bool iselenoiso90=false;
    bool iselecutbased=false;
    float elemva=-2;
    float elemva_HZZ=-2;
    bool isconversionveto=false;
    int elemissinghits = 999;
    bool iselechargeconsistent=false;

    int decay=-1, MVADM2017v1=-99, genmatch = -1;
    float ieta=-1,full5x5_ieta=-1,hOverE=-1,etasuperatvtx=-1,phisuperatvtx=-1,IoEmIoP=-999.,IoEmIoP_ttH=-999.,depositTracker=-1,depositEcal=-1,depositHcal=-1;//,SCeta=-999.;
    int decayModeFindingOldDMs=-1, decayModeFindingNewDMs=-1; // tau 13 TeV ID
    float byCombinedIsolationDeltaBetaCorrRaw3Hits=-1., chargedIsoPtSum=-1., neutralIsoPtSum=-1., puCorrPtSum=-1.; // tau 13 TeV RAW iso info
    int numChargedParticlesSignalCone=-1, numNeutralHadronsSignalCone=-1, numPhotonsSignalCone=-1, numParticlesSignalCone=-1, numChargedParticlesIsoCone=-1, numNeutralHadronsIsoCone=-1, numPhotonsIsoCone=-1, numParticlesIsoCone=-1;
    float leadChargedParticlePt=-1., trackRefPt=-1.;
    int typeOfMuon=0;
    float byIsolationMVArun2v1DBoldDMwLTraw=-1, byIsolationMVArun2017v2DBoldDMwLTraw2017=-1, byIsolationMVArun2017v1DBoldDMwLTraw2017=-1, byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017=-1; //FRA
    int  byVVLooseIsolationMVArun2017v2DBoldDMwLT2017=-1; //FRA
    float byDeepTau2017v2p1VSjetraw=-1, byDeepTau2017v2p1VSeraw=-1, byDeepTau2017v2p1VSmuraw=-1;
    Long64_t tauIDflag = 0;
    float footprintCorrection, neutralIsoPtSumWeight, photonPtSumOutsideSignalCone;
    std::vector<double> SVPos;
    std::vector<double> SVCov;

    double GEOMFlightLenght(-999);
    double GEOMFlightLenghtSignificance(-999);
    double iFLSign(-999.);


    std::vector<double> SVChi2NDof;
    std::vector<std::vector<double>  >PionsP4;
    std::vector<double> PionsCharge;
    std::vector<std::vector<double>  >RefitPionsP4;
    std::vector<double> RefitPionsCharge;
    float iPFTauTrack_deltaR=-999;
    std::vector<double> iPFTauTrackLV;
    float ia1_M=-999.;
    float ia1_B=-999.;
    int ia1_pdgid=-999;
    int ia1_Charge=-999;
    std::vector<double> ia1_cov;
    std::vector<double> ia1_par;



    float dxy_innerTrack = -1., dz_innerTrack = -1., sip = -1., error_trackpt=-1.;
    //    int jetNDauChargedMVASel = -1;
    float miniRelIsoCharged = -1., miniRelIsoNeutral = -1.;;
    float jetPtRel = -1., jetPtRatio = -1.;
    float jetBTagCSV=-1., jetBTagDeepCSV=-1., jetBTagDeepFlavor=-1.;
    //    float lepMVA_mvaId = -1.;
    float iMuon_M=-999.;
    float iMuon_B=-999.;
    int iMuon_pdgid=-999;
    int iMuon_trackCharge=-999;
    std::vector<double> iMuon_cov;
    std::vector<double> iMuon_par;

    float iPFTau_Track_M=-999;
    float iPFTau_Track_B=-999;
    float iPFTau_Track_pdgid=-999;
    float iPFTau_Track_trackCharge=-999;

    std::vector<double> iPFTau_Track_cov;
    std::vector<double> iPFTau_Track_par;
    //
    GlobalPoint aPVPoint(_pv_x, _pv_y, _pv_z);
    GlobalPoint aPVRefitPoint(_pvRefit_x, _pvRefit_y, _pvRefit_z);
    GlobalPoint aPVGenPoint(_pvGen_x, _pvGen_y, _pvGen_z);
    TVector3 pcaPV = getPCA(event, setup, cand->bestTrack(), aPVPoint);
    TVector3 pcaRefitPV = getPCA(event, setup, cand->bestTrack(), aPVRefitPoint);
    TVector3 pcaGenPV;
    std::vector<float> pcaRefitPVBS_vecx;
    std::vector<float> pcaRefitPVBS_vecy;
    std::vector<float> pcaRefitPVBS_vecz;
    if(theisMC) pcaGenPV = getPCA(event, setup, cand->bestTrack(), aPVGenPoint);
    if(type==ParticleType::MUON){
      for(unsigned int i=0; i<_RefitPVBS_x.size(); i++){
        GlobalPoint aPVBSPoint(_RefitPVBS_x.at(i), _RefitPVBS_y.at(i), _RefitPVBS_z.at(i));
        TVector3 pcaRefitPVBS = getPCA(event, setup, cand->bestTrack(), aPVBSPoint);
	pcaRefitPVBS_vecx.push_back(pcaRefitPVBS.X());
        pcaRefitPVBS_vecy.push_back(pcaRefitPVBS.Y());
        pcaRefitPVBS_vecz.push_back(pcaRefitPVBS.Z());
      }
      //_MVADM2016v1.push_back(-99);
      //_MVADM2017v1.push_back(-99);
      muIDflag=userdatahelpers::getUserInt(cand,"muonID");
      //discr = (float) muIDflag; // not really needed, will use the muonID branch in ntuples...

      if(userdatahelpers::getUserFloat(cand,"isPFMuon"))typeOfMuon |= 1 << 0;
      if(userdatahelpers::getUserFloat(cand,"isGlobalMuon"))typeOfMuon |= 1 << 1;
      if(userdatahelpers::getUserFloat(cand,"isTrackerMuon"))typeOfMuon |= 1 << 2;
      depositTracker=userdatahelpers::getUserFloat(cand,"DepositR03TrackerOfficial");
      depositEcal=userdatahelpers::getUserFloat(cand,"DepositR03ECal");
      depositHcal=userdatahelpers::getUserFloat(cand,"DepositR03Hcal");

      dxy_innerTrack = userdatahelpers::getUserFloat(cand,"dxy_innerTrack");
      dz_innerTrack = userdatahelpers::getUserFloat(cand,"dz_innerTrack");
      error_trackpt = userdatahelpers::getUserFloat(cand,"rel_error_trackpt");
      sip = userdatahelpers::getUserFloat(cand,"SIP");

      //      jetNDauChargedMVASel= LeptonIsoHelper::jetNDauChargedMVASel(cand, closest_jet);
      std::pair<float,float> miniRelIso = LeptonIsoHelper::miniRelIso_ChargedNeutral(cand, pfCands_charged, pfCands_neutral, rho_miniRelIso, theYear);

      miniRelIsoCharged = miniRelIso.first;
      miniRelIsoNeutral = miniRelIso.second;
      jetPtRel = LeptonIsoHelper::jetPtRel(*cand, closest_jet,theJECName);
      jetPtRatio = LeptonIsoHelper::jetPtRatio(*cand, closest_jet,theJECName);
      jetBTagDeepCSV = closest_jet.bDiscriminator("pfDeepCSVJetTags:probb") + closest_jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      jetBTagDeepFlavor = closest_jet.bDiscriminator("pfDeepFlavourJetTags:probb") + closest_jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + closest_jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
      jetBTagCSV = closest_jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

      //      lepMVA_mvaId  = userdatahelpers::getUserFloat(cand,"segmentCompatibility");

      const pat::Muon *patmuon  = dynamic_cast<const pat::Muon*>(cand);

      iMuon_M=userdatahelpers::getUserFloat(cand,"Muon_M");
      iMuon_B=userdatahelpers::getUserFloat(cand,"Muon_B");
      iMuon_pdgid=userdatahelpers::getUserInt(cand,"Muon_pdgid");
      iMuon_trackCharge=userdatahelpers::getUserInt(cand,"Muon_trackCharge");


      if(patmuon->hasUserData("Muon_cov")){
	for(unsigned int i =0; i < patmuon->userData<std::vector<double>  >("Muon_cov")->size(); i++) {iMuon_cov.push_back(patmuon->userData<std::vector<double>  >("Muon_cov")->at(i)); }
      }
      if(patmuon->hasUserData("Muon_par")){
	for(unsigned int i =0; i < patmuon->userData<std::vector<double> >("Muon_par")->size(); i++){iMuon_par.push_back(patmuon->userData<std::vector<double> >("Muon_par")->at(i));  }
      }


    }else if(type==ParticleType::ELECTRON){

      ieta=userdatahelpers::getUserFloat(cand,"sigmaIetaIeta");
      full5x5_ieta=userdatahelpers::getUserFloat(cand,"full5x5_sigmaIetaIeta");
      hOverE=userdatahelpers::getUserFloat(cand,"hOverE");
      etasuperatvtx=userdatahelpers::getUserFloat(cand,"deltaEtaSuperClusterTrackAtVtx");
      phisuperatvtx=userdatahelpers::getUserFloat(cand,"deltaPhiSuperClusterTrackAtVtx");
      IoEmIoP=userdatahelpers::getUserFloat(cand,"IoEmIoP");
      IoEmIoP_ttH=userdatahelpers::getUserFloat(cand,"IoEmIoP_ttH");

      if(userdatahelpers::getUserFloat(cand,"isEleIDLoose") == 1) iseleLoose=true;
      if(userdatahelpers::getUserFloat(cand,"isEleID80") == 1) isele80=true;
      if(userdatahelpers::getUserFloat(cand,"isEleID90") == 1) isele90=true;
      if(userdatahelpers::getUserFloat(cand,"isEleNoIsoIDLoose") == 1) iselenoisoLoose=true;
      if(userdatahelpers::getUserFloat(cand,"isEleNoIsoID80") == 1) iselenoiso80=true;
      if(userdatahelpers::getUserFloat(cand,"isEleNoIsoID90") == 1) iselenoiso90=true;
      if(userdatahelpers::getUserFloat(cand,"isEleCutBased") == 1) iselecutbased=true;
      elemva=(userdatahelpers::getUserFloat(cand,"mvaValue_Iso"));
      elemva_HZZ=(userdatahelpers::getUserFloat(cand,"mvaValue_HZZ"));

      if(userdatahelpers::getUserInt(cand,"isConversionVeto") == 1)isconversionveto=true;
      error_trackpt = userdatahelpers::getUserFloat(cand,"rel_error_trackpt");
      elemissinghits = userdatahelpers::getUserInt(cand,"missingHit");
      if((userdatahelpers::getUserInt(cand,"isGsfCtfScPixChargeConsistent") + userdatahelpers::getUserInt(cand,"isGsfScPixChargeConsistent"))>1)iselechargeconsistent=true;
      sip = userdatahelpers::getUserFloat(cand,"SIP");
      std::pair<float,float> miniRelIso = LeptonIsoHelper::miniRelIso_ChargedNeutral(cand, pfCands_charged, pfCands_neutral, rho_miniRelIso, theYear);
      miniRelIsoCharged = miniRelIso.first;
      miniRelIsoNeutral = miniRelIso.second;
      jetPtRel = LeptonIsoHelper::jetPtRel(*cand, closest_jet,theJECName);
      jetPtRatio = LeptonIsoHelper::jetPtRatio(*cand, closest_jet,theJECName);
      jetBTagCSV = closest_jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      jetBTagDeepCSV = closest_jet.bDiscriminator("pfDeepCSVJetTags:probb") + closest_jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      jetBTagDeepFlavor = closest_jet.bDiscriminator("pfDeepFlavourJetTags:probb") + closest_jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + closest_jet.bDiscriminator("pfDeepFlavourJetTags:problepb");

    }else if(type==ParticleType::TAU){
      decay = userdatahelpers::getUserFloat(cand,"decayMode");
      genmatch = userdatahelpers::getUserFloat(cand,"genmatch");
      MVADM2017v1= userdatahelpers::getUserFloat(cand,"MVADM2017v1");
      decayModeFindingOldDMs = userdatahelpers::getUserInt (cand, "decayModeFinding");
      decayModeFindingNewDMs = userdatahelpers::getUserInt (cand, "decayModeFindingNewDMs");
      for (uint itau =0; itau<ntauIds; itau++){
	int id = userdatahelpers::getUserInt (cand,  tauIDStrings[itau]);
	if(id>0){
	  tauIDflag |= (Long64_t(1) << itau);
	  hTauIDs->Fill(id);
	}
      }
      //      if(decayModeFindingOldDMs>0.5 && (tauIDflag & (1<<3) == (1<<3)) )   ntausc++;
      //   std::cout<<" tauIDflag    " <<tauIDflag <<  "[passed   " <<ntausc<<std::endl;
      //againstElectronMVA5category = userdatahelpers::getUserFloat (cand, "againstElectronMVA5category");
      footprintCorrection = userdatahelpers::getUserFloat (cand, "footprintCorrection");
      neutralIsoPtSumWeight = userdatahelpers::getUserFloat (cand, "neutralIsoPtSumWeight");
      photonPtSumOutsideSignalCone = userdatahelpers::getUserFloat (cand, "photonPtSumOutsideSignalCone");
      byCombinedIsolationDeltaBetaCorrRaw3Hits = userdatahelpers::getUserFloat (cand, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
      byIsolationMVArun2v1DBoldDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2v1DBoldDMwLTraw");
      byIsolationMVArun2017v2DBoldDMwLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v2DBoldDMwLTraw2017"); //FRA
      byIsolationMVArun2017v1DBoldDMwLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v1DBoldDMwLTraw2017"); //FRA
      byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017"); //FRA
      byIsolationMVArun2v1DBoldDMwLTraw=userdatahelpers::getUserFloat (cand, "byIsolationMVArun2v1DBoldDMwLTraw");
      byVVLooseIsolationMVArun2017v2DBoldDMwLT2017= userdatahelpers::getUserInt (cand, "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"); //FRA
      byDeepTau2017v2p1VSjetraw=userdatahelpers::getUserFloat(cand,"byDeepTau2017v2p1VSjetraw");
      byDeepTau2017v2p1VSeraw=userdatahelpers::getUserFloat(cand,"byDeepTau2017v2p1VSeraw");
      byDeepTau2017v2p1VSmuraw=userdatahelpers::getUserFloat(cand,"byDeepTau2017v2p1VSmuraw");
      chargedIsoPtSum = userdatahelpers::getUserFloat (cand, "chargedIsoPtSum");
      neutralIsoPtSum = userdatahelpers::getUserFloat (cand, "neutralIsoPtSum");
      puCorrPtSum = userdatahelpers::getUserFloat (cand, "puCorrPtSum");
      numChargedParticlesSignalCone = userdatahelpers::getUserInt (cand, "numChargedParticlesSignalCone");
      numNeutralHadronsSignalCone = userdatahelpers::getUserInt (cand, "numNeutralHadronsSignalCone");
      numPhotonsSignalCone = userdatahelpers::getUserInt (cand, "numPhotonsSignalCone");
      numParticlesSignalCone = userdatahelpers::getUserInt (cand, "numParticlesSignalCone");
      numChargedParticlesIsoCone = userdatahelpers::getUserInt (cand, "numChargedParticlesIsoCone");
      numNeutralHadronsIsoCone = userdatahelpers::getUserInt (cand, "numNeutralHadronsIsoCone");
      numPhotonsIsoCone = userdatahelpers::getUserInt (cand, "numPhotonsIsoCone");
      numParticlesIsoCone = userdatahelpers::getUserInt (cand, "numParticlesIsoCone");
      leadChargedParticlePt = userdatahelpers::getUserFloat (cand, "leadChargedParticlePt");
      trackRefPt = userdatahelpers::getUserFloat (cand, "trackRefPt");
      iPFTauTrack_deltaR = userdatahelpers::getUserFloat (cand, "PFTauTrack_deltaR");
      iFLSign = userdatahelpers::getUserFloat (cand, "FLSig");
      ia1_M=userdatahelpers::getUserFloat(cand,"a1_M");
      ia1_B=userdatahelpers::getUserFloat(cand,"a1_B");
      ia1_pdgid=userdatahelpers::getUserInt(cand,"a1_pdgid");
      ia1_Charge=userdatahelpers::getUserInt(cand,"a1_charge");

      GEOMFlightLenght= userdatahelpers::getUserFloat (cand, "GEOMFlightLenght");
      GEOMFlightLenghtSignificance= userdatahelpers::getUserFloat (cand, "GEOMFlightLenghtSignificance");
      const pat::Tau *taon  = dynamic_cast<const pat::Tau*>(cand);
      if(taon){
	pcaPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVPoint);
	pcaRefitPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVRefitPoint);
	if(theisMC) pcaGenPV = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVGenPoint);
        for(unsigned int i=0; i<_RefitPVBS_x.size(); i++){
	  GlobalPoint aPVBSPoint(_RefitPVBS_x.at(i), _RefitPVBS_y.at(i), _RefitPVBS_z.at(i));
	  TVector3 pcaRefitPVBS = getPCA(event, setup, taon->leadChargedHadrCand()->bestTrack(), aPVBSPoint);
          pcaRefitPVBS_vecx.push_back(pcaRefitPVBS.X());
          pcaRefitPVBS_vecy.push_back(pcaRefitPVBS.Y());
          pcaRefitPVBS_vecz.push_back(pcaRefitPVBS.Z());
	}

	reco::CandidatePtrVector chCands = taon->signalChargedHadrCands();
	reco::CandidatePtrVector neCands = taon->signalGammaCands();
	chargedP4 = math::XYZTLorentzVector();
	neutralP4 = math::XYZTLorentzVector();
	for(reco::CandidatePtrVector::const_iterator id=chCands.begin();id!=chCands.end(); ++id) chargedP4 += (*id)->p4();
	for(reco::CandidatePtrVector::const_iterator id=neCands.begin();id!=neCands.end(); ++id) neutralP4 += (*id)->p4();

	iPFTau_Track_M=userdatahelpers::getUserFloat(cand,"TauTrackFiller_M");
	iPFTau_Track_B=userdatahelpers::getUserFloat(cand,"TauTrackFiller_B");
	iPFTau_Track_pdgid=userdatahelpers::getUserInt(cand,"TauTrackFiller_pdgid");
	iPFTau_Track_trackCharge=userdatahelpers::getUserInt(cand,"TauTrackFiller_trackCharge");

	if(taon->hasUserData("TauTrackFiller_cov")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double>  >("TauTrackFiller_cov")->size(); i++) {iPFTau_Track_cov.push_back(taon->userData<std::vector<double>  >("TauTrackFiller_cov")->at(i)); }
	}
	if(taon->hasUserData("TauTrackFiller_par")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("TauTrackFiller_par")->size(); i++){iPFTau_Track_par.push_back(taon->userData<std::vector<double> >("TauTrackFiller_par")->at(i));  }

	}

	if(taon->hasUserData("PFTauTrackLV")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("PFTauTrackLV")->size(); i++){iPFTauTrackLV.push_back(taon->userData<std::vector<double> >("PFTauTrackLV")->at(i));  }
	}
	if(taon->hasUserData("SVPos")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("SVPos")->size(); i++){SVPos.push_back(taon->userData<std::vector<double> >("SVPos")->at(i));  }
	}
	if(taon->hasUserData("SVCov")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("SVCov")->size(); i++){SVCov.push_back(taon->userData<std::vector<double> >("SVCov")->at(i));  }
	}
	if(taon->hasUserData("SVChi2NDofMatchingQual")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("SVChi2NDofMatchingQual")->size(); i++){SVChi2NDof.push_back(taon->userData<std::vector<double> >("SVChi2NDofMatchingQual")->at(i));  }
	}
	if(taon->hasUserData("iPionP4")){
	  for(unsigned int i =0; i < taon->userData<std::vector<std::vector<double> > >("iPionP4")->size(); i++) {PionsP4.push_back(taon->userData<std::vector<std::vector<double> > >("iPionP4")->at(i));
	  }
	}
	if(taon->hasUserData("iPionCharge")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("iPionCharge")->size(); i++){PionsCharge.push_back(taon->userData<std::vector<double> >("iPionCharge")->at(i));  }
	}
	if(taon->hasUserData("iRefitPionP4")){
	  for(unsigned int i =0; i < taon->userData<std::vector<std::vector<double> > >("iRefitPionP4")->size(); i++){RefitPionsP4.push_back(taon->userData<std::vector<std::vector<double> > >("iRefitPionP4")->at(i));
	  }
	}
	if(taon->hasUserData("iRefitPionCharge")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("iRefitPionCharge")->size(); i++){RefitPionsCharge.push_back(taon->userData<std::vector<double> >("iRefitPionCharge")->at(i)); }
	}
	if(taon->hasUserData("PFTau_a1_cov")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double>  >("PFTau_a1_cov")->size(); i++) {ia1_cov.push_back(taon->userData<std::vector<double>  >("PFTau_a1_cov")->at(i)); }
	}
	if(taon->hasUserData("PFTau_a1_lvp")){
	  for(unsigned int i =0; i < taon->userData<std::vector<double> >("PFTau_a1_lvp")->size(); i++){ia1_par.push_back(taon->userData<std::vector<double> >("PFTau_a1_lvp")->at(i));  }
	}
      }
      ntaus++;
    }
    //  std::cout<<"   "<< iFLSign <<std::endl;

    TauFLSignificance.push_back(iFLSign);
    PFTauLeadTrackLV.push_back(iPFTauTrackLV);
    PFTauTrack_deltaR.push_back(iPFTauTrack_deltaR);
    PFTauGEOMFlightLenght.push_back(GEOMFlightLenght);
    PFTauGEOMFlightLenghtSignificance.push_back(GEOMFlightLenghtSignificance);
    _PFTauSVPos.push_back(SVPos);
    _PFTauSVCov.push_back(SVCov);
    //  std::cout<<" size  "<<TauFLSignificance.size() <<"   " <<PFTauGEOMFlightLenght.size()  <<std::endl;


    _PFTauSVChi2NDofMatchingQuality.push_back(SVChi2NDof);
    _PFTauPionsP4.push_back(PionsP4);
    _PFTauRefitPionsP4.push_back(RefitPionsP4);
    _PFTauPionsCharge.push_back(PionsCharge);
    _PFTauRefitPionsCharge.push_back(RefitPionsCharge);
    PFTau_a1_M.push_back(ia1_M);
    PFTau_a1_B.push_back(ia1_B);
    PFTau_a1_pdgid.push_back(ia1_pdgid);
    PFTau_a1_charge.push_back(ia1_Charge);
    PFTau_a1_cov.push_back(ia1_cov);
    PFTau_a1_lvp.push_back(ia1_par);

    if(ia1_par.size() != 0){
      TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
      TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
      int ltau=0;
      if(ia1_par.size()==LorentzVectorParticle::NLorentzandVertexPar){
        for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
          a1_par(k,0)=ia1_par.at(k);
          for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
            a1_cov(k,j)=ia1_cov.at(ltau);
            a1_cov(j,k)=ia1_cov.at(ltau);
            ltau++;
          }
        }
      }
      A1LVP.push_back(LorentzVectorParticle(a1_par,a1_cov,ia1_pdgid,ia1_Charge,ia1_B));
    }
    else A1LVP.push_back(LorentzVectorParticle());

    //_discriminator.push_back(discr);
    _daughters_typeOfMuon.push_back(typeOfMuon);
    _daughters_muonID.push_back(muIDflag);
    _daughters_tauID.push_back(tauIDflag);
    _daughters_footprintCorrection.push_back(footprintCorrection);
    _daughters_neutralIsoPtSumWeight.push_back(neutralIsoPtSumWeight);
    _daughters_photonPtSumOutsideSignalCone.push_back(photonPtSumOutsideSignalCone);

    _daughters_charge.push_back(cand->charge());
    //    _daughters_iseleBDT.push_back(isgood);
    _daughters_iseleWPLoose.push_back(iseleLoose);
    _daughters_iseleWP80.push_back(isele80);
    _daughters_iseleWP90.push_back(isele90);

    _daughters_iseleNoIsoWPLoose.push_back(iselenoisoLoose);
    _daughters_iseleNoIsoWP80.push_back(iselenoiso80);
    _daughters_iseleNoIsoWP90.push_back(iselenoiso90);
    _daughters_iseleCutBased.push_back(iselecutbased);

    _daughters_eleMVAnt.push_back(elemva);
    //    _daughters_eleMVAntNoIso.push_back(elemvanoiso);
    _daughters_eleMVA_HZZ.push_back(elemva_HZZ);
    _daughters_passConversionVeto.push_back(isconversionveto);
    _daughters_eleMissingHits.push_back(elemissinghits);
    //_daughters_eleMissingLostHits.push_back(elemissinglosthits);
    _daughters_iseleChargeConsistent.push_back(iselechargeconsistent);

    //_daughters_iseleCUT.push_back(userdatahelpers::getUserInt(cand,"isCUT"));
    _decayType.push_back(decay);
    _genmatch.push_back(genmatch);
    _MVADM2017v1.push_back(MVADM2017v1);
    _daughters_IetaIeta.push_back(ieta);
    _daughters_full5x5_IetaIeta.push_back(full5x5_ieta);
    _daughters_hOverE.push_back(hOverE);
    _daughters_deltaEtaSuperClusterTrackAtVtx.push_back(etasuperatvtx);
    _daughters_deltaPhiSuperClusterTrackAtVtx.push_back(phisuperatvtx);
    _daughters_IoEmIoP.push_back(IoEmIoP);
    _daughters_IoEmIoP_ttH.push_back(IoEmIoP_ttH);
    //_daughters_SCeta.push_back(SCeta);
    _daughters_depositR03_tracker.push_back(depositTracker);
    _daughters_depositR03_ecal.push_back(depositEcal);
    _daughters_depositR03_hcal.push_back(depositHcal);
    _daughters_decayModeFindingOldDMs.push_back(decayModeFindingOldDMs);
    _daughters_decayModeFindingNewDMs.push_back(decayModeFindingNewDMs);
    _daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(byCombinedIsolationDeltaBetaCorrRaw3Hits);
    _daughters_chargedIsoPtSum.push_back(chargedIsoPtSum);
    _daughters_neutralIsoPtSum.push_back(neutralIsoPtSum);
    _daughters_puCorrPtSum.push_back(puCorrPtSum);
    _daughters_byIsolationMVArun2v1DBoldDMwLTraw.push_back(byIsolationMVArun2v1DBoldDMwLTraw);
    _daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017.push_back(byIsolationMVArun2017v2DBoldDMwLTraw2017); //FRA
    _daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017.push_back(byIsolationMVArun2017v1DBoldDMwLTraw2017); //FRA
    _daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017.push_back(byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017); //FRA
    _daughters_byIsolationMVArun2v1DBoldDMwLTraw.push_back(byIsolationMVArun2v1DBoldDMwLTraw);
    _daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(byVVLooseIsolationMVArun2017v2DBoldDMwLT2017); //FRA
    _daughters_byDeepTau2017v2p1VSjetraw.push_back(byDeepTau2017v2p1VSjetraw);
    _daughters_byDeepTau2017v2p1VSeraw.push_back(byDeepTau2017v2p1VSeraw);
    _daughters_byDeepTau2017v2p1VSmuraw.push_back(byDeepTau2017v2p1VSmuraw);
    _daughters_numChargedParticlesSignalCone.push_back(numChargedParticlesSignalCone);
    _daughters_numNeutralHadronsSignalCone.push_back(numNeutralHadronsSignalCone);
    _daughters_numPhotonsSignalCone.push_back(numPhotonsSignalCone);
    _daughters_numParticlesSignalCone.push_back(numParticlesSignalCone);
    _daughters_numChargedParticlesIsoCone.push_back(numChargedParticlesIsoCone);
    _daughters_numNeutralHadronsIsoCone.push_back(numNeutralHadronsIsoCone);
    _daughters_numPhotonsIsoCone.push_back(numPhotonsIsoCone);
    _daughters_numParticlesIsoCone.push_back(numParticlesIsoCone);
    _daughters_leadChargedParticlePt.push_back(leadChargedParticlePt);
    _daughters_trackRefPt.push_back(trackRefPt);
    Muon_M.push_back(iMuon_M);
    Muon_B.push_back(iMuon_B);
    Muon_pdgid.push_back(iMuon_pdgid);
    Muon_trackCharge.push_back(iMuon_trackCharge);
    Muon_cov.push_back(iMuon_cov);
    Muon_par.push_back(iMuon_par);

    if(iMuon_par.size() != 0) {
      TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1);
      TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar);
      int lmu=0;
      for(int k=0; k<TrackParticle::NHelixPar; k++){
        mu_par(k,0)=iMuon_par.at(k);
        for(int j=k; j<TrackParticle::NHelixPar; j++){
          mu_cov(k,j)=iMuon_cov.at(lmu);
          lmu++;
        }
      }
      MuonTrack.push_back(TrackParticle(mu_par,mu_cov,iMuon_pdgid,iMuon_M,iMuon_trackCharge,iMuon_B));
    }
    else MuonTrack.push_back(TrackParticle());

    PFTau_Track_M.push_back(iPFTau_Track_M);
    PFTau_Track_B.push_back(iPFTau_Track_B);
    PFTau_Track_pdgid.push_back(iPFTau_Track_pdgid);
    PFTau_Track_charge.push_back(iPFTau_Track_trackCharge);

    PFTau_Track_cov.push_back(iPFTau_Track_cov);
    PFTau_Track_par.push_back(iPFTau_Track_par);


    _dxy_innerTrack.push_back(dxy_innerTrack);
    _dz_innerTrack.push_back(dz_innerTrack);
    _daughters_rel_error_trackpt.push_back(error_trackpt);
    _SIP.push_back(sip);

    //    _daughters_jetNDauChargedMVASel.push_back(jetNDauChargedMVASel);
    _daughters_miniRelIsoCharged.push_back(miniRelIsoCharged);
    _daughters_miniRelIsoNeutral.push_back(miniRelIsoNeutral);
    _daughters_jetPtRel.push_back(jetPtRel);
    _daughters_jetPtRatio.push_back(jetPtRatio);
    _daughters_jetBTagCSV.push_back(jetBTagCSV);
    _daughters_jetBTagDeepCSV.push_back(jetBTagDeepCSV);
    _daughters_jetBTagDeepFlavor.push_back(jetBTagDeepFlavor);
    //    _daughters_lepMVA_mvaId.push_back(lepMVA_mvaId);

    _daughters_pca_x.push_back(pcaPV.X());
    _daughters_pca_y.push_back(pcaPV.Y());
    _daughters_pca_z.push_back(pcaPV.Z());

    _daughters_pcaRefitPV_x.push_back(pcaRefitPV.X());
    _daughters_pcaRefitPV_y.push_back(pcaRefitPV.Y());
    _daughters_pcaRefitPV_z.push_back(pcaRefitPV.Z());

    _daughters_pcaRefitPVBS_x.push_back(pcaRefitPVBS_vecx);
    _daughters_pcaRefitPVBS_y.push_back(pcaRefitPVBS_vecy);
    _daughters_pcaRefitPVBS_z.push_back(pcaRefitPVBS_vecz);

    _daughters_pcaGenPV_x.push_back(pcaGenPV.X());
    _daughters_pcaGenPV_y.push_back(pcaGenPV.Y());
    _daughters_pcaGenPV_z.push_back(pcaGenPV.Z());

    _daughters_charged_px.push_back(chargedP4.X());
    _daughters_charged_py.push_back(chargedP4.Y());
    _daughters_charged_pz.push_back(chargedP4.Z());
    _daughters_charged_e.push_back(chargedP4.T());

    _daughters_neutral_px.push_back(neutralP4.X());
    _daughters_neutral_py.push_back(neutralP4.Y());
    _daughters_neutral_pz.push_back(neutralP4.Z());
    _daughters_neutral_e.push_back(neutralP4.T());

    //TRIGGER MATCHING
    Long64_t LFtriggerbit=0,L3triggerbit=0,filterFired=0;
    Long64_t trgMatched = 0;
    Long64_t triggertypeIsGood = 0;
    float hltpt=0;
    //cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    // list of indexes of all TO standalone that are matched to this specific daughter and pass HLT filter(s)
    // use as: toStandaloneMatched.at(idxHLT).at(idx tostadalone)
    vector<vector<int>> toStandaloneMatched (myTriggerHelper->GetNTriggers(), vector<int>(0));

    // for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    for (size_t idxto = 0; idxto < triggerObjects->size(); ++idxto)
      {
	pat::TriggerObjectStandAlone obj = triggerObjects->at(idxto);
	obj.unpackFilterLabels(event,*triggerBits);
	//check if the trigger object matches cand
	bool triggerType=false;

	if(deltaR2(obj,*cand)<0.25){

	  //cout << "######### NEW OBJECT MATCHED to offline " << cand->pdgId() << " of pt = " << cand->pt() << " HLT obj pt " << obj.pt() << endl;

	  if (type==ParticleType::TAU && (obj.hasTriggerObjectType(trigger::TriggerTau)|| obj.hasTriggerObjectType(trigger::TriggerL1TauJet)))triggerType=true;
	  if (type==ParticleType::ELECTRON && (obj.hasTriggerObjectType(trigger::TriggerElectron) || obj.hasTriggerObjectType(trigger::TriggerPhoton)))triggerType=true;
	  if (type==ParticleType::MUON && (obj.hasTriggerObjectType(trigger::TriggerMuon)))triggerType=true;
	  //check fired paths
	  obj.unpackPathNames(names);
	  std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	  std::vector<std::string> pathNamesLast = obj.pathNames(true);

	  for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

	    int triggerbit = myTriggerHelper->FindTriggerNumber(pathNamesAll[h],true);
	    if (triggerbit < 0) continue ; // not a path I want to save
	    bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
	    bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );

	    triggerMapper trgmap = myTriggerHelper->GetTriggerMap(pathNamesAll[h]);
	    bool isfilterGood = true;
	    int IDsearch = 0;
	    if (type==ParticleType::ELECTRON) IDsearch = 11;
	    else if (type==ParticleType::MUON) IDsearch = 13;
	    else if(type==ParticleType::TAU) IDsearch = 15;
	    int legPosition = trgmap.GetLegFromID(IDsearch);

	    if (legPosition == 1)
	      {
		for(int ifilt=0;ifilt<trgmap.GetNfiltersleg1();ifilt++)
		  {
		    string label = trgmap.Getfilter(true,ifilt);
		    if (label.empty()) continue;
		    if(! obj.hasFilterLabel(label.c_str()))isfilterGood=false;
		  }
	      }
	    else if (legPosition == 2)
	      {
		for(int ifilt=0;ifilt<trgmap.GetNfiltersleg2();ifilt++)
		  {
		    string label = trgmap.Getfilter(false,ifilt);
		    if (label.empty()) continue;
		    if(! obj.hasFilterLabel(label.c_str()))isfilterGood=false;
		  }
	      }
	    else isfilterGood = false;

	    //_isFilterFiredLast;
	    if(isfilterGood)filterFired |= (unsigned long long)(1) <<triggerbit;
	    if(triggerType) triggertypeIsGood |= (unsigned long long)(1) << triggerbit;
	    if(isLF)LFtriggerbit |= (unsigned long long)(1) <<triggerbit;
	    if(isL3)L3triggerbit |= (unsigned long long)(1) <<triggerbit;
	  } // loop on all trigger paths

	    // -------------- now do matching "filter-wise" to do x-check
	    // trigger matching checking labels
	  const std::vector<std::string>& vLabels = obj.filterLabels();

	  for (int triggerbit = 0; triggerbit < myTriggerHelper->GetNTriggers(); ++triggerbit)
	    {
	      triggerMapper trgmap = myTriggerHelper->GetTriggerMap(triggerbit);
	      bool istrgMatched = true;
	      int IDsearch = 0;
	      if (type==ParticleType::ELECTRON)  IDsearch = 11;
	      else if (type==ParticleType::MUON) IDsearch = 13;
	      else if(type==ParticleType::TAU)   IDsearch = 15;
	      int legPosition = trgmap.GetLegFromID(IDsearch);

	      if (legPosition == 1)
		{
		  for(int ifilt=0;ifilt<trgmap.GetNfiltersleg1();ifilt++)
		    {
		      string label = trgmap.Getfilter(true,ifilt);
		      //if(type==ParticleType::TAU)cout << " @@ leg 1 looking for " << label << endl;
		      if (label.empty()) continue;
		      if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) istrgMatched=false;

		    }
		}

	      else if (legPosition == 2)
		{
		  for(int ifilt=0;ifilt<trgmap.GetNfiltersleg2();ifilt++)
		    {
		      string label = trgmap.Getfilter(false,ifilt);
		      //if(type==ParticleType::TAU)cout << " @@ leg 2 looking for " << label << endl;
		      if (label.empty()) continue;
		      if (find(vLabels.begin(), vLabels.end(), label) == vLabels.end()) istrgMatched=false;
		    }
		}
	      else istrgMatched = false;
	      //if(type==ParticleType::TAU)cout<<istrgMatched<<endl;
	      // FIXME: should I check type? --> no, multiple filters should be enough
	      if(istrgMatched)
		{
		  //cout<<"triggerbit: "<<triggerbit<<endl;
		  //trgMatched |= (long(1) <<triggerbit);
		  trgMatched |= ((unsigned long long)(1) <<triggerbit);
		  //cout<<"triggerbit: "<<triggerbit<<endl;
		  toStandaloneMatched.at(triggerbit).push_back(idxto);
		}
	      // cout << "istrgMatched ? " << istrgMatched << endl;

	    } // loop on triggerbit from 0 to GetNTriggers()

	} // if dR < 0.25
      } // loop on all trigger candidates

    _daughters_isGoodTriggerType.push_back(triggertypeIsGood);
    _daughters_FilterFired.push_back(filterFired);
    _daughters_L3FilterFired.push_back(LFtriggerbit);
    _daughters_L3FilterFiredLast.push_back(L3triggerbit);
    _daughters_trgMatched.push_back(trgMatched);
    _daughters_HLTpt.push_back(hltpt);

    vector<int> vTrgMatchedIdx;
    for (int idxHLT = 0; idxHLT < myTriggerHelper->GetNTriggers(); ++idxHLT)
      {
	// if I have more than 1 match, I will be sure that different hlt objects are matched in a pair,
	// so I don't care and I put a value of -1
	// if I have no match, I don't care so I put -1 as well
	if (toStandaloneMatched.at(idxHLT).size() != 1)
	  vTrgMatchedIdx.push_back(-1);
	// if I have exactly 1 match, I store the Trigger Object index.
	// I will compare it later for the two legs of the pair and see if I matched separate objects
	else
	  vTrgMatchedIdx.push_back(toStandaloneMatched.at(idxHLT).at(0));
      }
    vTrgMatchedToDau_idx.push_back(vTrgMatchedIdx);


    // L1 candidate matching -- to correct for the missing seed
    bool isL1IsoTauMatched = false;

    std::vector<Float_t> L1IsoTau_et;
    for (int ibx = l1taus->getFirstBX(); ibx <= l1taus->getLastBX(); ++ibx)
      {
	for (BXVector<l1t::Tau>::const_iterator it=l1taus->begin(ibx); it!=l1taus->end(ibx); it++)
	  {
	    if (it->et() > 0 && ibx ==0){
	      if (it->hwIso() > 0.5){
		TLorentzVector tlv_L1Tau;
		TLorentzVector tlv_Tau;
		tlv_L1Tau.SetPtEtaPhiM(it->et(),
				       it->eta(),
				       it->phi(),
				       0.);

		tlv_Tau.SetPtEtaPhiM(cand->pt(),
				     cand->eta(),
				     cand->phi(),
				     0.);

		if ((tlv_L1Tau.DeltaR(tlv_Tau)*tlv_L1Tau.DeltaR(tlv_Tau)) < 0.25 && tlv_L1Tau.Pt()>32) { //GB
		  isL1IsoTauMatched = true;
		  L1IsoTau_et.push_back(it->et());
		}

	      }

	    }
	  }

      }

    if(isL1IsoTauMatched) {
      std::sort(L1IsoTau_et.begin(), L1IsoTau_et.end());
      Float_t L1IsoTau_etMax = *L1IsoTau_et.rbegin();
      _daughters_highestEt_L1IsoTauMatched.push_back(L1IsoTau_etMax) ;
    }else{
      _daughters_highestEt_L1IsoTauMatched.push_back(-1) ;
    }
  }
  //std::cout<<" size  "<<TauFLSignificance.size() <<"   " <<PFTauGEOMFlightLenght.size()  <<std::endl;

}

void HTauTauNtuplizer::FillGenInfo(const  edm::Event& event)
{
  edm::Handle<edm::View<pat::GenericParticle> > candHandle;

  //event.getByLabel ("genInfo", candHandle);
  event.getByToken (theGenericTag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();

  // Handle<edm::View<reco::GenParticle> > prunedHandle;
  // event.getByToken("prunedGenParticles", prunedHandle);

  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen)
    {
      // fill gen particle branches
      _genpart_px.push_back(igen->px());
      _genpart_py.push_back(igen->py());
      _genpart_pz.push_back(igen->pz());
      _genpart_e.push_back(igen->energy());
      _genpart_pdg.push_back(igen->pdgId());
      _genpart_status.push_back(igen->status());

      int HMIndex = -1;
      int MSSMHMIndex = -1;
      int TopMIndex = -1;
      int TauMIndex = -1;
      int ZMIndex = -1;
      int WMIndex = -1;
      int bMIndex = -1;
      int HZDecayMode = -1;
      int TopDecayMode = -1;
      int WDecayMode = -1;
      int TauGenDecayMode = -1;
      int TauGenDetailedDecayMode = -1;
      TVector3 pca(99,99,99);

      //	if (igen->hasUserInt("DataMC_Type_idx"))   _DataMC_Type = igen->userInt("DataMC_Type_idx");

      if (igen->hasUserInt("HMothIndex"))   HMIndex = igen->userInt("HMothIndex");
      if (igen->hasUserInt("MSSMHMothIndex"))   MSSMHMIndex = igen->userInt("MSSMHMothIndex");
      if (igen->hasUserInt("TopMothIndex")) TopMIndex = igen->userInt("TopMothIndex");
      if (igen->hasUserInt("TauMothIndex")) TauMIndex = igen->userInt("TauMothIndex");
      if (igen->hasUserInt("ZMothIndex"))   ZMIndex = igen->userInt("ZMothIndex");
      if (igen->hasUserInt("WMothIndex")) WMIndex = igen->userInt("WMothIndex");
      if (igen->hasUserInt("bMothIndex")) bMIndex = igen->userInt("bMothIndex");
      if (igen->hasUserInt("HZDecayMode"))   HZDecayMode = igen->userInt("HZDecayMode");
      if (igen->hasUserInt("TopDecayMode"))   TopDecayMode = igen->userInt("TopDecayMode");
      if (igen->hasUserInt("WDecayMode"))   WDecayMode = igen->userInt("WDecayMode");
      if (igen->hasUserInt("tauGenDecayMode"))   TauGenDecayMode = igen->userInt("tauGenDecayMode");
      if (igen->hasUserInt("tauGenDetailedDecayMode"))   TauGenDetailedDecayMode = igen->userInt("tauGenDetailedDecayMode");
      if (igen->hasUserFloat("pca_x")) pca = TVector3(igen->userFloat("pca_x"),igen->userFloat("pca_y"),igen->userFloat("pca_z"));

      _genpart_HMothInd.push_back(HMIndex);
      _genpart_MSSMHMothInd.push_back(MSSMHMIndex);
      _genpart_TopMothInd.push_back(TopMIndex);
      _genpart_TauMothInd.push_back(TauMIndex);
      _genpart_ZMothInd.push_back(ZMIndex);
      _genpart_WMothInd.push_back(WMIndex);
      _genpart_bMothInd.push_back(bMIndex);
      _genpart_HZDecayMode.push_back(HZDecayMode);
      _genpart_TopDecayMode.push_back(TopDecayMode);
      _genpart_WDecayMode.push_back(WDecayMode);
      _genpart_TauGenDecayMode.push_back(TauGenDecayMode);
      _genpart_TauGenDetailedDecayMode.push_back(TauGenDetailedDecayMode);
      _genpart_pca_x.push_back(pca.X());
      _genpart_pca_y.push_back(pca.Y());
      _genpart_pca_z.push_back(pca.Z());

      //const pat::GenericParticle* genClone = &(*igen);
      //int flags = CreateFlagsWord (genClone);
      int flags = igen -> userInt ("generalGenFlags");
      _genpart_flags.push_back(flags);
      if(igen->hasUserInt("HZDecayMode") || igen->hasUserInt("WDecayMode") || igen->hasUserInt("TopDecayMode")){
	_pvGen_x = igen->vx();
	_pvGen_y = igen->vy();
	_pvGen_z = igen->vz();
      }
    }
}


void HTauTauNtuplizer::FillGenJetInfo(const edm::Event& event)
{
  edm::Handle<edm::View<reco::GenJet>> genJetHandle;
  //event.getByLabel ("slimmedGenJets", genJetHandle);
  event.getByToken (theGenJetTag, genJetHandle);
  unsigned int genJetSize = genJetHandle->size();
  // to retrieve gen jet flavour from matched gen jet
  edm::Handle<edm::View<pat::Jet>> patjetHandle;
  edm::Handle<edm::View<pat::Jet>> patSmearedjetHandle;
  edm::Handle<edm::View<pat::Jet>> patSmearedjetDownHandle;
  edm::Handle<edm::View<pat::Jet>> patSmearedjetUpHandle;
  //event.getByLabel("jets", patjetHandle);
  event.getByToken(theJetTag, patjetHandle);
  event.getByToken(theSmearedJetTag, patSmearedjetHandle);
  event.getByToken(theSmearedJetDownTag, patSmearedjetDownHandle);
  event.getByToken(theSmearedJetUpTag, patSmearedjetUpHandle);
  unsigned int jetSize = patjetHandle->size();

  for (unsigned int igj = 0; igj < genJetSize; igj++)
    {
      const reco::GenJet& genJet = (*genJetHandle)[igj];
      _genjet_px.push_back ( genJet.px() );
      _genjet_py.push_back ( genJet.py() );
      _genjet_pz.push_back ( genJet.pz() );
      _genjet_e .push_back ( genJet.energy() );

      // jet flavour
      int partFlav = -999;
      int hadrFlav = -999;

      for (unsigned int ijet = 0; ijet < jetSize; ijet++)
	{
	  const pat::Jet& patjet = (*patjetHandle)[ijet];
	  const reco::GenJet * thismatchedGenJet =  patjet.genJet();

	  if (thismatchedGenJet == &genJet)
	    {
	      // no error, checked :-)
	      //if (partFlav != -999 && patjet.partonFlavour() != partFlav) cout << igj << " MISMATCH! Part flav: " << partFlav << " " << patjet.partonFlavour() << endl;
	      //if (hadrFlav != -999 && patjet.hadronFlavour() != hadrFlav) cout << igj << " MISMATCH! Hadr flav: " << hadrFlav << " " << patjet.hadronFlavour() << endl;
	      partFlav = patjet.partonFlavour();
	      hadrFlav = patjet.hadronFlavour();
	      break;
	    }
	}

      _genjet_partonFlavour.push_back(partFlav);
      _genjet_hadronFlavour.push_back(hadrFlav);


    }

  return;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void HTauTauNtuplizer::fillMCTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (!iEvent.isRealData()  || IsEmbed) {
    TauDecay_CMSSW myTauDecay;
    DataMCType DMT;
    DataMC_Type_idx = DMT.GetType();
    Handle <edm::View<reco::GenParticle> > genHandle;
    iEvent.getByToken (ThePrunedGenTag_, genHandle);
    myTauDecay.CheckForSignal(DataMC_Type_idx, genHandle);
    _DataMC_Type=DataMC_Type_idx;
    //    std::cout<<"HTauTauNtuplizer:  _DataMC_Type "<< _DataMC_Type <<std::endl;
    if (do_MCComplete_ || IsEmbed) {
      for (unsigned int iGenPartilce = 0; iGenPartilce < genHandle->size(); iGenPartilce++){
	const GenParticle& genP = (*genHandle)[iGenPartilce];
	if ( !isGoodGenParticle(genP) ) continue;
	MC_pdgid.push_back(genP.pdgId());
	MC_charge.push_back(genP.charge());
	std::vector<float> iMC_p4;
	iMC_p4.push_back(genP.p4().E());
	iMC_p4.push_back(genP.p4().Px());
	iMC_p4.push_back(genP.p4().Py());
	iMC_p4.push_back(genP.p4().Pz());

	MC_p4.push_back(iMC_p4);
	MC_midx.push_back(-1);
	MC_status.push_back(genP.status());
	MC_childpdgid.push_back(std::vector<int>());
	MC_childidx.push_back(std::vector<int>());
      }
      unsigned int i = 0;
      for (unsigned int iGenPartilce = 0; iGenPartilce < genHandle->size(); iGenPartilce++){

	const GenParticle& genP = (*genHandle)[iGenPartilce];
	if ( !isGoodGenParticle(genP) ) continue;
	for(unsigned int d = 0; d <(*genHandle)[iGenPartilce].numberOfDaughters(); d++){

	  const reco::Candidate *dau=(*genHandle)[iGenPartilce].daughter(d);
	  unsigned int j = 0;
	  for (unsigned int jGenPartilce = 0; jGenPartilce < genHandle->size(); jGenPartilce++){
	    const GenParticle& jgenP = (*genHandle)[jGenPartilce];

	    if ( !isGoodGenParticle(jgenP) ) continue;
	    if (dau->status() ==jgenP.status() && dau->p4() == jgenP.p4() && dau->pdgId() == jgenP.pdgId() && dau->numberOfMothers() == jgenP.numberOfMothers()
		&& dau->numberOfDaughters() == jgenP.numberOfDaughters()) {
	      MC_midx.at(j) = i;
	      MC_childidx.at(i).push_back(j);
	      MC_childpdgid.at(i).push_back(dau->pdgId());
	    }
	    j++;
	  }
	}
	i++;
      }
    }
    if (do_MCSummary_ || IsEmbed) {
      DataMCType DMT;
      for (unsigned int iGenPartilce = 0; iGenPartilce < genHandle->size(); iGenPartilce++){

	const GenParticle& genP = (*genHandle)[iGenPartilce];
	if(DMT.isSignalParticle((*genHandle)[iGenPartilce].pdgId()) && genP.numberOfDaughters() > 1){

	  MCSignalParticle_childpdgid.push_back(std::vector<int>());
	  MCSignalParticle_pdgid.push_back(genP.pdgId());
	  MCSignalParticle_charge.push_back(genP.charge());
	  MCSignalParticle_Tauidx.push_back(std::vector<unsigned int>());
	  std::vector<double> iSig_Poca;
	  iSig_Poca.push_back(genP.vx());
	  iSig_Poca.push_back(genP.vy());
	  iSig_Poca.push_back(genP.vz());
	  MCSignalParticle_Poca.push_back(iSig_Poca);

	  std::vector<double> iSig_p4;
	  iSig_p4.push_back(genP.p4().E());
	  iSig_p4.push_back(genP.p4().Px());
	  iSig_p4.push_back(genP.p4().Py());
	  iSig_p4.push_back(genP.p4().Pz());
	  MCSignalParticle_p4.push_back(iSig_p4);
	  // look for daughter tau

	  for(unsigned int i = 0; i <(*genHandle)[iGenPartilce].numberOfDaughters(); i++){

	    const reco::Candidate *dau=(*genHandle)[iGenPartilce].daughter(i);
	    //	  MCSignalParticle_childpdgid.at(MCSignalParticle_childpdgid.size() - 1).push_back(dau->pdgId());
	    if (abs(dau->pdgId()) == PDGInfo::tau_minus) {

	      unsigned int tauidx = MCTauandProd_p4.size();
	      MCSignalParticle_Tauidx.at(MCSignalParticle_Tauidx.size() - 1).push_back(tauidx);
	      // Analysis the tau decay
	      unsigned int JAK_ID, TauBitMask;
	      myTauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau), JAK_ID, TauBitMask);
	      std::vector<const reco::GenParticle*> TauDecayProducts = myTauDecay.Get_TauDecayProducts();
	      MCTauandProd_midx.push_back(myTauDecay.Get_MotherIdx());
	      MCTau_JAK.push_back(JAK_ID);
	      //		cout<<"check JAK ID " <<  JAK_ID <<endl;
	      MCTau_DecayBitMask.push_back(TauBitMask);
	      MCTauandProd_pdgid.push_back(std::vector<int>());
	      MCTauandProd_charge.push_back(std::vector<int>());
	      MCTauandProd_p4.push_back(std::vector<std::vector<double> >());
	      MCTauandProd_Vertex.push_back(std::vector<std::vector<double> >());
	      // if(_DataMC_Type==10130533 || _DataMC_Type==10230533)std::cout<<_DataMC_Type << std::endl;
	      // double apx(0),apy(0), apz(0), ae(0);
	      for (unsigned int i = 0; i < TauDecayProducts.size(); i++) {
		//cout<<TauDecayProducts.at(i)->pdgId()<<endl;
		MCTauandProd_pdgid.at(tauidx).push_back(TauDecayProducts.at(i)->pdgId());
		MCTauandProd_charge.at(tauidx).push_back(TauDecayProducts.at(i)->charge());

		std::vector<double> iTauandProd_p4;
		std::vector<double> iTauandProd_vertex;

		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().E());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Px());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Py());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Pz());

		iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vx());
		iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vy());
		iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vz());

		MCTauandProd_p4.at(tauidx).push_back(iTauandProd_p4);
		MCTauandProd_Vertex.at(tauidx).push_back(iTauandProd_vertex);
	      }
	      // if((_DataMC_Type==10130533 || _DataMC_Type==10230533 )&& apx!=0){

	      // TLorentzVector a1(apx,apy,apz,ae);
	      // std::cout<<" gen a1 Mass and Px "<< a1.M() <<"   "<< a1.Px() << std::endl;
	      // }
	    }
	  }
	}
      }
    }
  }
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool HTauTauNtuplizer::EleVeto(const reco::Candidate* cand)
{
  //cout<<endl<<cand->isElectron()<<" ";
  if(cand->isElectron()){
    //cout<<endl<<"Electron:";
    //cout<<" "<<TauPt<<" "<<std::abs(cand->eta())<<" ";
    if(cand->pt()>10 && std::abs(cand->eta())<2.5){
      //  cout<<" "<<std::abs(userdatahelpers::getUserFloat(cand,"dxy"))<<" "<<std::abs(userdatahelpers::getUserFloat(cand,"dz"))<<" ";
      if(std::abs(userdatahelpers::getUserFloat(cand,"dxy")) < 0.045 && std::abs(userdatahelpers::getUserFloat(cand,"dz")) < 0.2){
	//	cout<<" "<<userdatahelpers::getUserFloat(cand,"isEleID90")<<" ";
	if(userdatahelpers::getUserFloat(cand,"isEleNoIsoID90") == 1){
	  //	  cout<<" "<<userdatahelpers::getUserInt(cand,"isConversionVeto")<<" ";
	  if(userdatahelpers::getUserInt(cand,"isConversionVeto") == 1){
	    //	    cout<<" "<<userdatahelpers::getUserInt(cand,"missingHit")<<" ";
	    if(userdatahelpers::getUserInt(cand,"missingHit")<=1){
	      //	      cout<<" "<<userdatahelpers::getUserFloat(cand,"combRelIsoPF")<<" ";
	      if(userdatahelpers::getUserFloat(cand,"combRelIsoPF")<(0.3*cand->pt())){
		return true;
	      }
	    }
	  }
	}
      }
    }
  }
  return false;
}
bool HTauTauNtuplizer::MuVeto(const reco::Candidate* cand)
{
  //  cout<<endl<<cand->isMuon()<<" ";
  if(cand->isMuon()){
    //cout<<endl<<"Muon:";
    //  cout<<" "<<TauPt<<" "<<std::abs(cand->eta())<<" ";
    if(cand->pt()>10 && std::abs(cand->eta())<2.4){
      //      cout<<" "<<std::abs(userdatahelpers::getUserFloat(cand,"dxy"))<<" "<<std::abs(userdatahelpers::getUserFloat(cand,"dz"))<<" ";
      if(std::abs(userdatahelpers::getUserFloat(cand,"dxy")) < 0.045 && std::abs(userdatahelpers::getUserFloat(cand,"dz")) < 0.2){
	//	cout<<" "<<CHECK_BIT(userdatahelpers::getUserInt(cand,"muonID"),2)<<" ";
	if(CHECK_BIT(userdatahelpers::getUserInt(cand,"muonID"),2)){
	  //	  cout<<" "<<userdatahelpers::getUserFloat(cand,"combRelIsoPF")<<" ";
	  if(userdatahelpers::getUserFloat(cand,"combRelIsoPF")<(0.3*cand->pt())){
	    return true;
	  }
	}
      }
    }
  }
  return false;
}

bool HTauTauNtuplizer::DiEle(const reco::Candidate* cand1,const reco::Candidate* cand2)
{
  bool kin=false, vertex=false, isele=false, iso=false;

  if(cand1->isElectron() && cand2->isElectron() &&  ((cand1->charge()/abs(cand1->charge())) != (cand2->charge()/abs(cand2->charge())))){
    {
      if(deltaR(cand1->p4(),cand2->p4())>0.15)
	{
	  kin = (cand1->pt()>15 && cand2->pt()>15 && fabs(cand1->eta())<2.5 && fabs(cand2->eta())<2.5);
	  vertex = (fabs(userdatahelpers::getUserFloat(cand1,"dxy"))<0.045 && fabs(userdatahelpers::getUserFloat(cand2,"dxy"))< 0.045 && fabs(userdatahelpers::getUserFloat(cand1,"dz")) < 0.2 && fabs(userdatahelpers::getUserFloat(cand2,"dz")) < 0.2);
	  isele=(userdatahelpers::getUserFloat(cand1,"isEleNoIsoID90") && userdatahelpers::getUserFloat(cand2,"isEleNoIsoID90"));
	  iso= (userdatahelpers::getUserFloat(cand1,"combRelIsoPF") < (0.3*cand1->pt()) && userdatahelpers::getUserFloat(cand2,"combRelIsoPF")<(0.3*cand2->pt()));
	  if((kin && vertex && isele && iso)==1)return true;
	}
    }
  }
  return false;
}

bool HTauTauNtuplizer::DiMuon(const reco::Candidate* cand1,const reco::Candidate* cand2)
{
  bool kin=false, vertex=false, isele=false, iso=false;

  if(cand1->isMuon() && cand2->isMuon() &&  ((cand1->charge()/abs(cand1->charge())) != (cand2->charge()/abs(cand2->charge())))){
    {
      if(deltaR(cand1->p4(),cand2->p4())>0.15)
	{
	  kin = (cand1->pt()>15 && cand2->pt()>15 && fabs(cand1->eta())<2.4 && fabs(cand2->eta())<2.4);
	  vertex = (fabs(userdatahelpers::getUserFloat(cand1,"dxy"))<0.045 && fabs(userdatahelpers::getUserFloat(cand2,"dxy"))< 0.045 && fabs(userdatahelpers::getUserFloat(cand1,"dz")) < 0.2 && fabs(userdatahelpers::getUserFloat(cand2,"dz")) < 0.2);
	  isele=((userdatahelpers::getUserFloat(cand1,"isPFMuon") && userdatahelpers::getUserFloat(cand2,"isPFMuon")) && (userdatahelpers::getUserFloat(cand1,"isGlobalMuon") && userdatahelpers::getUserFloat(cand2,"isGlobalMuon")) && (userdatahelpers::getUserFloat(cand1,"isTrackerMuon") && userdatahelpers::getUserFloat(cand2,"isTrackerMuon")));
	  iso= (userdatahelpers::getUserFloat(cand1,"combRelIsoPF") < (0.3*cand1->pt()) && userdatahelpers::getUserFloat(cand2,"combRelIsoPF")<(0.3*cand2->pt()));
	  if((kin && vertex && isele && iso)==1)return true;
	}
    }
  }
  return false;
}

TLorentzVector HTauTauNtuplizer::getVisMomentumNoLep(const reco::Candidate* genL/*int tauindex,std::vector<std::vector<std::vector<double> > > MCTauandProd_p4,std::vector<std::vector<int> > MCTauandProd_pdgid*/)
{
  TLorentzVector lvLept(0,0,0,0);
  TLorentzVector lvZero(0,0,0,0);
  int NMCTauDecayProducts=0;
  //int status = userdatahelpers::getUserInt(genL,"status");
  double dR=0.5;
  double itau=9999;
  math::XYZTLorentzVector MCTau;
  for(unsigned int i_tau=0; i_tau < MCTauandProd_p4.size();i_tau++)
    {

      MCTau=math::XYZTLorentzVector(MCTauandProd_p4.at(i_tau).at(0).at(1),MCTauandProd_p4.at(i_tau).at(0).at(2),MCTauandProd_p4.at(i_tau).at(0).at(3),MCTauandProd_p4.at(i_tau).at(0).at(0));

      if (deltaR(genL->p4(),MCTau)<dR)
	{
	  dR=deltaR(genL->p4(),MCTau);
	  NMCTauDecayProducts=MCTauandProd_p4.at(i_tau).size();
	  itau=i_tau;
	}
    }
  if (dR==0.5) return lvZero;
  for(int i_dau=1; i_dau < NMCTauDecayProducts; i_dau++) {
    //if (status!=-1 && status != genP.status())continue;
    if( abs(MCTauandProd_pdgid.at(itau).at(i_dau)) == PDGInfo::nu_e ||
	abs(MCTauandProd_pdgid.at(itau).at(i_dau)) == PDGInfo::nu_mu ||
	abs(MCTauandProd_pdgid.at(itau).at(i_dau)) == 11 ||
	abs(MCTauandProd_pdgid.at(itau).at(i_dau)) == 13){

      lvLept += TLorentzVector(MCTauandProd_p4.at(itau).at(i_dau).at(1),MCTauandProd_p4.at(itau).at(i_dau).at(2),MCTauandProd_p4.at(itau).at(i_dau).at(3),MCTauandProd_p4.at(itau).at(i_dau).at(0));
    }
  }
  TLorentzVector lvvis = TLorentzVector(MCTauandProd_p4.at(itau).at(0).at(1),MCTauandProd_p4.at(itau).at(0).at(2),MCTauandProd_p4.at(itau).at(0).at(3),MCTauandProd_p4.at(itau).at(0).at(0));

  lvvis -= lvLept;
  return lvvis;
}

int HTauTauNtuplizer::GenMatch(const reco::Candidate* genL/*, int tauindex*/, const edm::Event& event/*,std::vector<std::vector<std::vector<double> > > MCTauandProd_p4,std::vector<std::vector<int> > MCTauandProd_pdgid*/)
{
  edm::Handle<edm::View<reco::GenParticle> > genHandle;
  event.getByToken (ThePrunedGenTag_, genHandle);

  TauDecay_CMSSW myTauDecay;
  int GenMatch = 6;
  //int status = userdatahelpers::getUserInt(genL,"status");
  //int id = userdatahelpers::getUserInt(genL,"id");
  float minDR = 0.2;

  //if (status == 99999 && id == 99999) return -1; // default values in *Filler.cc when no matched gen found
  //cout<<genHandle->size()<<endl;
  // get matched particle by looking at pdgId, px, py, pz, e --> not the most elegant, but should work
  for (unsigned int iGen = 0; iGen < genHandle->size(); iGen++)
    //for(std::vector<reco::GenParticle>::const_iterator genP = genHandle.begin();genP != genHandle.end(); ++genP )
    {
      const GenParticle/*pat::GenericParticle*/ &genP = (*genHandle)[iGen];
      const GenStatusFlags statusFlags = genP.statusFlags();
      // first check on pdgId only, then status to remove most of checks with a minimum amount of comparison

      //if (id == genP.pdgId())
      if (deltaR(genP.p4(),genL->p4())<minDR) {
	minDR = deltaR(genP.p4(),genL->p4());

	bool type1 = abs(genP.pdgId())==11 && statusFlags.isPrompt() && genP.pt()>8;
	bool type2 = abs(genP.pdgId())==13 && statusFlags.isPrompt() && genP.pt()>8;
	bool type3 = abs(genP.pdgId())==11 && statusFlags.isDirectPromptTauDecayProduct() && genP.pt()>8;
	bool type4 = abs(genP.pdgId())==13 && statusFlags.isDirectPromptTauDecayProduct() && genP.pt()>8;
	//cout<<abs(genP.pdgId())<< "  "<<statusFlags.isPrompt() <<"  "<<getVisMomentumNoLep(genL).Pt()<<endl;
	bool type5 = abs(genP.pdgId())==15 &&  statusFlags.isPrompt() &&  getVisMomentumNoLep(genL/*tauindex,MCTauandProd_p4, MCTauandProd_pdgid*/).Pt()>15;
	//if(abs(genP.pdgId())==15) cout<<(getVisMomentumNoLep(genL).Pt()>15)<<" "<<statusFlags.isPrompt()<<endl;
	if (type1) GenMatch = 1;
	else if (type2) GenMatch = 2;
	else if (type3) GenMatch = 3;
	else if (type4) GenMatch = 4;
	else if (type5) GenMatch = 5;
      }
    }
  //}
  return GenMatch;
}
// return index of gen matched to reco lepton, and -1 if not existing or not found
int HTauTauNtuplizer::GetMatchedGen (const reco::Candidate* genL, const edm::Event& event)
{
  //cout.precision(15); // just to check real precision

  edm::Handle<edm::View<pat::GenericParticle>>candHandle;
  //event.getByLabel ("genInfo", candHandle);
  event.getByToken (theGenericTag, candHandle);

  int index = -1;
  int status = userdatahelpers::getUserInt(genL,"status");
  int id = userdatahelpers::getUserInt(genL,"id");
  if (status == 99999 && id == 99999) return -1; // default values in *Filler.cc when no matched gen found

  // get matched particle by looking at pdgId, px, py, pz, e --> not the most elegant, but should work
  for (unsigned int iGen = 0; iGen < candHandle->size(); iGen++)
    {
      const pat::GenericParticle& genP = (*candHandle)[iGen];
      // first check on pdgId only, then status to remove most of checks with a minimum amount of comparison
      if (id == genP.pdgId())
	{
	  if (status == genP.status())
	    {
	      float px = userdatahelpers::getUserFloat(genL,"genPx");
	      float py = userdatahelpers::getUserFloat(genL,"genPy");
	      float pz = userdatahelpers::getUserFloat(genL,"genPz");
	      float e = userdatahelpers::getUserFloat(genL,"genE");
	      //cout << "     ==> I'm in with " << fixed<< px << " " << fixed<<py << " " << fixed<<pz << " " << fixed<<e << "  ||| " << fixed<<genP.px() << " " << fixed<<genP.py() << " " << fixed<<genP.pz() << " " << fixed<<genP.energy() << endl;
	      // cast to float helps in the EPSIL comparison, but for safety EPSIL = 1.e-5 is used (seems reasonably small for the value range we use in 0.001 -- 1000) with "meaningful" values O(10) or larger
	      if ( fabs((float)genP.px() - px) < EPSIL && fabs((float)genP.py() - py) < EPSIL && fabs((float)genP.pz() - pz) < EPSIL && fabs((float)genP.energy() - e) < EPSIL )
		{
		  index = iGen;
		  break; // found, no other comparisons needed
		}
	    }
	}
    }

  return index;
}
//Fill L1 objects
void HTauTauNtuplizer::FillL1Obj(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event){

  for (int ibx = taus->getFirstBX(); ibx <= taus->getLastBX(); ++ibx)
    {
      for (BXVector<l1t::Tau>::const_iterator it=taus->begin(ibx); it!=taus->end(ibx); it++)
	{
	  if (it->et() > 0 && ibx ==0){

	    _L1_tauEt .push_back(it->et());
	    _L1_tauEta.push_back(it->eta());
	    _L1_tauPhi.push_back(it->phi());
	    _L1_tauIso.push_back(it->hwIso());
	  }
	}
    }

  for (int ibx = jets->getFirstBX(); ibx <= jets->getLastBX(); ++ibx)
    {
      for (BXVector<l1t::Jet>::const_iterator it=jets->begin(ibx); it!=jets->end(ibx); it++)
	{
	  if (it->et() > 0&& ibx ==0){

	    _L1_jetEt .push_back(it->et());
	    _L1_jetEta.push_back(it->eta());
	    _L1_jetPhi.push_back(it->phi());

	  }
	}
    }
}

void HTauTauNtuplizer::endJob(){
  hCounter->SetBinContent(1,Nevt_Gen);
  hCounter->SetBinContent(2,Nevt_PassTrigger);
  hCounter->SetBinContent(3,Npairs);

  hCounter->GetXaxis()->SetBinLabel(1,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2,"Nevt_PassTrigger");
  hCounter->GetXaxis()->SetBinLabel(3,"Npairs");

  for(int i=0;i<myTriggerHelper->GetNTriggers();i++){
    hCounter->GetXaxis()->SetBinLabel(i+4,(myTriggerHelper->printTriggerName(i)).c_str());
    std::cout<<"i: "<<i<<" "<<"name: "<<myTriggerHelper->printTriggerName(i)<<endl;
  }
  for(int i=1;i<=ntauIds;i++){
    hTauIDs->GetXaxis()->SetBinLabel(i,tauIDStrings[i-1].Data());
  }

  //cout<<"trigincr: "<<trigincr<<endl;
  //cout<<"filtersincr: "<<filtersincr<<endl;


}


void HTauTauNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){

  // For 2018 data (run < 315974) the muon filters of the MuTau triggers (leg1) need to be customized
  if (!theisMC && theYear==2018 && iRun.run()<315974)
    {
      // HLT paths with "HPS"
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltHpsOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});

      // HLT paths without "HPS"
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"});
      if ( myTriggerHelper->HasTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v") )
	myTriggerHelper->ChangeTriggerMap("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",{"hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"});
    }

  Bool_t changedConfig = false;

  //if(!hltConfig_.init(iRun, iSetup, triggerResultsLabel.process(), changedConfig)){
  if(!hltConfig_.init(iRun, iSetup, processName.process(), changedConfig)){
    edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
    return;
  }

  if(changedConfig || foundPaths.size()==0){
    //cout<<"The present menu is "<<hltConfig.tableName()<<endl;
    indexOfPath.clear();
    foundPaths.clear();
    //for(size_t i=0; i<triggerPaths.size(); i++){
    // bool foundThisPath = false;
    for(size_t j=0; j<hltConfig_.triggerNames().size(); j++){
      string pathName = hltConfig_.triggerNames()[j];

      indexOfPath.push_back(j);
      foundPaths.push_back(pathName);

      cout << j << " - TTT: " << pathName << endl;
      //	  edm::LogInfo("AnalyzeRates")<<"Added path "<<pathName<<" to foundPaths";
    }
  }

}


void HTauTauNtuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}
void HTauTauNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup){
  if (theisMC)
    {
      // try {event.getByToken(theLHEPTag, lheEventProduct);} catch (...) {;}
      // if (lheEventProduct.isValid())
      edm::Handle<GenLumiInfoHeader> gen_header;
      iLumi.getByToken(genLumiHeaderTag, gen_header);
      if (gen_header.isValid())
	{
	  string model = gen_header->configDescription();
	  // cout<<model<<endl;  // prints, e.g. T1tttt_1500_100
	  _susyModel = model;
	}
    }
}
void HTauTauNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  iLumi.getByToken(theTotTag, nEventsTotalCounter);
  Nevt_Gen += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsPassTrigCounter;
  iLumi.getByToken(thePassTag, nEventsPassTrigCounter);
  Nevt_PassTrigger += nEventsPassTrigCounter->value;
}


bool HTauTauNtuplizer::CompareLegs(const reco::Candidate *i, const reco::Candidate *j){
  int iType=2,jType=2;

  if(i->isElectron())iType=1;
  else if(i->isMuon())iType=0;

  if(j->isElectron())jType=1;
  else if(j->isMuon())jType=0;

  if(iType>jType) return false;
  else if(iType==jType && i->pt()<j->pt()) return false;

  return true;
}

// implements operator "<" (return i < j)
bool HTauTauNtuplizer::ComparePairsbyIso(pat::CompositeCandidate i, pat::CompositeCandidate j){

  //Second criteria: ISO
  float isoi=999,isoj=999;
  int cand1j=-1,cand1i=-1;

  if(CompareLegs(i.daughter(0),i.daughter(1)))cand1i=0;
  else cand1i=1;
  if(CompareLegs(j.daughter(0),j.daughter(1)))cand1j=0;
  else cand1j=1;

  //step 1, leg 1 ISO
  //byCombinedIsolationDeltaBetaCorrRaw3Hits
  isoi=userdatahelpers::getUserFloat(i.daughter(cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(cand1j),"combRelIsoPF");
  if (!i.daughter(cand1i)->isMuon() && !i.daughter(cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(cand1i),"byDeepTau2017v2p1VSjetraw");
  if (!j.daughter(cand1j)->isMuon() && !j.daughter(cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(cand1j),"byDeepTau2017v2p1VSjetraw");
  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 2, leg 1 Pt
  if(i.daughter(cand1i)->pt()>j.daughter(cand1j)->pt()) return true;
  else if(i.daughter(cand1i)->pt()<j.daughter(cand1j)->pt()) return false;

  //step 3, leg 2 ISO
  isoi=userdatahelpers::getUserFloat(i.daughter(1-cand1i),"combRelIsoPF");
  isoj=userdatahelpers::getUserFloat(j.daughter(1-cand1j),"combRelIsoPF");
  if (!i.daughter(1-cand1i)->isMuon() && !i.daughter(1-cand1i)->isElectron()) isoi= -userdatahelpers::getUserFloat(i.daughter(1-cand1i),"byDeepTau2017v2p1VSjetraw");
  if (!j.daughter(1-cand1j)->isMuon() && !j.daughter(1-cand1j)->isElectron()) isoj= -userdatahelpers::getUserFloat(j.daughter(1-cand1j),"byDeepTau2017v2p1VSjetraw");

  if (isoi<isoj)return true;
  else if(isoi>isoj)return false;

  //step 4, leg 2 Pt
  if(i.daughter(1-cand1i)->pt()>j.daughter(1-cand1j)->pt()) return true;
  //else if(i.daughter(1-cand1i)->pt()<j.daughter(1-cand1j)->pt()) return false;

  return false;

}

bool HTauTauNtuplizer::ComparePairsbyPt(pat::Jet i, pat::Jet j){

  //First criteria: OS<SS
  //if ( j.charge()==0 && i.charge()!=0) return false;
  //else if ( i.charge()==0 && j.charge()!=0) return true;

  //Second criteria: legs pt
  if(i.pt()<j.pt()) return false;
  if(i.pt()>j.pt()) return true;

  return true;
}

bool HTauTauNtuplizer::CHECK_BIT(unsigned long long var, int pos){
  unsigned long long CHECK1=(unsigned long long)(var & ((unsigned long long)(1) << pos));
  unsigned long long CHECK2=(unsigned long long)((unsigned long long)(1) << pos);
  //cout<<"Check method: "<<"var: "<<var<<"  pos: "<<pos<<"   "<<CHECK1<<"  "<< CHECK2<<endl;
  return ( CHECK1==CHECK2 ); }

bool HTauTauNtuplizer::isGoodGenParticle(const reco::GenParticle &GenPar){
  //	if (GenPar.p4().Pt() > MCCompletePtCut_) return true;  // in case of you want to apply some cuts on it
  int id = abs(GenPar.pdgId());
  if (id == PDGInfo::Z0 || id == PDGInfo::W_plus) return true;
  if (id == PDGInfo::Higgs0) return true;
  if (id >= PDGInfo::Z_prime0 && id <= PDGInfo::Higgs_plus) return true; //BSM resonances
  if (id == PDGInfo::t) return true;
  return false;
}
float HTauTauNtuplizer::ComputeMT (math::XYZTLorentzVector visP4, float METx, float METy)
{
  math::XYZTLorentzVector METP4 (METx, METy, 0, 0); // I only care about transverse plane
  double dphi = deltaPhi(visP4.phi(), METP4.phi());
  return sqrt(2.*visP4.Pt()*METP4.Pt()*(1.-cos(dphi)));

  // this code is equivalent, but can give NaN and rounding errors when met and lepton are very collinear
  // float MET = TMath::Sqrt (METx*METx + METy*METy);
  // math::XYZTLorentzVector METP4 (METx, METy, 0, MET);

  // float scalSum = MET + visP4.pt();
  // math::XYZTLorentzVector vecSum (visP4);
  // vecSum += METP4;
  // float vecSumPt = vecSum.pt();

  // return TMath::Sqrt (scalSum*scalSum - vecSumPt*vecSumPt);
}
float HTauTauNtuplizer::ComputePUPPIMT (math::XYZTLorentzVector visP4,float PUPPImetPhi, float PUPPImetPt)
{
  double dphi = deltaPhi(visP4.phi(), PUPPImetPhi);
  return sqrt(2.*visP4.Pt()*PUPPImetPt*(1.-cos(dphi)));
}
float HTauTauNtuplizer::ComputeMTtot (math::XYZTLorentzVector visP41, math::XYZTLorentzVector visP42,float METx, float METy)
{
  math::XYZTLorentzVector METP4 (METx, METy, 0, 0); // I only care about transverse plane
  double dphi1 = deltaPhi(visP41.phi(), METP4.phi());
  double dphi2 = deltaPhi(visP42.phi(), METP4.phi());
  double dphi12 = deltaPhi(visP41.phi(), visP42.phi());
  return sqrt(2.*visP41.Pt()*METP4.Pt()*(1.-cos(dphi1))+2.*visP42.Pt()*METP4.Pt()*(1.-cos(dphi2))+2.*visP41.Pt()*visP42.Pt()*(1.-cos(dphi12)));
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool HTauTauNtuplizer::refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup){
  edm::Handle<edm::View<pat::PackedCandidate> >pfCandHandle;
  iEvent.getByToken(thePFCandTag,pfCandHandle);
  const edm::View<pat::PackedCandidate>* cands = pfCandHandle.product();

  edm::Handle<pat::TauCollection> tauHandle;
  iEvent.getByToken(theNewTauTag, tauHandle);

  Handle<vector<reco::Vertex> >  vertices;
  iEvent.getByToken(theVtxTag,vertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotTag, beamSpot);

  Handle <RefitVertexCollection> RefitPVBSHandle;
  iEvent.getByToken(RefitVtxBSTag,RefitPVBSHandle);

  Handle <RefitVertexCollection> RefitPVNoBSHandle;
  iEvent.getByToken(RefitVtxNoBSTag,RefitPVNoBSHandle);

  const RefitVertexCollection* RefitPVBS = RefitPVBSHandle.product();
  const RefitVertexCollection* RefitPVNoBS = RefitPVNoBSHandle.product();

  boost::hash<const reco::Candidate*> hasher;

  //Get tracks associated with pfPV
  reco::TrackCollection pvTracks;
  reco::TrackCollection allTracks;

  //Old
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
    if(!(*cands)[i].bestTrack()) continue;

    unsigned int key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();
    allTracks.push_back(*((*cands)[i].bestTrack()));
    if(key!=0 ||
       (quality!=pat::PackedCandidate::UsedInFitTight
	&& quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

    pvTracks.push_back(*((*cands)[i].bestTrack()));
  }

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);

  std::vector<reco::TransientTrack> transTracks;
  TransientVertex transVtx;
  for(auto iter: pvTracks) transTracks.push_back(transTrackBuilder->build(iter));   // <-  this contains tracks from tau decays and they should be removed

  bool fitOk = false;
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.001);
    try {
      transVtx = avf.vertex(transTracks, *beamSpot);
      fitOk = true;
    } catch (...) {
      fitOk = false;
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  }

  fitOk = fitOk && transVtx.isValid() && fabs(transVtx.position().x())<1 && fabs(transVtx.position().y())<1;

  if(fitOk) {
    _pvRefit_x = transVtx.position().x();
    _pvRefit_y = transVtx.position().y();
    _pvRefit_z = transVtx.position().z();

    reco::Vertex primarVertex = transVtx;

    TMatrixTSym<double> apvcov(3);
    //      svcov.ResizeTo(3,3);
    math::Error<3>::type pvCov;
    primarVertex.fill(pvCov);
    for (int i = 0; i <3; i++){
      for (int j = 0; j < 3; j++) {
	apvcov(i, j) = pvCov(i, j);
	apvcov(j, i) = pvCov(i, j);
      }
    }


    for (int i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
	_pvRefit_cov.push_back(apvcov(i, j));
      }
    }
  }
  else {
    _pvRefit_x = (*vertices)[0].x();
    _pvRefit_y = (*vertices)[0].y();
    _pvRefit_z = (*vertices)[0].z();
  }


  if(RefitPVBSHandle.isValid())
    {
      for(unsigned int irefit=0;irefit<RefitPVBS->size();++irefit)
	{
	  _RefitPVBS_x.push_back(RefitPVBS->at(irefit).x());
	  _RefitPVBS_y.push_back(RefitPVBS->at(irefit).y());
	  _RefitPVBS_z.push_back(RefitPVBS->at(irefit).z());

          vector<vector<double>> pvbs_cov(3,vector<double>(3));
	  math::Error<3>::type pvbsCov;
	  RefitPVBS->at(irefit).fill(pvbsCov);
	  for(int i=0; i<3; i++){
	    for(int j=i; j<3; j++){
	      pvbs_cov[i][j]=pvbsCov(i,j);
	      if(i!=j){
	        pvbs_cov[j][i]=pvbsCov(i,j);
	      }
	    }
          }
	  _RefitPVBS_Cov.push_back(pvbs_cov);

	  std::vector<size_t> hashesBS;
	  for(auto name: RefitPVBS->at(irefit).userCandNames())
	    {
	      edm::Ptr<reco::Candidate> aRecoCand = RefitPVBS->at(irefit).userCand( name );
	      size_t hash = hasher(aRecoCand.get());
	      hashesBS.push_back(hash);
	    }

	  _VertexHashBS1.push_back(hashesBS[0]);
	  _VertexHashBS2.push_back(hashesBS[1]);
	  //cout<<"Vertex1: "<<hashesBS[0]<<endl;
	  //cout<<"Vertex2: "<<hashesBS[1]<<endl;
	}
    }

  if(RefitPVNoBSHandle.isValid())
    {
      for(unsigned int jrefit=0;jrefit<RefitPVNoBS->size();++jrefit)
	{
	  _RefitPVNoBS_x.push_back(RefitPVNoBS->at(jrefit).x());
	  _RefitPVNoBS_y.push_back(RefitPVNoBS->at(jrefit).y());
	  _RefitPVNoBS_z.push_back(RefitPVNoBS->at(jrefit).z());

	  std::vector<size_t> hashesNoBS;
	  for(auto name: RefitPVNoBS->at(jrefit).userCandNames())
	    {
	      edm::Ptr<reco::Candidate> aRecoCand = RefitPVNoBS->at(jrefit).userCand( name );
	      size_t hash = hasher(aRecoCand.get());
	      hashesNoBS.push_back(hash);
	    }
	  _VertexHashNoBS1.push_back(hashesNoBS[0]);
	  _VertexHashNoBS2.push_back(hashesNoBS[1]);
	}
    }

  return RefitPVNoBSHandle.isValid();
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

bool HTauTauNtuplizer::findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  Handle<vector<reco::Vertex> >  vertices;
  iEvent.getByToken(theVtxTag,vertices);

  if(vertices->size()==0) return false;   //at least one vertex

  _pv_x = (*vertices)[0].x();
  _pv_y = (*vertices)[0].y();
  _pv_z = (*vertices)[0].z();
  _isRefitPV = refitPV(iEvent, iSetup);
  return true;
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
TVector3 HTauTauNtuplizer::getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
				  const reco::Track *aTrack,
				  const GlobalPoint & aPoint){
  TVector3 aPCA;
  if(!doCPVariables || !aTrack ||  _npv==0 || aTrack->pt()<2) return aPCA;

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  if(!transTrackBuilder.isValid()){
    std::cout<<"Problem with TransientTrackBuilder"<<std::endl;
    return aPCA;
  }

  reco::TransientTrack transTrk=transTrackBuilder->build(aTrack);

  AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
  GlobalPoint pos  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).globalPosition();

  aPCA.SetX(pos.x() - aPoint.x());
  aPCA.SetY(pos.y() - aPoint.y());
  aPCA.SetZ(pos.z() - aPoint.z());

  return aPCA;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


// // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void HTauTauNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

//define this as a plug-in
DEFINE_FWK_MODULE(HTauTauNtuplizer);
