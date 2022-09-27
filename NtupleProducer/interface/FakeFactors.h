#ifndef FakeFactors_h
#define FakeFactors_h

#include <TH2D.h>
#include <string>
#include <TMVA/Reader.h>
#include <RooFunctor.h>
#include <RooWorkspace.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/SelectionTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/DaughterDataHelpers.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/LorentzVectorParticle.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/TrackParticle.h>
#include <HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h>

class FakeFactors {
 public:
  FakeFactors(Int_t theYear, TH2D* ff_fracs_qcd, TH2D* ff_fracs_wjets, TH2D* ff_fracs_qcd_ss, TH2D* ff_fracs_wjets_ss, TH2D* ff_fracs_qcd_aiso, TH2D* ff_fracs_wjets_aiso, TH2D* ff_fracs_qcd_highmt, TH2D* ff_fracs_wjets_highmt, std::shared_ptr<RooWorkspace> ff_ws_);
  ~FakeFactors(){};
  void Initialize(LorentzVectorParticle TauLVP, const reco::Candidate *mu, int tauDM, float met_px, float met_py, int Njets, double dijetMass, TVector3 pv, std::vector<std::vector<double>> pvcov);
  std::map<std::string, double> GetFakeFactors(std::string sysType);

 private:
  std::shared_ptr<RooWorkspace> w_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  std::map<std::string, double> fake_factors_;

  TH2D *ff_fracs_qcd_;
  TH2D *ff_fracs_wjets_;
  TH2D *ff_fracs_qcd_ss_;
  TH2D *ff_fracs_wjets_ss_;
  TH2D *ff_fracs_qcd_aiso_;
  TH2D *ff_fracs_wjets_aiso_;
  TH2D *ff_fracs_qcd_highmt_;
  TH2D *ff_fracs_wjets_highmt_;

  TString xml_file, year_;
  TMVA::Reader *reader_;

  float met_, pt_1_, pt_2_, mva_dm_2_, mt_1_, m_vis_, pt_tt_, mjj_, n_jets_;

  std::vector<std::string> systs_mvadm_;
  std::vector<double> args_, args_qcd_, args_w_, args_ttbar_;
  //
  double GetIPsig(LorentzVectorParticle TauLVP, TVector3 pv, std::vector<std::vector<double>> pvcov);
};

#endif
