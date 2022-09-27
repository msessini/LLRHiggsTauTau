#include <LLRHiggsTauTau/NtupleProducer/interface/FakeFactors.h>

FakeFactors::FakeFactors(Int_t theYear, TH2D* ff_fracs_qcd, TH2D* ff_fracs_wjets, TH2D* ff_fracs_qcd_ss, TH2D* ff_fracs_wjets_ss, TH2D* ff_fracs_qcd_aiso, TH2D* ff_fracs_wjets_aiso, TH2D* ff_fracs_qcd_highmt, TH2D* ff_fracs_wjets_highmt, std::shared_ptr<RooWorkspace> ff_ws_) {

  auto year_ = std::to_string(theYear);
  ff_fracs_qcd_ = ff_fracs_qcd;
  ff_fracs_wjets_ = ff_fracs_wjets;
  ff_fracs_qcd_ss_ = ff_fracs_qcd_ss;
  ff_fracs_wjets_ss_ = ff_fracs_wjets_ss;
  ff_fracs_qcd_aiso_ = ff_fracs_qcd_aiso;
  ff_fracs_wjets_aiso_ = ff_fracs_wjets_aiso;
  ff_fracs_qcd_highmt_ = ff_fracs_qcd_highmt;
  ff_fracs_wjets_highmt_ = ff_fracs_wjets_highmt;

  systs_mvadm_ = {"","_wjets_syst_up","_wjets_syst_down","_wjets_met_up","_wjets_met_down","_wjets_l_pt_up","_wjets_l_pt_down","_wjets_stat_unc1_njet0_mvadm10_up","_wjets_stat_unc2_njet0_mvadm10_up","_wjets_stat_unc1_njet0_mvadm10_down","_wjets_stat_unc2_njet0_mvadm10_down","_wjets_stat_unc1_njet1_mvadm10_up","_wjets_stat_unc2_njet1_mvadm10_up","_wjets_stat_unc1_njet1_mvadm10_down","_wjets_stat_unc2_njet1_mvadm10_down","_wjets_stat_unc1_njet2_mvadm10_up","_wjets_stat_unc2_njet2_mvadm10_up","_wjets_stat_unc1_njet2_mvadm10_down","_wjets_stat_unc2_njet2_mvadm10_down","_qcd_syst_up","_qcd_syst_down","_qcd_met_up","_qcd_met_down","_qcd_l_pt_up","_qcd_l_pt_down","_qcd_stat_unc1_njet0_mvadm10_up","_qcd_stat_unc2_njet0_mvadm10_up","_qcd_stat_unc1_njet0_mvadm10_down","_qcd_stat_unc2_njet0_mvadm10_down","_qcd_stat_unc1_njet1_mvadm10_up","_qcd_stat_unc2_njet1_mvadm10_up","_qcd_stat_unc1_njet1_mvadm10_down","_qcd_stat_unc2_njet1_mvadm10_down","_qcd_stat_unc1_njet2_mvadm10_up","_qcd_stat_unc2_njet2_mvadm10_up","_qcd_stat_unc1_njet2_mvadm10_down","_qcd_stat_unc2_njet2_mvadm10_down","_ttbar_syst_up","_ttbar_syst_down","_ttbar_met_up","_ttbar_met_down"};

  for(auto s : systs_mvadm_) {
    fns_["ff_lt_medium_mvadmbins"+s] = std::shared_ptr<RooFunctor>(
								   ff_ws_->function(("ff_mt_medium_mvadmbins"+s).c_str())->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,os,met_var_qcd,met_var_w,mt,m_iso,pass_single,mvis,WpT,wjets_frac,qcd_frac,ttbar_frac")));
  }

  fns_["ff_lt_medium_mvadmbins_qcd"] = std::shared_ptr<RooFunctor>(
            ff_ws_->function("ff_mt_medium_mvadmbins_qcd")->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,os,met_var_qcd,m_iso,pass_single")));
  fns_["ff_lt_medium_mvadmbins_wjets"] = std::shared_ptr<RooFunctor>(
            ff_ws_->function("ff_mt_medium_mvadmbins_wjets")->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,met_var_w,mt,pass_single,mvis,WpT")));
  //fns_["ff_lt_medium_mvadmbins_ttbar"] = std::shared_ptr<RooFunctor>(
  //          ff_ws_->function("ff_mt_medium_mvadmbins_ttbar")->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,met_var_w")));

}

void FakeFactors::Initialize(LorentzVectorParticle TauLVP, const reco::Candidate *mu, int tauDM, float met_px, float met_py, int Njets, double dijetMass, TVector3 pv, std::vector<std::vector<double>> pvcov) {

  math::XYZTLorentzVector taup4(TauLVP.LV().Px(), TauLVP.LV().Py(), TauLVP.LV().Pz(), TauLVP.LV().E());
  math::XYZTLorentzVector mup4 = mu->p4();
  math::XYZTLorentzVector metp4(met_px, met_py, 0., 0.);
  math::XYZTLorentzVector metp4w(metp4.px() + mup4.px(), metp4.py() + mup4.py(), 0., mup4.pt());
  //
  pt_tt_ = sqrt(pow((taup4.px() + mup4.px() + metp4.px()),2) + pow((taup4.py() + mup4.py() + metp4.py()),2));
  pt_1_ = mup4.pt();
  pt_2_ = taup4.pt();
  met_ = metp4.pt();
  m_vis_ = (taup4 + mup4).mass();
  n_jets_ = Njets;
  mjj_ = dijetMass;
  mva_dm_2_ = tauDM;
  mt_1_ = seltools::ComputeMT(mup4, metp4.px(), metp4.py());
  //
  double iso_1_ = userdatahelpers::getUserFloat(mu,"combRelIsoPF");
  double met_var_qcd = (metp4.pt()/taup4.pt())*cos(deltaPhi(metp4.phi(), taup4.phi()));
  double met_var_w = (metp4w.pt()/taup4.pt())*cos(deltaPhi(metp4w.phi(), taup4.phi()));
  double WpT = metp4w.pt();
  //
  double singlemupt = 25.;
  if(year_ == "2016") singlemupt = 23.;
  double pass_single = 1.;
  if(mup4.pt() < singlemupt) pass_single = 0.;
  //
  bool isOS = (TauLVP.Charge()*mu->charge())<0 ? true : false;
  // load MVA scroes reader for fractions
  reader_ = new TMVA::Reader();
  reader_->AddVariable("pt_tt", &pt_tt_);
  reader_->AddVariable("pt_1", &pt_1_);
  reader_->AddVariable("pt_2", &pt_2_);
  reader_->AddVariable("met", &met_);
  reader_->AddVariable("m_vis", &m_vis_);
  reader_->AddVariable("n_jets", &n_jets_);
  reader_->AddVariable("mjj", &mjj_);
  reader_->AddVariable("mva_dm_2", &mva_dm_2_);
  reader_->AddVariable("mt_1", &mt_1_);
  xml_file="src/LLRHiggsTauTau/NtupleProducer/data/fake_factors/fractions_2018_mt.xml";
  //xml_file="/opt/sbg/cms/safe1/cms/msessini/MuTauProducer/CMSSW_10_2_23/src/LLRHiggsTauTau/NtupleProducer/data/fake_factors/fractions_2018_mt.xml";
  reader_->BookMVA("BDT method", xml_file);
  //
  std::vector<float> scores = reader_->EvaluateMulticlass("BDT method");
  double qcd_score = scores[1];
  double w_score = scores[0];
  //
  double w_frac = ff_fracs_wjets_->GetBinContent(ff_fracs_wjets_->FindBin(qcd_score,w_score));
  double qcd_frac = ff_fracs_qcd_->GetBinContent(ff_fracs_qcd_->FindBin(qcd_score,w_score));
  //
  if(!isOS) {
    w_frac = ff_fracs_wjets_ss_->GetBinContent(ff_fracs_wjets_ss_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_ss_->GetBinContent(ff_fracs_qcd_ss_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) qcd_frac = 1.;
  }
  if(iso_1_>0.15) {
    w_frac = ff_fracs_wjets_aiso_->GetBinContent(ff_fracs_wjets_aiso_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_aiso_->GetBinContent(ff_fracs_qcd_aiso_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) qcd_frac = 1.;
  }
  if(mt_1_>70) {
    w_frac = ff_fracs_wjets_highmt_->GetBinContent(ff_fracs_wjets_highmt_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_highmt_->GetBinContent(ff_fracs_qcd_highmt_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) w_frac = 1.;
  }
  double ttbar_frac = 1. - w_frac - qcd_frac;
  double os = 1.;
  if(!isOS) os = 0.;
  //
  double ipsig = GetIPsig(TauLVP, pv, pvcov);
  //
  args_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,os,met_var_qcd,met_var_w,mt_1_,iso_1_,pass_single,m_vis_,WpT,w_frac,qcd_frac,ttbar_frac};
  args_qcd_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,os,met_var_qcd,iso_1_,pass_single};
  args_w_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,met_var_w,mt_1_,pass_single,m_vis_,WpT};
  args_ttbar_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,met_var_w};
}

std::map<std::string, double> FakeFactors::GetFakeFactors(std::string sysType) {

  fake_factors_["ff_nominal"] = fns_["ff_lt_medium_mvadmbins"]->eval(args_.data());
  fake_factors_["ff_nominal_qcd"] = fns_["ff_lt_medium_mvadmbins_qcd"]->eval(args_qcd_.data()); 
  fake_factors_["ff_nominal_w"] = fns_["ff_lt_medium_mvadmbins_wjets"]->eval(args_w_.data());
  //fake_factors_["ff_nominal_ttbar"] = fns_["ff_lt_medium_mvadmbins_ttbar"]->eval(args_ttbar_.data());
  //
  if(sysType == "Nominal") {
    for(auto s : systs_mvadm_) {
      if(s == "") continue;
      fake_factors_["ff"+s] = fns_["ff_lt_medium_mvadmbins"+s]->eval(args_.data());
    }
  }
  return fake_factors_;
}

//////////////////////////////////////////

double FakeFactors::GetIPsig(LorentzVectorParticle TauLVP, TVector3 pv, std::vector<std::vector<double>> pvcov) {

  TrackParticle TauTrack(TauLVP.getParMatrix(), TauLVP.getCovMatrix(), TauLVP.PDGID(), TauLVP.Mass(), TauLVP.Charge(), TauLVP.BField());
 
  std::vector<float> h_param = {float(TauTrack.Parameter(TrackParticle::kappa)), 
				float(TauTrack.Parameter(TrackParticle::lambda)),
                                float(TauTrack.Parameter(TrackParticle::phi)),
                                float(TauTrack.Parameter(TrackParticle::dxy)),
                                float(TauTrack.Parameter(TrackParticle::dz))};

  ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(TauLVP.Parameter(0),TauLVP.Parameter(1),TauLVP.Parameter(2));
  ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pvtx(pv.X(),pv.Y(),pv.Z());

  ImpactParameter IP;
  TVector3 ip = IP.CalculatePCA(TauTrack.BField(), h_param, ref, pvtx);

  ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov;
  TMatrixTSym<double> cov = TauLVP.getCovMatrix();
  SMatrixSym3D SigmaPrV;
  for(int i=0; i<5; i++) {
    for(int j=0; j<5; j++) {
      helixCov(i,j) = cov(i,j);
      if(i<3 && j<3) SigmaPrV(i,j) = pvcov[i][j];
    }
  }

  ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

  double mag = ip.Mag();
  ROOT::Math::SVector<double, 3> ip_svec;
  ip_svec(0) = ip.X();
  ip_svec(1) = ip.Y();
  ip_svec(2) = ip.Z();

  ip_svec = ip_svec.Unit();

  double uncert = sqrt(ROOT::Math::Dot( ip_svec, ip_cov * ip_svec));
  double sig = mag/uncert;
  return sig;
}
