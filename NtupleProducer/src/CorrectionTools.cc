#include <LLRHiggsTauTau/NtupleProducer/interface/SelectionTools.h>
#include <LLRHiggsTauTau/NtupleProducer/interface/CorrectionTools.h>

void corrector::IPCorrection(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, TVector3 IP, double Eta, IpCorrection ipCorrector) {

  edm::Handle<edm::View<pat::GenericParticle> > candHandle;
  event.getByToken(generictag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();
  TVector3 muVtx(0., 0., 0.);
  TVector3 tauVtx(0., 0., 0.);
  TVector3 p(0., 0., 0.);
  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen) {
    unsigned pdgid = abs(igen->pdgId());
    unsigned status = abs(igen->status());
    unsigned flag = igen->userInt ("generalGenFlags");
    if(pdgid == 13 && status == 1){
      muVtx.SetXYZ(igen->vx(), igen->vy(), igen->vz());
      p.SetXYZ(igen->px(), igen->py(), igen->pz());
    }
    if(pdgid == 15 && seltools::CHECK_BIT(flag,0) && status == 2) tauVtx.SetXYZ(igen->vx(), igen->vy(), igen->vz());
  }
  TVector3 k = muVtx - tauVtx;
  TVector3 genIP = k - (p.Dot(k) / p.Mag2()) * p;
  //
  TVector3 uncorrIP = IP;
  IP = ipCorrector.correctIp(IP, genIP, Eta);
}

void corrector::METRecoilCorrection(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, const pat::MET& PUPPImet, int Njets, float shiftedMETx, float shiftedMETy, std::string sysType, std::string var, RecoilCorrector* recoilPuppiMetCorrector, MEtSys* recoilPuppiMetShifter) {

  edm::Handle<edm::View<pat::GenericParticle> > candHandle;
  event.getByToken(generictag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();
  math::XYZTLorentzVector visgenP4, genP4;
  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen) {
    unsigned pdgid = abs(igen->pdgId());
    unsigned status = abs(igen->status());
    unsigned flag = igen->userInt ("generalGenFlags");
    math::XYZTLorentzVector iGenP4(igen->px(), igen->py(), igen->pz(), igen->energy());
    if((pdgid>=11 && pdgid<=16 && seltools::CHECK_BIT(flag,8) && status == 1) || seltools::CHECK_BIT(flag,10)) genP4+=iGenP4;
    if(((pdgid == 11 || pdgid == 13 || pdgid == 15) && seltools::CHECK_BIT(flag,8) && status==1) || (seltools::CHECK_BIT(flag,10) && !(pdgid==12||pdgid==14||pdgid==16))) visgenP4+=iGenP4;
  }
  //
  if(sysType == "METReso") {
    if(var == "Up") {
      recoilPuppiMetShifter->ApplyMEtSys(
					 (float)PUPPImet.px(), // uncorrected type I puppi met px (float)
					 (float)PUPPImet.py(), // uncorrected type I puppi met py (float)
					 (float)genP4.px(), // generator Z/W/Higgs px (float)
					 (float)genP4.py(), // generator Z/W/Higgs py (float)
					 (float)visgenP4.px(), // generator visible Z/W/Higgs px (float)
					 (float)visgenP4.py(), // generator visible Z/W/Higgs py (float)
					 Njets,  // number of jets (hadronic jet multiplicity) (int)
					 1, // shift for hadronic recoil response
					 0, // upward shift
					 shiftedMETx, // shifted type I puppi met px (float)
					 shiftedMETy  // shifted type I puppi met py (float)
					 );
    }
    if(var == "Down") {
      recoilPuppiMetShifter->ApplyMEtSys(
					 (float)PUPPImet.px(), // uncorrected type I puppi met px (float)
					 (float)PUPPImet.py(), // uncorrected type I puppi met py (float)
					 (float)genP4.px(), // generator Z/W/Higgs px (float)
					 (float)genP4.py(), // generator Z/W/Higgs py (float)
					 (float)visgenP4.px(), // generator visible Z/W/Higgs px (float)
					 (float)visgenP4.py(), // generator visible Z/W/Higgs py (float)
					 Njets,  // number of jets (hadronic jet multiplicity) (int)
					 1, // shift for hadronic recoil response
					 1, // downward shift
					 shiftedMETx, // shifted type I puppi met px (float)
					 shiftedMETy  // shifted type I puppi met py (float)
					 );
    }
  }
  else if(sysType == "METScale") {
    if(var == "Up") {
      recoilPuppiMetShifter->ApplyMEtSys(
					 (float)PUPPImet.px(), // uncorrected type I puppi met px (float)
					 (float)PUPPImet.py(), // uncorrected type I puppi met py (float)
					 (float)genP4.px(), // generator Z/W/Higgs px (float)
					 (float)genP4.py(), // generator Z/W/Higgs py (float)
					 (float)visgenP4.px(), // generator visible Z/W/Higgs px (float)
					 (float)visgenP4.py(), // generator visible Z/W/Higgs py (float)
					 Njets,  // number of jets (hadronic jet multiplicity) (int)
					 0, // shift for hadronic recoil response
					 0, // upward shift
					 shiftedMETx, // shifted type I puppi met px (float)
					 shiftedMETy  // shifted type I puppi met py (float)
					 );
    }
    if(var == "Down") {
      recoilPuppiMetShifter->ApplyMEtSys(
					 (float)PUPPImet.px(), // uncorrected type I puppi met px (float)
					 (float)PUPPImet.py(), // uncorrected type I puppi met py (float)
					 (float)genP4.px(), // generator Z/W/Higgs px (float)
					 (float)genP4.py(), // generator Z/W/Higgs py (float)
					 (float)visgenP4.px(), // generator visible Z/W/Higgs px (float)
					 (float)visgenP4.py(), // generator visible Z/W/Higgs py (float)
					 Njets,  // number of jets (hadronic jet multiplicity) (int)
					 0, // shift for hadronic recoil response
					 1, // downward shift
					 shiftedMETx, // shifted type I puppi met px (float)
					 shiftedMETy  // shifted type I puppi met py (float)
					 );
    }
  }
  else {
    recoilPuppiMetCorrector->CorrectWithHist(
					     (float)PUPPImet.px(), // uncorrected type I pf met px (float)
					     (float)PUPPImet.py(), // uncorrected type I pf met py (float)
					     (float)genP4.px(), // generator Z/W/Higgs px (float)
					     (float)genP4.py(), // generator Z/W/Higgs py (float)
					     (float)visgenP4.px(), // generator visible Z/W/Higgs px (float)
					     (float)visgenP4.py(), // generator visible Z/W/Higgs py (float)
					     Njets,  // number of jets (hadronic jet multiplicity) (int)
					     shiftedMETx, // corrected type I pf met px (float)
					     shiftedMETy  // corrected type I pf met py (float)
					     );
  }
}

math::XYZTLorentzVector corrector::P4Corrected(math::XYZTLorentzVector p4, int genmatch, int DM, std::string Unc, TH1* histTES, TGraph* histFES){

  double ShiftP = 1.;
  double ShiftM = 1.;
  //
  if(genmatch == 5) { // For genuine taus

    double NomTES=0;
    double UncTES=0;
    //
    Int_t binTES = histTES->GetXaxis()->FindBin(DM);
    NomTES = histTES->GetBinContent(binTES);
    UncTES = histTES->GetBinError(binTES);

    if(Unc=="Nom") {
      ShiftP = NomTES;
      ShiftM = NomTES;
    }
    else if(Unc=="Up") {
      ShiftP = NomTES + UncTES;
      ShiftM = NomTES + UncTES;
    }
    else if(Unc=="Down") {
      ShiftP = NomTES - UncTES;
      ShiftM = NomTES - UncTES;
    }

    if(DM==0) ShiftM = 1.; //Pion mass is not shifted
  }

  else if(genmatch == 1 || genmatch == 3) { // For genuine electrons

    double NomFES=0;
    double UncFES=0;
    //
    if(p4.eta() <= 1.479) { //Barrel
      if(DM==0) {
        NomFES = histFES->GetY()[0];
        if(Unc=="Up") UncFES = std::abs(histFES->GetErrorYhigh(0));
        if(Unc=="Down") UncFES = -std::abs(histFES->GetErrorYlow(0));
      }
      if(DM==1) {
        NomFES = histFES->GetY()[1];
        if(Unc=="Up") UncFES = std::abs(histFES->GetErrorYhigh(1));
        if(Unc=="Down") UncFES = -std::abs(histFES->GetErrorYlow(1));
      }
    }
    if(p4.eta() > 1.479) { //Endcap
      if(DM==0) {
        NomFES = histFES->GetY()[2];
        if(Unc=="Up") UncFES = std::abs(histFES->GetErrorYhigh(2));
        if(Unc=="Down") UncFES = -std::abs(histFES->GetErrorYlow(2));
      }
      if(DM==1) {
        NomFES = histFES->GetY()[3];
        if(Unc=="Up") UncFES = std::abs(histFES->GetErrorYhigh(3));
        if(Unc=="Down") UncFES = -std::abs(histFES->GetErrorYlow(3));
      }
    }

    ShiftP = NomFES;
    ShiftM = NomFES;
    if(Unc == "Up" || Unc == "Down") {
      ShiftP += UncFES;
      ShiftM += UncFES;
    }

    if(DM==1) ShiftM = 1.; //Pion mass is not shifted
  }

  if(genmatch==2 || genmatch==4) { // For genuine muons

    if(Unc=="Up") {
      ShiftP = 1.01;
      ShiftM = 1.01;
    }
    if(Unc=="Down") {
      ShiftP = 0.99;
      ShiftM = 0.99;
    }
    if(DM==0) ShiftM = 1.; //Pion mass is not shifted
  }

  double px_scaled = p4.Px()*ShiftP;
  double py_scaled = p4.Py()*ShiftP;
  double pz_scaled = p4.Pz()*ShiftP;
  double mass_scaled = p4.M()*ShiftM;
  double en_scaled = std::sqrt(px_scaled*px_scaled + py_scaled*py_scaled + pz_scaled*pz_scaled + mass_scaled*mass_scaled);

  return math::XYZTLorentzVector(px_scaled, py_scaled, pz_scaled, en_scaled);
}

math::XYZTLorentzVector corrector::MuP4Corrected(math::XYZTLorentzVector p4, std::string Unc)
{
  double scale = 0.;
  double px_scaled = p4.Px();
  double py_scaled = p4.Py();
  double pz_scaled = p4.Pz();
  double m_scaled = p4.M();
  double E_scaled = p4.energy();

  if(abs(p4.eta())<1.2) scale = 0.004;                            //barrel
  else if(abs(p4.eta())>1.2 && abs(p4.eta())<2.1) scale = 0.009;  //near endcap
  else if(abs(p4.eta())>2.1 && abs(p4.eta())<2.4) scale = 0.027;  //far endcap

  if(Unc=="Up"){
    px_scaled *= (1+scale);
    py_scaled *= (1+scale);
    pz_scaled *= (1+scale);
    m_scaled *= (1+scale);
    E_scaled = std::sqrt(px_scaled*px_scaled + py_scaled*py_scaled + pz_scaled*pz_scaled + m_scaled*m_scaled);
  }
  else if(Unc=="Down"){
    px_scaled *= (1-scale);
    py_scaled *= (1-scale);
    pz_scaled *= (1-scale);
    m_scaled *= (1-scale);
    E_scaled = std::sqrt(px_scaled*px_scaled + py_scaled*py_scaled + pz_scaled*pz_scaled + m_scaled*m_scaled);
  }

  return math::XYZTLorentzVector(px_scaled, py_scaled, pz_scaled, E_scaled);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

std::map<std::string, double> weight::ZpTreweighting(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, std::string sysType, std::shared_ptr<RooWorkspace> w) {

  std::map<std::string, double> ZpTreweightingmap;
  //
  double zptw = 1.;
  edm::Handle<edm::View<pat::GenericParticle> > candHandle;

  event.getByToken(generictag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();
  math::XYZTLorentzVector genMomentum(0.,0.,0.,0.);
  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen) {
    unsigned pdgid = abs(igen->pdgId());
    unsigned status = abs(igen->status());
    unsigned flag = igen->userInt ("generalGenFlags");
    //
    if((pdgid == 15 && seltools::CHECK_BIT(flag,0) && status == 2) || ((pdgid == 11 || pdgid == 13) && status == 1)) {
      if(igen->pt() > 8) {
	genMomentum+=igen->p4();
      }
    }
  }
  //
  if(genMomentum.pt() != 0 && genMomentum.M() > 75 && genMomentum.M() < 120){
    auto arg = std::vector<double>{genMomentum.pt(),genMomentum.M()};
    zptw = std::shared_ptr<RooFunctor>(w->function("zptmass_weight_nom")->functor(w->argSet("z_gen_pt,z_gen_mass")))->eval(arg.data());
  }
  ZpTreweightingmap["wZpT"] = zptw;
  //
  if(sysType == "Nominal") {
    ZpTreweightingmap["wZpTUp"] = zptw*1.1;
    ZpTreweightingmap["wZpTDown"] = zptw*0.9;
  }
  //
  return ZpTreweightingmap;
}

std::map<std::string, double> weight::ToppTreweighting(const edm::Event& event, edm::EDGetTokenT<edm::View<pat::GenericParticle>> generictag, std::string sysType) {

  std::map<std::string, double> ToppTreweightingmap;
  //
  double topw = 1.;
  edm::Handle<edm::View<pat::GenericParticle> > candHandle;

  event.getByToken(generictag, candHandle);
  const edm::View<pat::GenericParticle>* gens = candHandle.product();
  math::XYZTLorentzVector genMomentum(0.,0.,0.,0.);
  for(edm::View<pat::GenericParticle>::const_iterator igen = gens->begin(); igen!=gens->end(); ++igen) {
    unsigned pdgid = abs(igen->pdgId());
    unsigned flag = igen->userInt ("generalGenFlags");
    //
    if(pdgid == 6 && flag == 8 && flag == 13){
      double pt = igen->pt();
      pt = std::min(pt, 472.);
      double a = 0.088, b = -0.00087, c = 9.2e-07;
      topw *= std::exp(a + b * pt + c * pt*pt);
    }
  }
  topw = std::sqrt(topw);
  ToppTreweightingmap["wToppT"] = topw;
  //
  if(sysType == "Nominal") {
    ToppTreweightingmap["wToppTUp"] = topw*topw;
    ToppTreweightingmap["wToppTDown"] = 1;
  }
  return ToppTreweightingmap;
}

std::map<std::string, double> weight::PrefiringWeight(double pref, double prefUp, double prefDown, Int_t theYear, std::string sysType) {

  std::map<std::string, double> PrefiringWeightmap;
  //
  double prefw = 1.;
  if(theYear==2016 || theYear==2017) {
    prefw = pref;
    if(sysType == "Nominal") {
      double prefwUp = prefUp;
      double prefwDown = prefDown;
      //
      if((prefwUp/prefw)<1.2) prefwUp = prefw*1.2;
      if((prefwDown/prefw)>0.8) prefwDown = prefw*0.8;
      PrefiringWeightmap["wPrefiringUp"] = prefwUp;
      PrefiringWeightmap["wPrefiringDown"] = prefwDown;
    }
  }
  PrefiringWeightmap["wPrefiringWT"] = prefw;
  //
  return PrefiringWeightmap;
}

std::map<std::string, double> weight::TauIDSF(int genmatch, float DM, math::XYZTLorentzVector p4, bool isEmbed, std::string sysType, std::shared_ptr<RooWorkspace> w, std::string Label) {

  std::map<std::string, double> TauIDSFmap;
  //
  auto decayMode = std::to_string(int(DM));
  //
  double t_pt=p4.Pt();
  double t_mvadm=DM;
  auto pt_mvadm = std::vector<double>{t_pt,t_mvadm};
  //
  double wIDvsJet = 1;
  double wIDvsEle = 1;
  double wIDvsMu = 1;
  //
  if(genmatch == 5) {
    if(!isEmbed) {
      wIDvsJet = std::shared_ptr<RooFunctor>(w->function("t_deeptauid_mvadm_medium")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
    }
    else wIDvsJet = std::shared_ptr<RooFunctor>(w->function("t_deeptauid_mvadm_embed_medium")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
  }
  else if(genmatch == 1 || genmatch == 3) {
    TauIDSFTool *tauSFvsEle=new TauIDSFTool(Label,"DeepTau2017v2p1VSe","VVLoose");
    wIDvsEle = tauSFvsEle->getSFvsEta(abs(p4.Eta()),genmatch);
    delete tauSFvsEle;
  }
  else if(genmatch == 2 || genmatch == 4) {
    TauIDSFTool *tauSFvsMu=new TauIDSFTool(Label,"DeepTau2017v2p1VSmu","Tight");
    wIDvsMu = tauSFvsMu->getSFvsEta(abs(p4.Eta()),genmatch);
    delete tauSFvsMu;
  }
  TauIDSFmap["wIDvsJet"] = wIDvsJet;
  TauIDSFmap["wIDvsEle"] = wIDvsEle;
  TauIDSFmap["wIDvsMu"] = wIDvsMu;

  if(sysType == "Nominal") {
    double wIDvsJetUp = 1;
    double wIDvsJetDown = 1;
    double wIDvsEleUp = 1;
    double wIDvsEleDown = 1;
    double wIDvsMuUp = 1;
    double wIDvsMuDown = 1;
    //
    if(genmatch == 5) {
      std::string ptlevel = (t_pt<40) ? "lowpt" : "highpt";
      if(!isEmbed) {
        wIDvsJetUp = std::shared_ptr<RooFunctor>(w->function(("t_deeptauid_mvadm_medium_"+ptlevel+"_mvadm"+decayMode+"_up").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
	wIDvsJetDown = std::shared_ptr<RooFunctor>(w->function(("t_deeptauid_mvadm_medium_"+ptlevel+"_mvadm"+decayMode+"_down").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
      }
      else {
        wIDvsJetUp = std::shared_ptr<RooFunctor>(w->function(("t_deeptauid_mvadm_embed_medium_"+ptlevel+"_mvadm"+decayMode+"_up").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
        wIDvsJetDown = std::shared_ptr<RooFunctor>(w->function(("t_deeptauid_mvadm_embed_medium_"+ptlevel+"_mvadm"+decayMode+"_down").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
      }
    }
    else if(genmatch == 1 || genmatch == 3) {
      TauIDSFTool *tauSFvsEle=new TauIDSFTool(Label,"DeepTau2017v2p1VSe","VVLoose");
      wIDvsEleUp = tauSFvsEle->getSFvsEta(abs(p4.Eta()),genmatch,"Up");
      wIDvsEleDown = tauSFvsEle->getSFvsEta(abs(p4.Eta()),genmatch,"Down");
      delete tauSFvsEle;
    }
    else if(genmatch == 2 || genmatch == 4) {
      TauIDSFTool *tauSFvsMu=new TauIDSFTool(Label,"DeepTau2017v2p1VSmu","Tight");
      wIDvsMuUp = tauSFvsMu->getSFvsEta(abs(p4.Eta()),genmatch,"Up");
      wIDvsMuDown = tauSFvsMu->getSFvsEta(abs(p4.Eta()),genmatch,"Down");
      delete tauSFvsMu;
    }
    TauIDSFmap["wIDvsJetUp"] = wIDvsJetUp;
    TauIDSFmap["wIDvsJetDown"] = wIDvsJetDown;
    TauIDSFmap["wIDvsEleUp"] = wIDvsEleUp;
    TauIDSFmap["wIDvsEleDown"] = wIDvsEleDown;
    TauIDSFmap["wIDvsMuUp"] = wIDvsMuUp;
    TauIDSFmap["wIDvsMuDOwn"] = wIDvsMuDown;
  }
  return TauIDSFmap;
}

std::map<std::string, double> weight::MuonIDTrkSF(math::XYZTLorentzVector p4, double m_iso, bool isEmbed, std::shared_ptr<RooWorkspace> w) {

  std::map<std::string, double> MuonIDTrkSFmap;
  //
  double m_pt = p4.Pt();
  double m_eta = p4.Eta();
  auto pt_eta_iso = std::vector<double>{m_pt,m_eta,m_iso};
  auto eta = std::vector<double>{m_eta};
  //
  double wID = 1;
  double wTrk = 1;
  //
  if(!isEmbed) {
    wID = std::shared_ptr<RooFunctor>(w->function("m_idiso_binned_ic_ratio")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  }
  else wID = std::shared_ptr<RooFunctor>(w->function("m_idiso_binned_ic_embed_ratio")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  //
  wTrk = std::shared_ptr<RooFunctor>(w->function("m_trk_ratio")->functor(w->argSet("m_eta")))->eval(eta.data());
  //
  MuonIDTrkSFmap["wID"] = wID;
  MuonIDTrkSFmap["wTrk"] = wTrk;
  //
  return MuonIDTrkSFmap;
}

std::map<std::string, double> weight::TriggerSF(int genmatch, float DM, math::XYZTLorentzVector taup4, double m_iso, math::XYZTLorentzVector mup4, Int_t theYear, bool isEmbed, std::string sysType, std::shared_ptr<RooWorkspace> w) {

  std::map<std::string, double> TriggerSFmap;
  //
  auto year = std::to_string(theYear);
  auto decayMode = std::to_string(int(DM));
  //
  double t_pt=taup4.Pt();
  double t_mvadm=DM;
  auto pt_mvadm = std::vector<double>{t_pt,t_mvadm};
  //
  double m_pt = mup4.Pt();
  double m_eta = mup4.Eta();
  auto pt_eta_iso = std::vector<double>{m_pt,m_eta,m_iso};
  //
  bool m_high_pT = (m_pt > 25.);
  if(theYear == 2016) m_high_pT = (m_pt > 23.);
  //
  double m_trg = 1;
  double m_trg_mc = 1;
  //
  m_trg = std::shared_ptr<RooFunctor>(w->function("m_trg_binned_ic_data")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  if(!isEmbed) {
    m_trg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_binned_ic_mc")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  }
  else m_trg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_binned_ic_embed")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  //
  double m_xtrg = 1;
  double m_xtrg_mc =1;
  if(theYear == 2016) {
    m_xtrg = std::shared_ptr<RooFunctor>(w->function("m_trg_19_binned_ic_data")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
    if(!isEmbed) {
      m_xtrg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_19_binned_ic_mc")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
    }
    else m_xtrg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_19_binned_ic_embed")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  }
  if (theYear == 2017 || theYear == 2018) {
    m_xtrg = std::shared_ptr<RooFunctor>(w->function("m_trg_20_binned_ic_data")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
    if(!isEmbed) {
      m_xtrg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_20_binned_ic_mc")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
    }
    else m_xtrg_mc = std::shared_ptr<RooFunctor>(w->function("m_trg_20_binned_ic_embed")->functor(w->argSet("m_pt,m_eta,m_iso")))->eval(pt_eta_iso.data());
  }
  //
  double t_trg = 1;
  if(!isEmbed) {
    t_trg = std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
  }
  else t_trg = std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
  //
  double single_m_sf = m_trg/m_trg_mc;
  double xtrg_mt_sf = (m_xtrg_mc>0) ? (m_xtrg*t_trg)/m_xtrg_mc : 0.0;
  //
  double xtrg_OR_sf = m_high_pT ? single_m_sf : xtrg_mt_sf;
  TriggerSFmap["wTrg"] = xtrg_OR_sf;

  if(sysType == "Nominal") {
    double t_trg_Up = 1;
    double t_trg_Down = 1;
    if(!m_high_pT) { // otherwise use single trigger weight which is 1 !
      if(!isEmbed && std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) > 0) {
	t_trg_Up = std::shared_ptr<RooFunctor>(w->function(("t_trg_ic_deeptau_medium_mvadm_mutau_ratio_mvadm"+decayMode+"_up").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) /
	  std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
	t_trg_Down = std::shared_ptr<RooFunctor>(w->function(("t_trg_ic_deeptau_medium_mvadm_mutau_ratio_mvadm"+decayMode+"_down").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) /
	  std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
      }
      else if (std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) > 0){
	t_trg_Up = std::shared_ptr<RooFunctor>(w->function(("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio_mvadm"+decayMode+"_up").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) /
	  std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
        t_trg_Down = std::shared_ptr<RooFunctor>(w->function(("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio_mvadm"+decayMode+"_down").c_str())->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data()) /
	  std::shared_ptr<RooFunctor>(w->function("t_trg_ic_deeptau_medium_mvadm_mutau_embed_ratio")->functor(w->argSet("t_pt,t_mvadm")))->eval(pt_mvadm.data());
      }
    }
    TriggerSFmap["wTrgUp"] = std::isnan(t_trg_Up) ? 0 : t_trg_Up;
    TriggerSFmap["wTrgDown"] = std::isnan(t_trg_Down) ? 0 : t_trg_Down;
  }
  //
  return TriggerSFmap;
}

std::map<std::string, double> weight::BTaggingSF(std::vector<const pat::Jet*> selectedJets, Int_t theYear, std::string sysType, TH2D* btagEfficiency, const BTagCalibration *calib, BTagCalibrationReader *reader) {

  std::map<std::string, double> BTaggingSFmap;
  //
  double CSVcut = 0.;
  //
  if(theYear == 2016) CSVcut = 0.6321;
  else if(theYear == 2017) CSVcut = 0.4941;
  else if(theYear == 2018) CSVcut = 0.4184;
  // 
  double P_mc = 1.;
  double P_data = 1.;
  double P_dataUp = 1.;
  double P_dataDown = 1.;
  for(std::vector<const pat::Jet*>::iterator iJet = selectedJets.begin(); iJet != selectedJets.end(); ++iJet) {
    //
    double SF = reader->eval_auto_bounds("central", BTagEntry::FLAV_B, std::abs((*iJet)->eta()), (*iJet)->pt());
    double epsilon = btagEfficiency->GetBinContent((*iJet)->pt(), (*iJet)->eta());
    if(((*iJet)->bDiscriminator("pfDeepCSVJetTags:probb") + (*iJet)->bDiscriminator("pfDeepCSVJetTags:probbb"))>CSVcut && (*iJet)->pt()>20 && abs((*iJet)->eta())<2.4) {
      P_mc *= epsilon;
      P_data *= SF*epsilon;
    }
    else {
      P_mc *= (1 - epsilon);
      P_data *= (1 - SF*epsilon);
    }
    //
    if(sysType == "Nominal") {
      //
      double SFUp = reader->eval_auto_bounds("up", BTagEntry::FLAV_B, std::abs((*iJet)->eta()), (*iJet)->pt());
      double SFDown = reader->eval_auto_bounds("down", BTagEntry::FLAV_B, std::abs((*iJet)->eta()), (*iJet)->pt());
      if(((*iJet)->bDiscriminator("pfDeepCSVJetTags:probb") + (*iJet)->bDiscriminator("pfDeepCSVJetTags:probbb"))>CSVcut && (*iJet)->pt()>20 && abs((*iJet)->eta())<2.4) {
        P_dataUp *= SFUp*epsilon;
        P_dataDown *= SFDown*epsilon;
      }
      else {
        P_dataUp *= (1 - SFUp*epsilon);
        P_dataDown *= (1 - SFDown*epsilon);
      }
    }
  }
  BTaggingSFmap["wBtag"] = P_data/P_mc;
  if(sysType == "Nominal") {
    BTaggingSFmap["wBtagUp"] = P_dataUp/P_mc;
    BTaggingSFmap["wBtagDown"] = P_dataDown/P_mc;
  }
  return BTaggingSFmap;
}

//////

double weight::SignalReweighting(Int_t theYear, int id) {

  double wSignal = 1.;
  //
  if(theYear == 2016) {
    if(id == DataMCType::H_tautau_ggF) wSignal = 0.2455;
    else if(id == DataMCType::H_tautau_VBF) wSignal = 0.2727;
    else if(id == DataMCType::ZH_tautau) wSignal = 0.2546;
    else if(id == DataMCType::WminusH_tautau) wSignal = 0.2596;
    else if(id == DataMCType::WplusH_tautau) wSignal = 0.2425;
  }
  if(theYear == 2016) {
    if(id == DataMCType::H_tautau_ggF) wSignal = 0.2447;
    else if(id == DataMCType::H_tautau_VBF) wSignal = 0.2697;
    else if(id == DataMCType::ZH_tautau) wSignal = 0.2514;
    else if(id == DataMCType::WminusH_tautau) wSignal = 0.2567;
    else if(id == DataMCType::WplusH_tautau) wSignal = 0.2394;
  }
  if(theYear == 2016) {
    if(id == DataMCType::H_tautau_ggF) wSignal = 0.2446;
    else if(id == DataMCType::H_tautau_VBF) wSignal = 0.2695;
    else if(id == DataMCType::ZH_tautau) wSignal = 0.2513;
    else if(id == DataMCType::WminusH_tautau) wSignal = 0.2563;
    else if(id == DataMCType::WplusH_tautau) wSignal = 0.2397;
  }
  return wSignal;
}

double weight::PileUpreweighting(float nPU, TH1D* PU_data, TH1D* PU_mc) {

  double wPU = 1.;
  PileUp *PUofficial = new PileUp();
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);
  wPU = double(PUofficial->get_PUweight(double(nPU)));
  delete PUofficial;
  //
  return wPU;
}

double Stitching(Int_t theYear, int id, int taugenmatch, int mugenmatch, double mvis, bool isZ, bool isW, int Npartons) {

  /*bool isTausample = true;
  if(taugenmatch == 1 || taugenmatch == 2 || mugenmatch == 1 || mugenmatch == 2) isTausample = false;
  // 
  double wStitch = 1.;
  //
  if(isZ && isTausample) {
    if(id == DataMCType::DY_ll_10to50) {
      if(!isTausample) {
        if(theYear == 2016) wStitch = 19.0080307;
        else if(theYear == 2017) wStitch = 19.5191962215717;
        else if(theYear == 2018) wStitch = 28.2040833505999;
      }
      else wStitch = -9999;
    }
    else {
      if(isTausample) {
        if(Npartons == 0 || Npartons >= 5) {
          if(mvis < 150) {
            if(theYear == 2016) wStitch = 1.49005321266736;
	    else if(theYear == 2017) wStitch = ;
            else if(theYear == 2018) wStitch = 28.2040833505999;*/
  return 1.;
}
