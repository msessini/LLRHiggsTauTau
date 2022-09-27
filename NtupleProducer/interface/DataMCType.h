#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1,
	     H_tautau_ggF=11,
	     H_tautau_VBF=12,
             ZH_tautau=13,
             WplusH_tautau=141,
             WminusH_tautau=142,
	     W_lnu=20,
	     W_taunu=200,
	     W_1qlnu=201,
	     W_2qlnu=202,
	     W_3qlnu=203,
	     W_4qlnu=204,
	     DY_ll=30, 
 	     DY_tautau=300,
	     DY_1qll=301,
	     DY_2qll=302,
	     DY_3qll=303,
	     DY_4qll=304,
	     DY_ll_10to50=31,
	     DY_etau_embedded=34,
	     DY_mutau_embedded=35,
	     DY_tautau_embedded=36,
	     WZ_1l3nu=47,
	     WZ_1l1nu2q=48,
	     WW_1l1nu2q=49,
	     ZZ_4l=53,
	     ZZ_2l2nu=54,
	     ZZ_2l2q=55,
	     WW_2l2nu=56,
	     WZ_2l2q=57,
	     WZ_3l1nu=58,
	     WH0L1_H_tautau=59,
	     QCD=60, 
	     ttbar_dilep=701,
	     ttbar_hadr=702,
	     ttbar_semilep=703,
	     tw=71,
	     tbarw=72,
	     ST_tchannel_top=73,
	     ST_tchannel_antitop=74,
	     H_WW_2l2nu_ggF=75,
	     H_WW_2l2nu_VBF=76,
	     unknown=999,
	     EWKWplus2Jets_Wlnu=210,
	     EWKWminus2Jets_Wlnu=211,
	     EWKZ2Jets_Zll=212
  };

  DataMCType();
  ~DataMCType();

  unsigned int GetType(TString name);
  unsigned int SignalCode(unsigned int type,unsigned int JAK_ID1, unsigned int nprong1,unsigned int JAK_ID2, unsigned int nprong2);
  void DecodeSignal(unsigned int code,unsigned int &type,unsigned int &JAK_ID1, unsigned int &nprong1,unsigned int &JAK_ID2, unsigned int &nprong2);
  bool isSignalParticle(int pdg_id);
  void StoreType(TString t){type=t;}
  unsigned int GetType(){return GetType(type);}

 private:
  static TString type;

};
#endif
