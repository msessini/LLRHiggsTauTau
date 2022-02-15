# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "H_tautau_VBF_uncorr_filtered_2017_V3bis"
config.General.workArea = "crab3_H_tautau_VBF_uncorr_filtered_2017_V3bis"
config.General.transferLogs = True


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer.py' # to produce LLR ntuples or EnrichedMiniAOD according to the RunNtuplizer bool
config.JobType.inputFiles = (['/opt/sbg/cms/safe1/cms/msessini/IPHCProductionTools/CMSSW_10_2_23/src/LLRHiggsTauTau/NtupleProducer/data','/opt/sbg/cms/safe1/cms/msessini/IPHCProductionTools/CMSSW_10_2_23/src/LLRHiggsTauTau/NtupleProducer/test/JECUncertaintySources'])#/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/TauAnalysisTools/TauTriggerSFs/data/2017_tauTriggerEff_DeepTau2017v2p1.root','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/data/TauES_dm_2017ReReco.root','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/data/data_pileup_pudistributions_mc_2017.root','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/data/data_pileup_pudistributions_data_2017.root','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/test/JECUncertaintySources/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/test/JECUncertaintySources/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK4PFchs.txt','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/test/Fall17_17Nov2017C_V32_DATA_UncertaintySources_AK4PFchs.txt','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/test/Fall17_17Nov2017DE_V32_DATA_UncertaintySources_AK4PFchs.txt','/home-pbs/gbourgat/cmssw/CMSSW_10_2_16/src/LLRHiggsTauTau/NtupleProducer/test/Fall17_17Nov2017F_V32_DATA_UncertaintySources_AK4PFchs.txt'])
config.JobType.maxMemoryMB=2500
#config.JobType.sendExternalFolder = True #Needed until the PR including the Spring16 ele MVA ID is integrated in CMSSW/cms-data.

config.section_("Data")
config.Data.inputDataset = "/VBFHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = "EventAwareLumiBased"
config.Data.unitsPerJob = 40000
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/msessini/Signal_TauhTauh/2017/'
#config.Data.publication = True
config.Data.outputDatasetTag = "H_tautau_VBF_uncorr_filtered_2017_V3bis"
#config.Data.inputDBS = 'phys03'

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'
