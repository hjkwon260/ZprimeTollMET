import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/u/user/hjkwon/testKW/CMSSW_9_4_13/src/ZprimeTollMET/NtupleMaker/python/94XMC/8C092947-1022-E911-8002-6CC2173DC2E0.root'
        'file:/u/user/hjkwon/SE_UserHome/CRAB_UserFiles/crab_MINIAOD_BMII_2p5_mu_PU/191205_151624/0000/ZpBSM_pythia8_MiniAOD_10.root'
        # 'file:/u/user/hjkwon/testKW/CMSSW_10_2_0/src/ZprimeTollMET/NtupleMaker/test/F864DF5E-948B-E811-9770-A4BF01025B08.root'
    )
)

process.ntuple = cms.EDAnalyzer('NtupleMaker',
	isMC = cms.bool(True),
	triggerresults = cms.untracked.InputTag("TriggerResults", "", "HLT"), # only hlt
	triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"), # only PAT
	# triggerobjects = cms.untracked.InputTag("slimmedPatTrigger"),
	triggerobjects = cms.untracked.InputTag("selectedPatTrigger"), # for signal MC
	# triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"),
	genparticles = cms.untracked.InputTag("prunedGenParticles"),
	muons = cms.untracked.InputTag("slimmedMuons"),
	pvs = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
	jets = cms.untracked.InputTag("slimmedJets"),
	MET = cms.untracked.InputTag("slimmedMETs")
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("ntuple.root"),
  closeFileFast = cms.untracked.bool(False),
  )

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("102X_mcRun2_asymptotic_v7")
# process.GlobalTag.globaltag = cms.string("102X_dataRun2_v12")



process.p = cms.Path(process.ntuple)
