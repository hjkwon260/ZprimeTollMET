import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('isMC',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "isMC")
options.parseArguments()

# print(options.isMC)

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/u/user/hjkwon/testKW/CMSSW_9_4_13/src/ZprimeTollMET/NtupleMaker/python/94XMC/8C092947-1022-E911-8002-6CC2173DC2E0.root'
        # 'file:/u/user/hjkwon/SE_UserHome/CRAB_UserFiles/crab_MINIAOD_BMII_2p5_mu_PU/191205_151624/0000/ZpBSM_pythia8_MiniAOD_10.root'
        # 'file:/u/user/hjkwon/testKW/CMSSW_10_2_0/src/ZprimeTollMET/NtupleMaker/test/F864DF5E-948B-E811-9770-A4BF01025B08.root'
        # '/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/D6B5847A-10C5-E811-9B22-A4BF0108B062.root'
      # '/store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/70000/2EE3B436-E445-A440-A629-89CE2962EC9B.root',
      '/store/mc/RunIISummer20UL17MiniAODv2/ZprimeTo2ChiTo2L_mZp-3300_mCH-1345_TuneCP2_13TeV-madgraph-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/2560000/61A65A8D-B9BB-E74E-A640-355954D2E77C.root'
    )
)


process.ntuple = cms.EDAnalyzer('NtupleMaker',
  isMC = cms.bool(bool(options.isMC)),
  pileupsummary = cms.untracked.InputTag("slimmedAddPileupInfo"),
  triggerresults = cms.untracked.InputTag("TriggerResults", "", "HLT"), # only hlt
  triggerresultsRECO = cms.untracked.InputTag("TriggerResults", "", "RECO"), # only PAT
  triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"), # only PAT
  # triggerobjects = cms.untracked.InputTag("slimmedPatTrigger"), # comment out, no susy
  # triggerobjects = cms.untracked.InputTag("selectedPatTrigger"), # for signal MC
  # triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"),
  prefiringweight = cms.untracked.InputTag("prefiringweight:nonPrefiringProb"),
  prefiringweightUp = cms.untracked.InputTag("prefiringweight:nonPrefiringProbUp"),
  prefiringweightDown = cms.untracked.InputTag("prefiringweight:nonPrefiringProbDown"),
  genparticles = cms.untracked.InputTag("prunedGenParticles"),
  genevtinfo = cms.untracked.InputTag("generator"),
  lheevtproduct = cms.untracked.InputTag("externalLHEProducer"),
  lheruninfo = cms.untracked.InputTag("externalLHEProducer"),
  muons = cms.untracked.InputTag("slimmedMuons"),
  electrons = cms.untracked.InputTag("slimmedElectrons"),
  pvs = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  # jets = cms.untracked.InputTag("slimmedJets"),
  jets = cms.untracked.InputTag("selectedUpdatedPatJetsNewDFTraining"),
  # DeepFlavour = cms.untracked.InputTag("selectedUpdatedPatJetsNewDFTraining"),
  MET = cms.untracked.InputTag("slimmedMETs"),
  puppiMET = cms.untracked.InputTag("slimmedMETsPuppi"),
  deepMET = cms.untracked.InputTag("deepMETProducer")
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("ntuple.root"),
  closeFileFast = cms.untracked.bool(False),
  )

#-- needed to rerun PAT jet --#
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi") # Exception Message: No "TransientTrackRecord" record found in the EventSetup.n
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if options.isMC==1:
    process.GlobalTag.globaltag = cms.string("106X_mcRun2_asymptotic_preVFP_v11")
else:
    process.GlobalTag.globaltag = cms.string("102X_dataRun2_v12")

# -- l1 prefireing (16, 17 ?)--#
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2016BtoH"), #Use 2017BtoF for 2017
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

# #-- Deep flavor rerun --#
# from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

# updateJetCollection(
#    process,
#    jetSource = cms.InputTag('slimmedJets'),
#    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
#    svSource = cms.InputTag('slimmedSecondaryVertices'),
#    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#    btagDiscriminators = [
#       'pfDeepFlavourJetTags:probb',
#       'pfDeepFlavourJetTags:probbb',
#       'pfDeepFlavourJetTags:problepb',
#       'pfDeepFlavourJetTags:probc',
#       'pfDeepFlavourJetTags:probuds',
#       'pfDeepFlavourJetTags:probg'
#       ],
#    postfix='NewDFTraining'
# )

# #-- electron --##
# from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# setupEgammaPostRecoSeq(process,
#                        runEnergyCorrections=False, #corrections by default are fine so no need to re-run
#                        era='2016-Legacy')  
# #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

# # -- Deep MET --#
# from RecoMET.METPUSubtraction.deepMETProducer_cfi import deepMETProducer
# process.deepMETProducer = deepMETProducer.clone()
# process.sequence = cms.Sequence(process.deepMETProducer)

# process.p = cms.Path(process.prefiringweight*process.egammaPostRecoSeq*process.sequence*process.ntuple)
process.p = cms.Path(process.prefiringweight*process.ntuple)
# process.p = cms.Path(process.prefiringweight*process.egammaPostRecoSeq*process.sequence)
# process.p.associate(process.patAlgosToolsTask)

#-- Local test --#
if options.isMC==1:
    # process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/D6B5847A-10C5-E811-9B22-A4BF0108B062.root')
    # process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/70000/2EE3B436-E445-A440-A629-89CE2962EC9B.root')
    process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer20UL16MiniAODAPVv2/ZprimeTo2ChiTo2L_mZp-3300_mCH-1345_TuneCP2_13TeV-madgraph-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2560000/0D9C6432-8901-F047-AA82-89872F93DEAF.root')
    # process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer20UL17MiniAODv2/ZprimeTo2ChiTo2L_mZp-3300_mCH-1345_TuneCP2_13TeV-madgraph-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/2560000/61A65A8D-B9BB-E74E-A640-355954D2E77C.root')
    # process.source.fileNames = cms.untracked.vstring('file:FCE4D694-BDD2-E911-BB3F-0CC47A78A340.root')
    # process.source.fileNames = cms.untracked.vstring(
    #   'file:Zprime_reMiniaod94X_2p5_PU_mu_1.root',
    #   'file:Zprime_reMiniaod94X_2p5_PU_mu_2.root',
    #   'file:Zprime_reMiniaod94X_2p5_PU_mu_3.root',
    #   'file:Zprime_reMiniaod94X_2p5_PU_mu_4.root',
      
    #   )
    print(process.source.fileNames)
else:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver1-v1/80000/F2494388-1D8C-E811-BB8E-0242AC1C0502.root')

# in case genevtinfo absent(?)
# process.options   = cms.untracked.PSet( 
#   wantSummary = cms.untracked.bool(True),
#   SkipEvent = cms.untracked.vstring('ProductNotFound') 
# )

#-- Write all instances --#
process.writeDataset = cms.OutputModule("PoolOutputModule",
  # splitLevel = cms.untracked.int32(0),
  # eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
  # outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
  fileName = cms.untracked.string('output_dataset.root'), ## ADAPT IT ##
  # dataset = cms.untracked.PSet(
  #     filterName = cms.untracked.string(''),
  #     dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
  # )
)

# process.pd = cms.EndPath(process.writeDataset)
