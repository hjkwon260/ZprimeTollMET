import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('isMC',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "isMC")
options.parseArguments()

print(options.isMC)

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/u/user/hjkwon/testKW/CMSSW_9_4_13/src/ZprimeTollMET/NtupleMaker/python/94XMC/8C092947-1022-E911-8002-6CC2173DC2E0.root'
        # 'file:/u/user/hjkwon/SE_UserHome/CRAB_UserFiles/crab_MINIAOD_BMII_2p5_mu_PU/191205_151624/0000/ZpBSM_pythia8_MiniAOD_10.root'
        # 'file:/u/user/hjkwon/testKW/CMSSW_10_2_0/src/ZprimeTollMET/NtupleMaker/test/F864DF5E-948B-E811-9770-A4BF01025B08.root'
        '/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/D6B5847A-10C5-E811-9B22-A4BF0108B062.root'
    )
)


process.ntuple = cms.EDAnalyzer('NtupleMaker',
	isMC = cms.bool(bool(options.isMC)),
	triggerresults = cms.untracked.InputTag("TriggerResults", "", "HLT"), # only hlt
	triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"), # only PAT
	triggerobjects = cms.untracked.InputTag("slimmedPatTrigger"),
	# triggerobjects = cms.untracked.InputTag("selectedPatTrigger"), # for signal MC
	# triggerresultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"),
  genparticles = cms.untracked.InputTag("prunedGenParticles"),
	genevtinfo = cms.untracked.InputTag("generator"),
  muons = cms.untracked.InputTag("slimmedMuons"),
	electrons = cms.untracked.InputTag("slimmedElectrons"),
	pvs = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
	# jets = cms.untracked.InputTag("slimmedJets"),
	jets = cms.untracked.InputTag("selectedUpdatedPatJetsNewDFTraining"),
	# DeepFlavour = cms.untracked.InputTag("selectedUpdatedPatJetsNewDFTraining"),
	MET = cms.untracked.InputTag("slimmedMETs")
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
    process.GlobalTag.globaltag = cms.string("102X_mcRun2_asymptotic_v7")
else:
    process.GlobalTag.globaltag = cms.string("102X_dataRun2_v12")

#-- Deep flavor rerun --#
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
   btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg'
      ],
   postfix='NewDFTraining'
)

#-- electron --##
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy')  
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

process.p = cms.Path(process.ntuple*process.egammaPostRecoSeq)
process.p.associate(process.patAlgosToolsTask)

#-- Local test --#
if options.isMC==1:
    # process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/D6B5847A-10C5-E811-9B22-A4BF0108B062.root')
    process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEFC23FC-37C7-E811-97DA-0CC47AA53D86.root')
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
