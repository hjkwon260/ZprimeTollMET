from CRABClient.UserUtilities import config
import datetime
now = datetime.datetime.now()
date = now.strftime('%Y%m%d')
config = config()

# version = 'v070420'

config.General.requestName = ''
config.General.workArea = 'ZprimellMET_v0701'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = './NtupleMaker_cfg_16.py'
config.JobType.pyCfgParams = ['isMC=1']
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
config.Data.unitsPerJob = 1
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/hkwon/ZprimellMET/%s' % (date)
config.Data.publication = False

# config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T3_KR_KNU'


### -- 'MultiCRAB' part -- ###

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # -- DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM -- #
    config.General.requestName = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' 
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    # # -- DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX -- #
    config.General.requestName = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX' 
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)
    
    # # -- TT_TuneCUETP8M2T4_13TeV-powheg -- #
    config.General.requestName = 'TT_TuneCUETP8M2T4_13TeV-powheg' 
    config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    # config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    # # -- TTTo2L2Nu_TuneCUETP8M2_ttHtranche3 -- #
    # config.General.requestName = 'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3' 
    # config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
    # crabCommand('submit', config = config)

    # # -- TTTo2L2Nu_TuneCP5_PSweights -- #
    # config.General.requestName = 'TTTo2L2Nu_TuneCP5_PSweights' 
    # config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    # crabCommand('submit', config = config)

    # -- ST_tW_antitop -- #
    config.General.requestName = 'ST_tW_antitop' 
    config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # -- ST_tW_top -- #
    config.General.requestName = 'ST_tW_top' 
    config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # -- WW -- #
    config.General.requestName = 'WW_TuneCUETP8M1_13TeV-pythia8' 
    config.Data.inputDataset = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    # -- ZZ -- #
    config.General.requestName = 'ZZ_TuneCUETP8M1_13TeV-pythia8' 
    config.Data.inputDataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    # -- WZ -- #
    config.General.requestName = 'WZ_TuneCUETP8M1_13TeV-pythia8' 
    config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    # -- Charginos inclusive -- #
    # config.General.requestName = 'SMS-TChipmWW_WWTo2LNu' 
    # config.Data.inputDataset = '/SMS-TChipmWW_WWTo2LNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
    # crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)