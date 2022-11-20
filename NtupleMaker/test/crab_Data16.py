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
config.JobType.pyCfgParams = ['isMC=0']
config.JobType.allowUndistributedCMSSW = True #CMSSW_10_2_10 on slc7_amd64_gcc700 is not among supported releases?

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/hkwon/ZprimellMET/%s' % (date)
config.Data.publication = False

#config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T3_KR_KNU'


### -- 'MultiCRAB' part -- ###

# GoldenJSON = './Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
GoldenJSON = './Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
# MuonPhysJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
StartRun = 271036
EndRun = 284044

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

### -- SingleMuon -- ###

    # -- Run2016B, ver1 -- #
#     config.General.requestName = 'SingleMuon_Run2016Bver1'
   # config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver1-v1/MINIAOD'
   # config.Data.lumiMask = GoldenJSON
   # config.Data.runRange = '%d-%d' % (StartRun, EndRun)
   # crabCommand('submit', config = config)
   # not in Golden json?

    # -- Run2016B, ver2-- #
    config.General.requestName = 'SingleMuon_Run2016Bver2'
    config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016C -- #
    config.General.requestName = 'SingleMuon_Run2016C'
    config.Data.inputDataset = '/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016D -- #
    config.General.requestName = 'SingleMuon_Run2016D'
    config.Data.inputDataset = '/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016E -- #
    config.General.requestName = 'SingleMuon_Run2016E'
    config.Data.inputDataset = '/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016F -- #
    config.General.requestName = 'SingleMuon_Run2016F'
    config.Data.inputDataset = '/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016G -- #
    config.General.requestName = 'SingleMuon_Run2016G'
    config.Data.inputDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
 
    # -- Run2016H -- #
    config.General.requestName = 'SingleMuon_Run2016H'
    config.Data.inputDataset = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)
#
    # ### -- SingleElectron -- ###

    # -- Run2016B, ver1 -- #
    # config.General.requestName = 'SingleElectron_Run2016Bver1'
    # config.Data.inputDataset = '/SingleElectron/Run2016B-17Jul2018_ver1-v1/MINIAOD'
    # config.Data.lumiMask = GoldenJSON
    # config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    # crabCommand('submit', config = config)
    # not in Golden json?

    # -- Run2016B, ver2-- #
    config.General.requestName = 'SingleElectron_Run2016Bver2'
    config.Data.inputDataset = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016C -- #
    config.General.requestName = 'SingleElectron_Run2016C'
    config.Data.inputDataset = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016D -- #
    config.General.requestName = 'SingleElectron_Run2016D'
    config.Data.inputDataset = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016E -- #
    config.General.requestName = 'SingleElectron_Run2016E'
    config.Data.inputDataset = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016F -- #
    config.General.requestName = 'SingleElectron_Run2016F'
    config.Data.inputDataset = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016G -- #
    config.General.requestName = 'SingleElectron_Run2016G'
    config.Data.inputDataset = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H -- #
    config.General.requestName = 'SingleElectron_Run2016H'
    config.Data.inputDataset = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)


    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
