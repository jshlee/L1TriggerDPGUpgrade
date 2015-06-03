runevents = 1000
runevents = -1;


import FWCore.ParameterSet.Config as cms
process = cms.Process('L1analysis')

process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuGMTScalesConfig_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#theInputFiles = ['file:/u/user/jlee/scratch/L1data0417/out_L1_cust_2019NewTF.root']
#out_L1muon2023.root
#out_L1muon2023GE11.root
#out_L1muon2023GE21.root

#fileOutputName = "SingleMu_SLHC12_PU0"
#theInputFiles = ['file:SingleMu_SLHC12_PU0.root']

fileOutputName = "test"
theInputFiles = ['file:/eos/uscms/store/user/lpcgem/dildick/DisplacedMuGEMCSCILT/dildick/DarkSUSY_mH_125_mGammaD_0400_ctau_00_14TeV_madgraph452_bridge224_LHE_pythia6_GEN_SIM/DarkSUSY_mH_125_mGammaD_0400_ctau_00_14TeV_madgraph452_bridge224_LHE_pythia6_DIGI_L1/75d0e84afc642f0aede2fcf263b4fa2e/out_L1_1_1_O8R.root']

#fileOutputName = "out_L1muon2023GE21"
#theInputFiles = ['file:out_L1muon2023GE21.root']

#fileOutputName = "out_L1muon2023GE11"
#theInputFiles = ['file:out_L1muon2023GE11.root']

#fileOutputName = "out_L1_cust_2019NewTF"
#theInputFiles = ['file:/u/user/jlee/scratch/L1data0423/out_L1_cust_2019NewTF.root']
#fileOutputName = "out_L1_cust_2019PtMethod32NewTF"
#theInputFiles = ['file:L1data0423/out_L1_cust_2019PtMethod32NewTF.root']
#fileOutputName = "out_L1_cust_2019UnFlipHsInLutWithGemNewTF"
#theInputFiles = ['file:L1data0423/out_L1_cust_2019UnFlipHsInLutWithGemNewTF.root']
#fileOutputName = "out_L1_cust_2019WithGemNewTF"
#theInputFiles = ['file:L1data0423/out_L1_cust_2019WithGemNewTF.root']

histofileName= "histo_"+fileOutputName+".root"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(runevents)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring("file:~/cmsrun/out_digi.root")
    fileNames = cms.untracked.vstring(*theInputFiles)
)
print "fileNames: ", process.source.fileNames
process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleMuPt100_cfi nevts:100'),
    name = cms.untracked.string('Applications')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
    histofileName
))

process.L1TAnalyser = cms.EDAnalyzer('L1TAnalyser',
    minPt = cms.double(2.0),
    maxPt = cms.double(100.0),
    minEta = cms.double(1.6),
    maxEta = cms.double(2.4),
    SRLUT = cms.PSet(
        Binary = cms.untracked.bool(False),
        ReadLUTs = cms.untracked.bool(False),
        LUTPath = cms.untracked.string('./'),
        UseMiniLUTs = cms.untracked.bool(True)
    ),
    debugTF = cms.bool(False)
)

process.pL1TAnalyser = cms.Path(process.L1TAnalyser)

process.schedule = cms.Schedule(process.pL1TAnalyser)
