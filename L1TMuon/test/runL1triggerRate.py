import FWCore.ParameterSet.Config as cms
process = cms.Process('L1')
# import of standard configuration
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load("Configuration.StandardSequences.L1Extra_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

from CalibMuon.CSCCalibration.CSCIndexer_cfi import CSCIndexerESProducer
process.CSCIndexerESProducer= CSCIndexerESProducer
from CalibMuon.CSCCalibration.CSCChannelMapper_cfi import CSCChannelMapperESProducer
process.CSCChannelMapperESProducer= CSCChannelMapperESProducer

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
fileOutputName = "out_L1.root"
histofileName  = "histo_L1.root"
# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring("file:out_digi_106.root")
    fileNames = cms.untracked.vstring("/store/user/calabria/NuGunGS/crab_SingleNu_TP2023HGCALDR_SLHC28/160708_135644/0000/out_digi_106.root")
)

FileList = 'SingleNu_TP2023HGCALDR_SLHC28.txt'
ff = open(FileList, "r")
files = ff.read().split('\n')
ff.close()
files = filter(lambda x: x.endswith('.root'),  files)
process.source.fileNames = files

# Output definition
outCommands = cms.untracked.vstring('drop *')
outCommands.append('keep *_genParticles_*_*')
outCommands.append('keep *_simCsctfDigis_*_*')
outCommands.append('keep *_simCsctfTrackDigis_*_*')
outCommands.append('keep *_simCscTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_g4SimHits_*_*')
outCommands.append('keep *_simMuonGEMCSCPadDigis_*_*')

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = outCommands,
    fileName = cms.untracked.string(fileOutputName),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('L1')
    )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(histofileName)
    )

# Other statements
process.L1TTriggerRate = cms.EDAnalyzer('L1TTriggerRate',
    minPt = cms.double(0.0),
    maxPt = cms.double(1000.0),
    minEta = cms.double(1.5),
    maxEta = cms.double(3.5),
    SRLUT = cms.PSet(
        Binary = cms.untracked.bool(False),
        ReadLUTs = cms.untracked.bool(False),
        LUTPath = cms.untracked.string('./'),
        UseMiniLUTs = cms.untracked.bool(True)
    ),
    debugTF = cms.bool(False)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.L1Extra_step = cms.Path(process.L1Extra)
process.pL1TTriggerRate = cms.Path(process.L1TTriggerRate)

#process.schedule = cms.Schedule(process.L1simulation_step,process.L1Extra_step,process.endjob_step,process.FEVTDEBUGHLToutput_step,process.pL1TTriggerRate)
process.schedule = cms.Schedule(process.endjob_step,process.FEVTDEBUGHLToutput_step,process.pL1TTriggerRate)

#from RecoParticleFlow.PandoraTranslator.customizeHGCalPandora_cff import cust_2023HGCalPandoraMuon 
#process = cust_2023HGCalPandoraMuon(process)
