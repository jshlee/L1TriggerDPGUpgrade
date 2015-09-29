runevents = 1000
runevents = -1;

import FWCore.ParameterSet.Config as cms
process = cms.Process('HitAnalysis')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

fileOutputName = "test"
histofileName= "histo_"+fileOutputName+".root"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(runevents)
)

# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
import os
datadir = '/pnfs/user/jlee/datagem/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023HGCALDR-HGCALForMUO_PU140BX25_newsplit_PH2_1K_FB_V6-v2/GEN-SIM-DIGI-RAW/'
for f in os.listdir(datadir):
    process.source.fileNames.append("file:"+datadir+f)

#print "fileNames: ", process.source.fileNames
process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
    histofileName
))

process.HitAnalysis = cms.EDAnalyzer('HitAnalysis',
    gemDigiInput      = cms.InputTag("simMuonGEMDigis"),
    gemPadDigiInput   = cms.InputTag("simMuonGEMCSCPadDigis"),
    gemCoPadDigiInput = cms.InputTag("simCscTriggerPrimitiveDigis"),
)

process.p = cms.Path(process.HitAnalysis)


## for debugging
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
