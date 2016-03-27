runevents = 1000
runevents = -1;

import FWCore.ParameterSet.Config as cms
process = cms.Process('HitRateAnalysis')
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
datadir = '/pnfs/user/jlee/data/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023HGCALDR-HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/GEN-SIM-DIGI-RAW/'
#datadir = '/pnfs/user/jlee/data/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023HGCALDR-HGCALForMUO_PU200BX25_newsplit_PH2_1K_FB_V6-v1/GEN-SIM-DIGI-RAW/'
for f in os.listdir(datadir):
    if '50983EA5-1FE2-E411-91B1-008CFA010D18.root' in f:
        continue
    if '.root' in f:
        process.source.fileNames.append("file:"+datadir+f)

#print "fileNames: ", process.source.fileNames
process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
    histofileName
))

process.HitRateAnalysis = cms.EDAnalyzer('HitRateAnalysis',
    gemDigiInput      = cms.InputTag("simMuonGEMDigis"),
    gemPadDigiInput   = cms.InputTag("simMuonGEMCSCPadDigis"),
    gemCoPadDigiInput = cms.InputTag("simCscTriggerPrimitiveDigis"),
)

process.p = cms.Path(process.HitRateAnalysis)


## for debugging
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
