import FWCore.ParameterSet.Config as cms

process = cms.Process('TEXTDUMP')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1))
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("Configuration.StandardSequences.Geometry_cff")



process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.L1TMuonText = cms.EDProducer(
    'L1TMuonTextDumper',
    doGen = cms.untracked.bool(True),
    genSrc = cms.untracked.InputTag("genParticles"),
    primitiveSrcs = cms.VInputTag(
    cms.InputTag('L1TMuonTriggerPrimitives','CSC'),
    cms.InputTag('L1TMuonTriggerPrimitives','DT'),
    cms.InputTag('L1TMuonTriggerPrimitives','RPC')
    ),
    converterSrcs = cms.VInputTag(
    cms.InputTag('L1CSCTFTrackConverter'),
    cms.InputTag('L1DTTFTrackConverter'),
    cms.InputTag('L1RPCbTFTrackConverter'),
    cms.InputTag('L1RPCfTFTrackConverter'),
    cms.InputTag('L1TMuonSimpleDeltaEtaHitMatcher')
    ),
    lutParam = cms.PSet(
    isBeamStartConf = cms.untracked.bool(True),
    ReadPtLUT = cms.bool(False),
    PtMethod = cms.untracked.uint32(32)
    )
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")


process.sptf = cms.EDProducer(
    'sptf',
    primitiveSrcs = cms.VInputTag(
    cms.InputTag('L1TMuonTriggerPrimitives','CSC'),
    cms.InputTag('L1TMuonTriggerPrimitives','DT'),
    cms.InputTag('L1TMuonTriggerPrimitives','RPC')
    ),
    converterSrcs = cms.VInputTag(
    cms.InputTag('L1CSCTFTrackConverter'),
    cms.InputTag('L1DTTFTrackConverter'),
    cms.InputTag('L1RPCbTFTrackConverter'),
    cms.InputTag('L1RPCfTFTrackConverter'),
    cms.InputTag('L1TMuonSimpleDeltaEtaHitMatcher')
    ),
    lutParam = cms.PSet(
    isBeamStartConf = cms.untracked.bool(True),
    ReadPtLUT = cms.bool(False),
    PtMethod = cms.untracked.uint32(32)
    )
)


##process.L1TMuonVerilogBasedMatcher = cms.EDProducer(
##    'L1TMuonVerilogBasedMatcher',
##   
##    primitiveSrcs = cms.VInputTag(
##    cms.InputTag('L1TMuonTriggerPrimitives','CSC')
#    cms.InputTag('L1TMuonTriggerPrimitives','CSC'),
#    cms.InputTag('L1TMuonTriggerPrimitives','DT'),
#    cms.InputTag('L1TMuonTriggerPrimitives','RPC')
##    ),
##    converterSrcs = cms.VInputTag(
##    cms.InputTag('L1CSCTFTrackConverter'),
##    cms.InputTag('L1DTTFTrackConverter'),
##    cms.InputTag('L1RPCbTFTrackConverter'),
##    cms.InputTag('L1RPCfTFTrackConverter'),
##    cms.InputTag('L1TMuonSimpleDeltaEtaHitMatcher')
##    ),
##    lutParam = cms.PSet(
##    isBeamStartConf = cms.untracked.bool(True),
##    ReadPtLUT = cms.bool(False),
##    PtMethod = cms.untracked.uint32(32)
    
##    )
##)


process.DiagMaker = cms.EDAnalyzer(
    'DiagMaker',
    testVar = cms.untracked.bool(False)
)

infile = 'file:L1TMuon2.root'
fileOutName = "test.root"

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    
    #eventsToProcess = cms.untracked.VEventRange('1:590','1:602','1:640','1:906','1:1029','1:1318','1:1622','1:1765','1:1817','1:2205','1:2206','1:2494','1:2541')
    #eventsToProcess = cms.untracked.VEventRange('1:1925-1:1927','1:2769-1:2771','1:3669-1:3671'),
    #eventsToSkip= cms.untracked.VEventRange('1:896')
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
	fileOutName
))

#process.L1TMuonSequence = cms.Sequence(process.L1TMuonVerilogBasedMatcher * process.DiagMaker)
#process.L1TMuonSequence = cms.Sequence(process.L1TMuonVerilogBasedMatcher)
process.L1TMuonSequence = cms.Sequence(process.L1TMuonText)
#process.L1TMuonSequence = cms.Sequence(process.L1TMuonText * process.sptf * process.DiagMaker)

process.L1TMuonPath = cms.Path(process.L1TMuonSequence)

process.schedule = cms.Schedule(process.L1TMuonPath)

