# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt100_cfi --conditions auto:upgrade2019 -n 100 --eventcontent FEVTDEBUGHLT -s L1 --datatier L1 --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2019WithGem --geometry Extended2019 --magField 38T_PostLS1 --filein file:out_digi.root --fileout file:out_L1_cust_2019WithGem.root

useUpdatedTF = True
withGEM = True
clctNplanesHitPattern3 = False
buildLCTfromALCTandGEMinOverlap = False

import FWCore.ParameterSet.Config as cms

process = cms.Process('L1')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

fileOutputName = "out_L1_cust_2019"
if withGEM:
    fileOutputName = fileOutputName+"WithGem"

if useUpdatedTF:
    process.load('L1TriggerDPGUpgrade.L1TMuon.L1TMuonTriggerPrimitiveProducer_cfi')
    process.load('L1TriggerDPGUpgrade.L1TMuon.L1CSCTFTrackConverter_cfi')
    process.load('L1TriggerDPGUpgrade.L1TMuon.L1DTTFTrackConverter_cfi')
    process.load('L1TriggerDPGUpgrade.L1TMuon.L1RPCTFTrackConverter_cfi')
    process.load('L1TriggerDPGUpgrade.L1TMuon.L1TMuonSimpleDeltaEtaHitMatcher_cfi')
    fileOutputName = fileOutputName+"NewTF"

histofileName= "histo_"+fileOutputName+".root"
fileOutputName = "file:"+fileOutputName+".root"
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

## input
directory = '/pnfs/user/trackfinder/tf/'
inputDir = [directory]
theInputFiles = []
import os
for d in range(len(inputDir)):
  my_dir = inputDir[d]
  if not os.path.isdir(my_dir):
    print "ERROR: This is not a valid directory: ", my_dir
    if d==len(inputDir)-1:
      print "ERROR: No input files were selected"
      exit()
    continue
  print "Proceed to next directory"
  ls = os.listdir(my_dir)
  ## this works only if you pass the location on pnfs - FIXME for files staring with store/user/...
#  theInputFiles.extend([my_dir[:16] + x for x in ls if x.endswith('root')])
  theInputFiles.extend(['file:' + my_dir[:50] + x for x in ls if x.endswith('.root')])

theInputFiles = theInputFiles[:60]
##inputFiles = ['file:out_SingleMuPt10Fwd_GEM2019_8PartIncRad_DIGI_L1.root']
print "InputFiles: ", theInputFiles

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring("file:" + directory + "tf/out_digi_test_73_1_M7u.root")
    fileNames = cms.untracked.vstring(*theInputFiles)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleMuPt100_cfi nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition
outCommands = cms.untracked.vstring('keep *')

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = outCommands,
    fileName = cms.untracked.string(fileOutputName),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('L1')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.L1simulation_step = cms.Path(process.SimL1Emulator)

if not useUpdatedTF:
    process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

if useUpdatedTF:
    process.L1NewTF_step = cms.Path( process.L1TMuonTriggerPrimitives +
                                     process.L1CSCTFTrackConverter    +
                                     process.L1DTTFTrackConverter     +
                                     process.L1RPCTFTrackConverters   +
                                     process.L1TMuonSimpleDeltaEtaHitMatcher )
    
    process.L1TMuonText = cms.EDProducer(
        'L1TMuonTextDumper',
        doGen = cms.untracked.bool(True),
        genSrc = cms.untracked.InputTag("genParticles"),
        primitiveSrcs = cms.VInputTag(
            cms.InputTag('L1TMuonTriggerPrimitives','CSC'),
            cms.InputTag('L1TMuonTriggerPrimitives','DT'),
            cms.InputTag('L1TMuonTriggerPrimitives','RPC'),
            cms.InputTag('L1TMuonTriggerPrimitives','GEM')
            ),
        converterSrcs = cms.VInputTag(
            cms.InputTag('L1CSCTFTrackConverter'),
            cms.InputTag('L1DTTFTrackConverter'),
            cms.InputTag('L1RPCbTFTrackConverter'),
            cms.InputTag('L1RPCfTFTrackConverter'),
            cms.InputTag('L1TMuonSimpleDeltaEtaHitMatcher')
            ),
        lutParam = process.simCsctfTrackDigis.SectorProcessor.PTLUT
                #cms.PSet(
                #    isBeamStartConf = cms.untracked.bool(True),
                #    ReadPtLUT = cms.bool(False),
                #    PtMethod = csctfTrackDigisUngangedME1a.SectorProcessor.PTLUT.PtMethod #cms.untracked.uint32(32)
                #    )
        )
    
    process.content = cms.EDAnalyzer("EventContentAnalyzer")
    
    
    process.sptf = cms.EDProducer('sptf',
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
    process.DiagMaker = cms.EDAnalyzer(
        'DiagMaker',
        testVar = cms.untracked.bool(False)
        )
    process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(
                                           histofileName
                                           ))
    
    #process.L1TMuonSequence = cms.Sequence(process.L1TMuonVerilogBasedMatcher * process.DiagMaker)
    #process.L1TMuonSequence = cms.Sequence(process.L1TMuonVerilogBasedMatcher)
    process.L1TMuonSequence = cms.Sequence(process.L1TMuonText)
    #process.L1TMuonSequence = cms.Sequence(process.L1TMuonText * process.sptf * process.DiagMaker)
    process.L1TMuonPath = cms.Path(process.L1TMuonSequence)
    process.schedule = cms.Schedule(process.L1simulation_step,process.L1NewTF_step,process.L1TMuonPath,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
if withGEM:
    from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019WithGem 
    process = cust_2019WithGem(process)

    if clctNplanesHitPattern3:
        process.simCscTriggerPrimitiveDigis.clctSLHC.clctNplanesHitPattern = cms.uint32(3)
    if buildLCTfromALCTandGEMinOverlap:
        process.simCscTriggerPrimitiveDigis.tmbSLHC.buildLCTfromALCTandGEMinOverlap = cms.untracked.bool(True)

else :
    from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019
    process = cust_2019(process)

print process.simCsctfTrackDigis.SectorProcessor.PTLUT
if useUpdatedTF:
    process.L1TMuonText.lutParam = process.simCsctfTrackDigis.SectorProcessor.PTLUT

# End of customisation functions
