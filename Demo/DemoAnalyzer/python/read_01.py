import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#               'file:step3.root'
#               'file:step_pi100GeV.root'
               'file:step02_pi50GeV.root'
#               'file:step6.root'
    )
)

process.demo = cms.EDAnalyzer('DemoAnalyzer')

process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('histo_example_02.root') )
    fileName = cms.string('example_pi50Gev_new.root') )

process.p = cms.Path(process.demo)
