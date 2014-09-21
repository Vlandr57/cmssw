# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_Tauola_13TeV_cfi --conditions auto:upgradePLS1 -n 10 --eventcontent FEVTDEBUGHLT -s GEN,SIM --datatier GEN-SIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry SimHCAL --magField 0T --fileout file:step1.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimHCAL_cff')
##process.load('Configuration.Geometry.GeometrySimECALHCAL_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
#process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('TTbar_Tauola_13TeV_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

myOutputCommands = cms.untracked.vstring()
myOutputCommands.extend([
    'drop *', 
    'keep *PHcalTB06Info_*_*_*'
    ])

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(1048576),
##    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = myOutputCommands,
    fileName = cms.untracked.string('file:step02_pi50GeV.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.g4SimHits.UseMagneticField = cms.bool(False)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')


# start particle gun generation
#-------------------------------
process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
##        PartID = cms.vint32(11),
##        PartID = cms.vint32(13),
        PartID = cms.vint32(-211),
#        PartID = cms.vint32(2212),
        MinPt = cms.double(10.00000),
        MaxPt = cms.double(10.00000),
        MinEta = cms.double( 2.400),
        MaxEta = cms.double( 2.400),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double( 3.14159265359)
#        MinPhi = cms.double( 0.200),
#        MaxPhi = cms.double( 0.200)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single pion pt 50'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)

# end generation single particle
#-----------------------------------


# Primary vertex smearing
#-------------------------
from IOMC.EventVertexGenerators.VtxSmearedParameters_cfi import *
process.VtxSmeared = cms.EDProducer("GaussEvtVtxGenerator",
    MeanX = cms.double(0.0),
    MeanY = cms.double(0.0),
    MeanZ = cms.double(0.0),
    SigmaY = cms.double(0.0),
    SigmaX = cms.double(0.0),
    SigmaZ = cms.double(0.0),
    TimeOffset = cms.double(0.0),
    src = cms.InputTag("generator")
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions

process.g4SimHits.Watchers = cms.VPSet(cms.PSet(
    HcalTestAnalysis = cms.PSet(
##//===> necessary for HcalTestAnalysis
         Eta0 = cms.double(0.0),
         Phi0 = cms.double(0.0),
## energy threshold in [MeV] for hit scoring 
##------------------------------------------
         Thresh= cms.double(0.10),
         LayerGrouping = cms.int32(2),
         CentralTower = cms.int32(0),
         FileName = cms.string("test_pi_gun_01.root"),
##//===>
        Names    = cms.vstring('HcalHits', 'EcalHitsEB', 'EcalHitsEE'),
##        Names    = cms.vstring('HcalHits'),
        EHCalMax = cms.untracked.double(2.0),
        ETtotMax = cms.untracked.double(20.0),
        Verbose  = cms.untracked.bool(True)
    ),

#//===> necessary for HcalTestAnalysis
    HcalQie = cms.PSet(
        NumOfBuckets  = cms.int32(10),
        BaseLine      = cms.int32(4),
        BinOfMax      = cms.int32(6),
        PreSamples    = cms.int32(0),
        EDepPerPE     = cms.double(0.0005),
        SignalBuckets = cms.int32(2),
        SigmaNoise    = cms.double(0.5),
        qToPE         = cms.double(4.0)
     ),
#//===>
    type = cms.string('HcalTestAnalysis')
))

