import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("AOD2NanoAOD")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("AOD2NanoAOD")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Set the maximum number of events to be processed (-1 processes all events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Define files of dataset
#files = FileUtils.loadListFromFile("data/CMS_MonteCarlo2012_Summer12_DR53X_TTbar_8TeV-Madspin_aMCatNLO-herwig_AODSIM_PU_S10_START53_V19-v2_00000_file_index.txt")
#files.extend(FileUtils.loadListFromFile("data/CMS_MonteCarlo2012_Summer12_DR53X_TTbar_8TeV-Madspin_aMCatNLO-herwig_AODSIM_PU_S10_START53_V19-v2_20000_file_index.txt"))

#process.source = cms.Source(
#   "PoolSource", fileNames=cms.untracked.vstring(*files))

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000A9D3F-CE4C-E311-84F8-001E673969D2.root'))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))

# Set global tag
# We don't have set the global tag for the educational samples. This simplifies running the code since we don't have to access the database.
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V27::All"

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.ak5PFJets.doAreaFastjet = True

process.ak5PFchsCorrectedJets = cms.EDProducer('CorrectedPFJetProducer',
        src = cms.InputTag("ak5PFJets"),
        correctors  = cms.VInputTag('ak5PFCHSL1FastL2L3Corrector')
        )

process.ak5PFchsCorrectedJetsSmeared = cms.EDProducer('SmearedPFJetProducer',
        src = cms.InputTag('ak5PFchsCorrectedJets'),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK5PFchs'),
        algopt = cms.string('AK5PFchs_pt')
        )


# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

# Register fileservice for output file
process.aod2nanoaod = cms.EDAnalyzer("AOD2NanoAOD", 
        plain_jets = cms.InputTag('slimmedJets'),
        corrected_jets = cms.InputTag('ak5PFchsCorrectedJets'),
        smeared_jets = cms.InputTag('ak5PFchsCorrectedJetsSmeared'), 
        isData = cms.bool(False)
        )
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("output.root"))

process.p = cms.Path(process.aod2nanoaod)
