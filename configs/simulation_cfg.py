import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from RecoMET.METProducers.METSigParams_cfi import *
import os 

relBase = os.environ['CMSSW_BASE']

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

process.ak5PFCorrectedJets = cms.EDProducer('PFJetCorrectionProducer',
        src = cms.InputTag("ak5PFJets"),
        correctors  = cms.vstring('ak5PFL2L3')
        )

process.ak5PFCorrectedJetsSmeared = cms.EDProducer('SmearedPFJetProducer',
        src = cms.InputTag('ak5PFCorrectedJets'),
        jetCorrLabel = cms.string("ak4PFL1FastL2L3Corrector"),
        dRmaxGenJetMatch = cms.string('min(0.5, 0.1 + 0.3*exp(-0.05*(genJetPt - 10.)))'),
        sigmaMaxGenJetMatch = cms.double(3.),
        inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
        lutName = cms.string('pfJetResolutionMCtoDataCorrLUT'),
        jetResolutions = METSignificance_params,
        skipRawJetPtThreshold = cms.double(10.), # GeV
        skipCorrJetPtThreshold = cms.double(1.e-2),
        srcGenJets = cms.InputTag('ak5GenJets'),
        shiftBy = cms.double(101.)
        )

#process.ak5PFCorrectedJetsSmeared = cms.EDProducer('PATPFJetMETcorrInputProducer',
#        src = cms.InputTag('ak5PFCorrectedJets'),
#        offsetCorrLabel = cms.string('L1FastJet'),
#        jetCorrLabel = cms.string("L3Absolute"),
#        type1JetPtThreshold = cms.double(10.0),
#        skipEM = cms.bool(True),
#        skipEMfractionThreshold = cms.double(0.90),
#        skipMuons = cms.bool(True),
#        skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
#        )

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

## Output Module Configuration (expects a path 'p')
process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('jet_corr_name.root'),
        #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
        outputCommands = cms.untracked.vstring('keep *')
        )

# Register fileservice for output file
process.aod2nanoaod = cms.EDAnalyzer("AOD2NanoAOD", 
        plain_jets = cms.InputTag('slimmedJets'),
        corrected_jets = cms.InputTag('ak5PFCorrectedJets'),
        smeared_jets = cms.InputTag('ak5PFCorrectedJetsSmeared'), 
        isData = cms.bool(False)
        )
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("output.root"))

#process.p = cms.Path(process.aod2nanoaod)
process.p = cms.Path(process.ak5PFCorrectedJets * process.ak5PFCorrectedJetsSmeared)# * process.aod2nanoaod)
process.ep = cms.EndPath(process.out)
