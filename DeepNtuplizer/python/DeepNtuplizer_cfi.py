import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                           vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                           jets       = cms.InputTag("slimmedJets"),
                           rho = cms.InputTag("fixedGridRhoFastjetAll"),
                           CSVcomputer = cms.string('candidateCombinedSecondaryVertexV2Computer'),
                           ipTagInfos = cms.string("pfImpactParameter"),
                           svTagInfos = cms.string("pfInclusiveSecondaryVertexFinder"),
                           bDiscriminators = cms.vstring(),
                           )
