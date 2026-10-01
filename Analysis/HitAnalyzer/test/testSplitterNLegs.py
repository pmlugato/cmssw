"""Minimal check that JpsiKCandidateSplitter emits N hadron legs.

Task 1.0b gate, openspec `add-second-jpsix-channel-end-to-end`. Runs the
splitter on one AlCaReco file for a configurable channel and writes only the
splitter's own products, so the leg counts can be read back with uproot
without a full Stage-2 job.

  cmsRun testSplitterNLegs.py channel=B0Kstar maxEvents=2000
"""
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

opts = VarParsing.VarParsing("analysis")
opts.register("channel", "BPlus", VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              "JpsiX channel suffix, e.g. BPlus / B0Kstar / BsPhi / Bc")
opts.register("nEvents", 2000, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, "events to read")
opts.register("infile",
              "/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/"
              "charmonium/Run2016G/TkAlJpsiX_00033B23-B475-E611-B79A-FA163ED0327F.root",
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, "AlCaReco input")
opts.register("out", "splitter_nlegs.root",
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, "output")
opts.parseArguments()

process = cms.Process("SPLITTEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
process.source = cms.Source(
    "PoolSource", fileNames=cms.untracked.vstring("file:" + opts.infile))

process.split = cms.EDProducer(
    "JpsiKCandidateSplitter",
    src=cms.InputTag("ALCARECOTkAlJpsiX%sResonances" % opts.channel),
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(opts.out),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_split_*_*",
    ),
)

process.p = cms.Path(process.split)
process.e = cms.EndPath(process.out)
