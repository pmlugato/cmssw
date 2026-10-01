// JpsiKCandidateSplitter
//
// Stage-2 Stage-1-output adapter for B+ -> J/psi K+ (and other one-track-
// bachelor channels). Consumes a `VertexCompositeCandidateCollection` of B+
// candidates (each with two daughters: a J/psi VCC and a bachelor
// RecoChargedCandidate) and emits, in the same per-B+ order:
//
//   - a `reco::VertexCompositeCandidateCollection` of dimuon (J/psi) sub-
//     candidates -- shallow copies of daughter(0). Their internal muon-RCC
//     daughters keep pointing at the persisted track collection, so the
//     downstream `ResidualGlobalCorrectionMakerTwoTrackG4e` (`srcCandidates=
//     <this>`, `src=ALCARECOTkAlJpsiX`) resolves tracks via the candidate's
//     fast path (RecoChargedCandidate::track()).
//   - a `reco::TrackCollection` of bachelor tracks -- value-copies of
//     daughter(1)::track() dereferenced. The Track value copy preserves its
//     internal TrackExtraRef into the persisted vector<reco::TrackExtra>, so
//     the single-track maker (`src=<this>`) reaches hits the usual way.
//   - two parallel `std::vector<int>` index branches, in the SAME order as the
//     emitted dimuon / bachelor collections, holding the source B+ candidate
//     index. Per-candidate joining downstream uses these.
//
// Species disambiguation is by the daughter `pdgId()` set at Stage 1 by
// `JpsiXCandidateProducer` (bachelor) / `TwoBodyDecayCandidateProducer`
// (J/psi muons and K*0/phi daughters). This splitter reads the bachelor
// pdgId and the two J/psi-daughter muon pdgIds off the candidate and
// propagates them onto its output as parallel `std::vector<int>` branches
// (`bachelorPdgId`, `muon0PdgId`, `muon1PdgId`). Downstream Stage-2
// mass-hypothesis lookup consumes those branches instead of hard-coding
// per-channel masses.
//
// Collection LAYOUT still defines the daughter ORDER
// (daughter(0) = J/psi VCC, daughter(1) = X): only the species
// TAG moves to the pdgId. See openspec change
// add-jpsi-x-muons-and-preprod-refinements.
//
// N HADRON LEGS (2026-09-30, openspec add-second-jpsix-channel-end-to-end
// task 1.0b). daughter(1) is a single RecoChargedCandidate only for the
// one-bachelor channels (B+, Bc). For B0 -> J/psi K*0(-> K pi) and
// Bs -> J/psi phi(-> K K) it is a nested VCC carrying TWO charged leaves.
// The old code dynamic_cast<RecoChargedCandidate>(daughter(1)) and skipped
// the candidate when that failed, so every four-body channel produced
// nothing at all. daughter(1) is now descended recursively to its charged
// leaves and each leaf is emitted as its own hadron leg.
//
// Legs are ordered BY CHARGE, positive first (`had0` = positive). Charge is
// unambiguous; pt ordering is not stable against resolution, and for
// B0 -> K*0 the K/pi assignment is the whole problem, so the ordering rule
// must not depend on it (design.md D4).
//
// BACKWARD COMPATIBILITY IS LOAD-BEARING. The original products
// (`bachelor`, `bachelorBCandIdx`, `bachelorPdgId`) are still emitted and
// still carry exactly what they carried before for one-leaf candidates, so
// the B+ path is bit-identical. Multi-leaf candidates contribute nothing to
// them. New consumers use the `hadron*` products, which cover both cases.
//
// No fit code, no candidate-level cuts. Stage-1 already selected.

#include <algorithm>
#include <memory>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

class JpsiKCandidateSplitter : public edm::stream::EDProducer<> {
public:
  explicit JpsiKCandidateSplitter(const edm::ParameterSet& iConfig);
  ~JpsiKCandidateSplitter() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> bplusToken_;
  // Diagnostic counters (one per produce() call). Total across the job is
  // assembled by the outer harness when needed; we just emit LogDebug here.
  // Not thread-shared; one instance per stream.
};

JpsiKCandidateSplitter::JpsiKCandidateSplitter(const edm::ParameterSet& iConfig)
    : bplusToken_(consumes<reco::VertexCompositeCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("src"))) {
  produces<reco::VertexCompositeCandidateCollection>("dimuon");
  produces<reco::TrackCollection>("bachelor");
  produces<std::vector<int>>("dimuonBCandIdx");
  produces<std::vector<int>>("bachelorBCandIdx");
  // Species tags read from daughter->pdgId() at Stage 1 (openspec change
  // add-jpsi-x-muons-and-preprod-refinements). Downstream Stage-2 lookup
  // consumes these instead of hard-coding per-channel masses.
  produces<std::vector<int>>("bachelorPdgId");
  produces<std::vector<int>>("muon0PdgId");
  produces<std::vector<int>>("muon1PdgId");
  // N-leg products (task 1.0b). Flat across candidates: one entry per
  // hadron leg, with `hadronBCandIdx` giving the source candidate and
  // `hadronLegIdx` the leg's position within it (0 = positive by charge).
  produces<reco::TrackCollection>("hadron");
  produces<std::vector<int>>("hadronBCandIdx");
  produces<std::vector<int>>("hadronLegIdx");
  produces<std::vector<int>>("hadronPdgId");
  produces<std::vector<int>>("nHadronLegs");  // one entry per emitted candidate
}

namespace {
// Collect the charged leaves under `cand`, descending composite daughters.
// A RecoChargedCandidate IS a leaf and is not descended further.
void collectChargedLeaves(const reco::Candidate* cand,
                          std::vector<const reco::RecoChargedCandidate*>& out,
                          int depth = 0) {
  if (cand == nullptr || depth > 4) {
    return;
  }
  if (const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(cand)) {
    out.push_back(rcc);
    return;
  }
  for (std::size_t i = 0; i < cand->numberOfDaughters(); ++i) {
    collectChargedLeaves(cand->daughter(i), out, depth + 1);
  }
}
}  // namespace

void JpsiKCandidateSplitter::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::VertexCompositeCandidateCollection> bplusH;
  iEvent.getByToken(bplusToken_, bplusH);

  auto outDimuon = std::make_unique<reco::VertexCompositeCandidateCollection>();
  auto outBachelor = std::make_unique<reco::TrackCollection>();
  auto outDimuonIdx = std::make_unique<std::vector<int>>();
  auto outBachelorIdx = std::make_unique<std::vector<int>>();
  auto outBachelorPdg = std::make_unique<std::vector<int>>();
  auto outMuon0Pdg = std::make_unique<std::vector<int>>();
  auto outMuon1Pdg = std::make_unique<std::vector<int>>();
  auto outHadron = std::make_unique<reco::TrackCollection>();
  auto outHadronIdx = std::make_unique<std::vector<int>>();
  auto outHadronLeg = std::make_unique<std::vector<int>>();
  auto outHadronPdg = std::make_unique<std::vector<int>>();
  auto outNLegs = std::make_unique<std::vector<int>>();

  int nSkip = 0;
  outDimuon->reserve(bplusH->size());
  outBachelor->reserve(bplusH->size());
  outDimuonIdx->reserve(bplusH->size());
  outBachelorIdx->reserve(bplusH->size());
  outBachelorPdg->reserve(bplusH->size());
  outMuon0Pdg->reserve(bplusH->size());
  outMuon1Pdg->reserve(bplusH->size());
  outHadron->reserve(2 * bplusH->size());
  outHadronIdx->reserve(2 * bplusH->size());
  outHadronLeg->reserve(2 * bplusH->size());
  outHadronPdg->reserve(2 * bplusH->size());
  outNLegs->reserve(bplusH->size());

  for (std::size_t ib = 0; ib < bplusH->size(); ++ib) {
    const reco::VertexCompositeCandidate& bplus = (*bplusH)[ib];
    if (bplus.numberOfDaughters() < 2) {
      ++nSkip;
      continue;
    }
    // daughter(0) is the J/psi VCC, daughter(1) is the bachelor RCC, per
    // `JpsiXCandidateProducer.cc`. LAYOUT is trusted; species/mass hypothesis
    // is looked up from daughter->pdgId() (Stage-1 non-zero invariant,
    // openspec add-jpsi-x-muons-and-preprod-refinements).
    const auto* dJpsi = dynamic_cast<const reco::VertexCompositeCandidate*>(
        bplus.daughter(0));
    if (!dJpsi) {
      ++nSkip;
      continue;
    }
    if (dJpsi->numberOfDaughters() != 2) {
      ++nSkip;
      continue;
    }

    // daughter(1) is the X: one charged track (B+, Bc) or a nested VCC with
    // two charged leaves (K*0 -> K pi, phi -> K K, Ks -> pi pi, Lambda -> p pi).
    std::vector<const reco::RecoChargedCandidate*> leaves;
    collectChargedLeaves(bplus.daughter(1), leaves);
    if (leaves.empty()) {
      ++nSkip;
      continue;
    }
    bool allHaveTracks = true;
    for (const auto* l : leaves) {
      if (l->track().isNull()) {
        allHaveTracks = false;
        break;
      }
    }
    if (!allHaveTracks) {
      ++nSkip;
      continue;
    }

    // Order legs by charge, positive first. Stable for equal charges, so a
    // same-sign pair keeps its Stage-1 order rather than being reordered
    // arbitrarily run to run.
    std::stable_sort(leaves.begin(), leaves.end(),
                     [](const reco::RecoChargedCandidate* a,
                        const reco::RecoChargedCandidate* b) {
                       return a->charge() > b->charge();
                     });

    outDimuon->push_back(*dJpsi);
    outDimuonIdx->push_back(static_cast<int>(ib));
    outMuon0Pdg->push_back(dJpsi->daughter(0)->pdgId());
    outMuon1Pdg->push_back(dJpsi->daughter(1)->pdgId());

    // Legacy one-bachelor products: populated only for single-leaf
    // candidates, with exactly the content they had before, so the B+ path
    // is bit-identical across this change.
    if (leaves.size() == 1) {
      outBachelor->push_back(*leaves[0]->track());
      outBachelorIdx->push_back(static_cast<int>(ib));
      outBachelorPdg->push_back(leaves[0]->pdgId());
    }

    // N-leg products: always populated, one entry per leg.
    for (std::size_t il = 0; il < leaves.size(); ++il) {
      outHadron->push_back(*leaves[il]->track());
      outHadronIdx->push_back(static_cast<int>(ib));
      outHadronLeg->push_back(static_cast<int>(il));
      outHadronPdg->push_back(leaves[il]->pdgId());
    }
    outNLegs->push_back(static_cast<int>(leaves.size()));
  }

  if (nSkip > 0) {
    LogDebug("JpsiKCandidateSplitter")
        << "skipped " << nSkip << " B+ candidates with unexpected daughter layout"
        << " out of " << bplusH->size() << " total";
  }

  iEvent.put(std::move(outDimuon), "dimuon");
  iEvent.put(std::move(outBachelor), "bachelor");
  iEvent.put(std::move(outDimuonIdx), "dimuonBCandIdx");
  iEvent.put(std::move(outBachelorIdx), "bachelorBCandIdx");
  iEvent.put(std::move(outBachelorPdg), "bachelorPdgId");
  iEvent.put(std::move(outMuon0Pdg), "muon0PdgId");
  iEvent.put(std::move(outMuon1Pdg), "muon1PdgId");
  iEvent.put(std::move(outHadron), "hadron");
  iEvent.put(std::move(outHadronIdx), "hadronBCandIdx");
  iEvent.put(std::move(outHadronLeg), "hadronLegIdx");
  iEvent.put(std::move(outHadronPdg), "hadronPdgId");
  iEvent.put(std::move(outNLegs), "nHadronLegs");
}

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiKCandidateSplitter;
DEFINE_FWK_MODULE(JpsiKCandidateSplitter);
