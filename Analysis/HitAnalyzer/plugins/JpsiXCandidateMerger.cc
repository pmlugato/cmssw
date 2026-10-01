// JpsiXCandidateMerger
//
// Concatenates several J/psi+X Stage-1 candidate collections into one, and
// emits a parallel `std::vector<int> channelId` saying which input each
// candidate came from.
//
// WHY (openspec add-second-jpsix-channel-end-to-end, design.md D9d).
// The Stage-2 chain is expensive almost entirely because of the Geant4e CVH
// refit, and the two muons of a J/psi are the SAME tracks for every channel
// that has a candidate in that event. Running one chain per channel refits
// them once per channel; running one chain over the merged collection refits
// them once per event. Muons are 2 of 3 legs for B+ and 2 of 4 for
// B0 -> J/psi K*0, so at two channels this is a 33-50% saving and it grows
// with the channel count.
//
// The stronger reason is consistency, not cost. The joint fit correlates
// momentum-scale nuisances ACROSS channels by name, so the muon treatment has
// to be identical between them. Merging makes that a tautology -- one module
// instance, one call, one refit that every channel reads -- rather than
// something that has to be kept true by configuration discipline.
//
// The alternative design (N parallel chains sharing a deduplicated track
// collection) needs every channel's (candIdx, leafIdx) remapped into the
// merged refit output. That is the CompositeDaughterTrackProducer +
// VertexCompositeCandidateRemapper pattern from the AlCaReco step, and it is
// the most intricate part of that step. This plugin replaces it.
//
// WHAT IS NOT DONE HERE.
//   * No deduplication. A candidate appearing in two channels is a genuinely
//     different candidate (different mass hypotheses, different daughter
//     structure), so both rows are kept. Deduplication at the TRACK level is
//     the refit's business, not the candidate collection's.
//   * No selection. "Every selection lives at histmaker level"
//     (add-nbody-cvh-maker). This plugin only concatenates.
//   * No reordering. Candidates appear in input order, inputs in the order
//     given, so `channelId` is non-decreasing and a channel's block is
//     contiguous.
//
// `channelId` indexes into `srcLabels` (a parallel, purely informational
// vector of names) so a reader can recover the channel without parsing
// InputTags. An input that is missing from the event contributes nothing and
// is not an error: channels are configured as a fixed list, but a given
// AlCaReco may legitimately lack one.

#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

class JpsiXCandidateMerger : public edm::stream::EDProducer<> {
public:
  explicit JpsiXCandidateMerger(const edm::ParameterSet&);
  ~JpsiXCandidateMerger() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  std::vector<edm::EDGetTokenT<reco::VertexCompositeCandidateCollection>> tokens_;
  std::vector<std::string> labels_;
};

JpsiXCandidateMerger::JpsiXCandidateMerger(const edm::ParameterSet& iConfig) {
  const auto srcs = iConfig.getParameter<std::vector<edm::InputTag>>("src");
  for (const auto& tag : srcs) {
    tokens_.push_back(consumes<reco::VertexCompositeCandidateCollection>(tag));
    labels_.push_back(tag.label() + (tag.instance().empty()
                                         ? std::string()
                                         : ":" + tag.instance()));
  }
  produces<reco::VertexCompositeCandidateCollection>();
  produces<std::vector<int>>("channelId");
  // Same information as a ValueMap keyed to the merged collection, because
  // the NanoAOD flat-table producer's externalVariables take ValueMaps, not
  // positional vectors. Both are emitted: the vector for positional
  // consumers (the makers' bCandIdxSrc idiom), the ValueMap for the table.
  produces<edm::ValueMap<int>>("channelIdMap");
  produces<std::vector<std::string>>("srcLabels");
}

void JpsiXCandidateMerger::produce(edm::Event& iEvent, const edm::EventSetup&) {
  auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();
  auto chan = std::make_unique<std::vector<int>>();

  for (std::size_t ic = 0; ic < tokens_.size(); ++ic) {
    edm::Handle<reco::VertexCompositeCandidateCollection> h;
    iEvent.getByToken(tokens_[ic], h);
    if (!h.isValid()) {
      continue;  // a channel absent from this AlCaReco is not an error
    }
    for (const auto& cand : *h) {
      out->push_back(cand);
      chan->push_back(static_cast<int>(ic));
    }
  }

  LogDebug("JpsiXCandidateMerger")
      << "merged " << out->size() << " candidates from " << tokens_.size()
      << " inputs";

  // Keep a copy before the unique_ptr is moved away: the ValueMap filler
  // needs both the values and the OrphanHandle of the collection they key to.
  const std::vector<int> chanVals = *chan;
  const auto outH = iEvent.put(std::move(out));

  auto vm = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler filler(*vm);
  filler.insert(outH, chanVals.begin(), chanVals.end());
  filler.fill();

  iEvent.put(std::move(vm), "channelIdMap");
  iEvent.put(std::move(chan), "channelId");
  iEvent.put(std::make_unique<std::vector<std::string>>(labels_), "srcLabels");
}

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiXCandidateMerger;
DEFINE_FWK_MODULE(JpsiXCandidateMerger);
