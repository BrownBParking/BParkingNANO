#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class DToKstarPiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DToKstarPiBuilder(const edm::ParameterSet &cfg):
    pi_selection_{cfg.getParameter<std::string>("piSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    Kstar_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("Kstar") )},
    Kstar_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("KstarTransientTracks") )},
    pions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions") )},
    pions_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("pionsTransientTracks") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~DToKstarPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> Kstar_;
  const edm::EDGetTokenT<TransientTrackCollection> Kstar_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  

};

void DToKstarPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> Kstar;
  evt.getByToken(Kstar_, Kstar);
  
  edm::Handle<TransientTrackCollection> Kstar_ttracks;
  evt.getByToken(Kstar_ttracks_, Kstar_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> D_out(new pat::CompositeCandidateCollection());
 
  for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
    edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);
    if( !pi_selection_(*pi_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector pi_p4(
      pi_ptr->pt(), 
      pi_ptr->eta(),
      pi_ptr->phi(),
      PI_MASS
      );

    for(size_t Kstar_idx = 0; Kstar_idx < Kstar->size(); ++Kstar_idx) {
      edm::Ptr<pat::CompositeCandidate> Kstar_ptr(Kstar, Kstar_idx);
      edm::Ptr<reco::Candidate> trk1_ptr = Kstar_ptr->userCand("trk1");
      edm::Ptr<reco::Candidate> trk2_ptr = Kstar_ptr->userCand("trk2");
      
      if (trk1_ptr->charge() == pi_ptr->charge()) continue;

      int trk1_idx = Kstar_ptr->userInt("trk1_idx");
      int trk2_idx = Kstar_ptr->userInt("trk2_idx");
      
      if (((int)pi_idx == trk1_idx  ) | ((int)pi_idx == trk2_idx)) continue;
      pat::CompositeCandidate D_cand;
      D_cand.setP4(Kstar_ptr->p4() + pi_p4);
      D_cand.setCharge(Kstar_ptr->charge() + pi_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the trkton passing the corresponding selection
      D_cand.addUserCand("trk1", trk1_ptr);
      D_cand.addUserCand("trk2", trk2_ptr);
      D_cand.addUserCand("Pi", pi_ptr);
      D_cand.addUserCand("Kstar", Kstar_ptr);

      D_cand.addUserInt("trk1_idx", trk1_idx);
      D_cand.addUserInt("trk2_idx", trk2_idx);
      D_cand.addUserInt("pi_idx", pi_idx);
   
      if( !pre_vtx_selection_(D_cand) ) continue;
    
      KinVtxFitter fitter(
        {Kstar_ttracks->at(trk1_idx), Kstar_ttracks->at(trk2_idx), pions_ttracks->at(pi_idx)},
        {K_MASS, PI_MASS, PI_MASS},
        {K_SIGMA, K_SIGMA, K_SIGMA} //some small sigma for the trk mass
        );
      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?

      auto fit_p4 = fitter.fitted_p4();

      D_cand.addUserFloat("sv_chi2", fitter.chi2());
      D_cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      D_cand.addUserFloat("sv_prob", fitter.prob());
      D_cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      D_cand.addUserFloat("fitted_eta" , fit_p4.eta());
      D_cand.addUserFloat("fitted_phi" , fit_p4.phi());
      D_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      

      D_cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, D_cand.p4())
        );
      D_cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
      auto lxy = l_xy(fitter, *beamspot);
      D_cand.addUserFloat("l_xy", lxy.value());
      D_cand.addUserFloat("l_xy_unc", lxy.error());

      if( !post_vtx_selection_(D_cand) ) continue;        
      D_out->push_back(D_cand);
    }
  }
  evt.put(std::move(D_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DToKstarPiBuilder);
