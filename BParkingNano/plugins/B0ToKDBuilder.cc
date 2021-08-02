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

class B0ToKDBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit B0ToKDBuilder(const edm::ParameterSet &cfg):
    K_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    D_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("D") )},
    D_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("DTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~B0ToKDBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> K_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> D_;
  const edm::EDGetTokenT<TransientTrackCollection> D_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;
 
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void B0ToKDBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> D;
  evt.getByToken(D_, D);
  
  edm::Handle<TransientTrackCollection> D_ttracks;
  evt.getByToken(D_ttracks_, D_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  edm::Handle<pat::MuonCollection> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);


  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_trk1_id, used_trk2_id, used_trk3_id, used_K_id;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  for(size_t K_idx = 0; K_idx < kaons->size(); ++K_idx) {
    edm::Ptr<pat::CompositeCandidate> K_ptr(kaons, K_idx);
    if( !K_selection_(*K_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector K_p4(
      K_ptr->pt(), 
      K_ptr->eta(),
      K_ptr->phi(),
      K_MASS
      );

    for(size_t D_idx = 0; D_idx < D->size(); ++D_idx) {
      edm::Ptr<pat::CompositeCandidate> D_ptr(D, D_idx);
      //edm::Ptr<reco::Candidate> trk1_ptr = D_ptr->userCand("trk1");
      //edm::Ptr<reco::Candidate> trk2_ptr = D_ptr->userCand("trk2");
      //edm::Ptr<reco::Candidate> trk3_ptr = D_ptr->userCand("trk3");
      unsigned int trk1_idx = D_ptr->userInt("trk1_idx");
      unsigned int trk2_idx = D_ptr->userInt("trk2_idx");
      unsigned int trk3_idx = D_ptr->userInt("trk3_idx");

      edm::Ptr<pat::CompositeCandidate> trk1_ptr(kaons, trk1_idx);
      edm::Ptr<pat::CompositeCandidate> trk2_ptr(kaons, trk2_idx);
      edm::Ptr<pat::CompositeCandidate> trk3_ptr(kaons, trk3_idx);
 
      if (trk3_ptr->charge() == K_ptr->charge()) continue;
 
      if ((K_idx == trk1_idx)|(K_idx == trk2_idx)|(K_idx == trk3_idx)) continue;   
      pat::CompositeCandidate cand;
      cand.setP4(D_ptr->p4() + K_p4);
      cand.setCharge(D_ptr->charge() + K_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the trkton passing the corresponding selection
      cand.addUserCand("trk1", trk1_ptr);
      cand.addUserCand("trk2", trk2_ptr);
      cand.addUserCand("trk3", trk3_ptr);
      cand.addUserCand("K", K_ptr);
      cand.addUserCand("D", D_ptr);

      cand.addUserInt("trk1_idx", trk1_idx);
      cand.addUserInt("trk2_idx", trk2_idx);
      cand.addUserInt("trk3_idx", trk3_idx);
      cand.addUserInt("K_idx", K_idx);
    
      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
    
      KinVtxFitter fitter(
        {D_ttracks->at(trk1_idx), D_ttracks->at(trk2_idx), D_ttracks->at(trk3_idx), kaons_ttracks->at(K_idx)},
        {K_MASS, PI_MASS, PI_MASS, K_MASS},
        {K_SIGMA, K_SIGMA, K_SIGMA, K_SIGMA} //some small sigma for the trk mass
        );
      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      used_trk1_id.emplace_back(trk1_idx);
      used_trk2_id.emplace_back(trk2_idx);
      used_trk3_id.emplace_back(trk3_idx);
      used_K_id.emplace_back(K_idx);

      auto fit_p4 = fitter.fitted_p4();

      float max_dr = 0.;
      for (const pat::Muon & mu: *trgMuons){
        //std::cout<<"dz wrt PV : "<<mu.vz() - PV.z()<<std::endl;

        //remove tracks inside trg muons jet
        float dR = reco::deltaR(cand,mu);
        if (dR > max_dr) max_dr = dR;
      }

      cand.addUserFloat("sv_dRtrgMu", max_dr);
      cand.addUserInt("sv_OK" , fitter.success());
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());
      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));

      cand.addUserFloat("fitted_trk1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_trk1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_trk1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("fitted_trk2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_trk2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_trk2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("fitted_trk3_pt" , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_trk3_eta", fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_trk3_phi", fitter.daughter_p4(2).phi());
      cand.addUserFloat("fitted_K_pt"  , fitter.daughter_p4(3).pt()); 
      cand.addUserFloat("fitted_K_eta" , fitter.daughter_p4(3).eta());
      cand.addUserFloat("fitted_K_phi" , fitter.daughter_p4(3).phi());
   
      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float trk1_iso03 = 0;
      float trk1_iso04 = 0;
      float trk2_iso03 = 0;
      float trk2_iso04 = 0;
      float trk3_iso03 = 0;
      float trk3_iso04 = 0;
      float K_iso03  = 0;
      float K_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
 
        //edm::Ptr<pat::CompositeCandidate> trk(iso_tracks, iTrk);  
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;

        bool skipTrack=true;
        for (const pat::Muon & mu: *trgMuons){
          //std::cout<<"dz wrt PV : "<<mu.vz() - PV.z()<<std::endl;

          //remove tracks inside trg muons jet
          if(reco::deltaR(trk, mu) < 0.4) 
            continue;
          //if dz is negative it is deactivated
          if((fabs(trk.vz() - mu.vz()) > 0.5))
            continue;
          skipTrack=false;
          break; // at least for one trg muon to pass this cuts
        }
        // if track is closer to at least a triggering muon keep it
        if (skipTrack) continue;


        // check if the track is the pion

        if (K_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two trks
        if (trk1_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk2_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk3_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        //if ((abs(trk.pdgId()) != 211) | (abs(trk.pdgId()) != 130)) continue;
        //if (K_idx == iTrk) continue;
        //if (trk1_idx == iTrk) continue;
        //if (trk2_idx == iTrk) continue;
        //if (trk3_idx == iTrk) continue;

        // add to final particle iso if dR < cone
        float dr_to_trk1 = deltaR(cand.userFloat("fitted_trk1_eta"), cand.userFloat("fitted_trk1_phi"), trk.eta(), trk.phi());
        float dr_to_trk2 = deltaR(cand.userFloat("fitted_trk2_eta"), cand.userFloat("fitted_trk2_phi"), trk.eta(), trk.phi());
        float dr_to_trk3 = deltaR(cand.userFloat("fitted_trk3_eta"), cand.userFloat("fitted_trk3_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_K_eta") , cand.userFloat("fitted_K_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_trk1 < 0.4){
          trk1_iso04 += trk.pt();
          if ( dr_to_trk1 < 0.3) trk1_iso03 += trk.pt();
        }
        if (dr_to_trk2 < 0.4){
          trk2_iso04 += trk.pt();
          if (dr_to_trk2 < 0.3)  trk2_iso03 += trk.pt();
        }
        if (dr_to_trk3 < 0.4){
          trk3_iso04 += trk.pt();
          if (dr_to_trk3 < 0.3)  trk3_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          K_iso04 += trk.pt();
          if (dr_to_k < 0.3) K_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      }
      cand.addUserFloat("trk1_iso03", trk1_iso03);
      cand.addUserFloat("trk1_iso04", trk1_iso04);
      cand.addUserFloat("trk2_iso03", trk2_iso03);
      cand.addUserFloat("trk2_iso04", trk2_iso04);
      cand.addUserFloat("trk3_iso03", trk3_iso03);
      cand.addUserFloat("trk3_iso04", trk3_iso04);
      cand.addUserFloat("K_iso03" , K_iso03 );
      cand.addUserFloat("K_iso04" , K_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );

      ret_val->push_back(cand);
    } // for(size_t D_idx = 0; D_idx < ditrktons->size(); ++D_idx) {
  } // for(size_t K_idx = 0; K_idx < kaons->size(); ++K_idx)

  for (auto & cand: *ret_val){
    cand.addUserInt("n_K_used", std::count(used_K_id.begin(),used_K_id.end(),cand.userInt("K_idx")));
    cand.addUserInt("n_trk1_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk1_idx"))+std::count(used_trk2_id.begin(),used_trk2_id.end(),cand.userInt("trk1_idx")));
    cand.addUserInt("n_trk2_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk2_idx"))+std::count(used_trk2_id.begin(),used_trk2_id.end(),cand.userInt("trk2_idx")));
    cand.addUserInt("n_trk3_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk3_idx"))+std::count(used_trk3_id.begin(),used_trk3_id.end(),cand.userInt("trk3_idx")));

  }

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(B0ToKDBuilder);
