#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

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

class BsToPiDsBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BsToPiDsBuilder(const edm::ParameterSet &cfg):
    pi_selection_{cfg.getParameter<std::string>("piSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    Ds_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("Ds") )},
    Ds_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("DsTransientTracks") )},
    pions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions") )},
    pions_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("pionsTransientTracks") )},
    //isotracksToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    //vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BsToPiDsBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> Ds_;
  const edm::EDGetTokenT<TransientTrackCollection> Ds_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;
  
  //const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  //const edm::EDGetTokenT<pat::CompositeCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BsToPiDsBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> Ds;
  evt.getByToken(Ds_, Ds);
  
  edm::Handle<TransientTrackCollection> Ds_ttracks;
  evt.getByToken(Ds_ttracks_, Ds_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  edm::Handle<pat::MuonCollection> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);

  //edm::Handle<reco::VertexCollection> vertexHandle;
  //evt.getByToken(vertexToken_, vertexHandle);
  //const reco::Vertex & PV = vertexHandle->front();
  //math::XYZPoint pv(PV.x(), PV.y(), PV.z());

  //edm::ESHandle<MagneticField> bFieldHandle;
  //stp.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  //for isolation
  //edm::Handle<pat::CompositeCandidateCollection> iso_tracks;
  //evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_trk1_id, used_trk2_id, used_trk3_id, used_pi_id;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
    edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);
    if( !pi_selection_(*pi_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector pi_p4(
      pi_ptr->pt(), 
      pi_ptr->eta(),
      pi_ptr->phi(),
      PI_MASS
      );

    for(size_t Ds_idx = 0; Ds_idx < Ds->size(); ++Ds_idx) {
      edm::Ptr<pat::CompositeCandidate> Ds_ptr(Ds, Ds_idx);

      //edm::Ptr<reco::Candidate> trk1_ptr = Ds_ptr->userCand("trk1");
      //edm::Ptr<reco::Candidate> trk2_ptr = Ds_ptr->userCand("trk2");
      //edm::Ptr<reco::Candidate> trk3_ptr = Ds_ptr->userCand("trk3");
      unsigned int trk1_idx = Ds_ptr->userInt("trk1_idx");
      unsigned int trk2_idx = Ds_ptr->userInt("trk2_idx");
      unsigned int trk3_idx = Ds_ptr->userInt("trk3_idx");

      edm::Ptr<pat::CompositeCandidate> trk1_ptr(pions, trk1_idx);
      edm::Ptr<pat::CompositeCandidate> trk2_ptr(pions, trk2_idx);
      edm::Ptr<pat::CompositeCandidate> trk3_ptr(pions, trk3_idx);

      if (trk3_ptr->charge() == pi_ptr->charge()) continue;
      if (( pi_idx == trk1_idx) | ( pi_idx == trk2_idx) | ( pi_idx == trk3_idx)) continue;
      
      pat::CompositeCandidate cand;
      cand.setP4(Ds_ptr->p4() + pi_p4);
      cand.setCharge(Ds_ptr->charge() + pi_ptr->charge());
      //std::cout<<"trk1 idx : "<<trk1_idx<<", charge : "<<trk1_ptr->charge()<<std::endl;
      //std::cout<<"trk2 idx : "<<trk2_idx<<", charge : "<<trk2_ptr->charge()<<std::endl;
      //std::cout<<"trk3 idx : "<<trk3_idx<<", charge : "<<trk3_ptr->charge()<<std::endl;
      //std::cout<<"pi idx : "<<pi_idx<<", charge : "<<pi_ptr->charge()<<std::endl;

      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the trkton passing the corresponding selection
      cand.addUserCand("trk1", trk1_ptr);
      cand.addUserCand("trk2", trk2_ptr);
      cand.addUserCand("trk3", trk3_ptr);
      cand.addUserCand("Pi", pi_ptr);
      cand.addUserCand("Ds", Ds_ptr);

      cand.addUserInt("trk1_idx", trk1_idx);
      cand.addUserInt("trk2_idx", trk2_idx);
      cand.addUserInt("trk3_idx", trk3_idx);
      cand.addUserInt("pi_idx", pi_idx);
    
      auto dr_info = min_max_dr({trk1_ptr, trk2_ptr, pi_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);

      //second mass hypothesis
      auto Ds_barP4 = Ds_ptr->polarP4();
      Ds_barP4.SetM(Ds_ptr->userFloat("barMass"));
      cand.addUserFloat("barMass",(pi_ptr->polarP4()+Ds_barP4).M() );

      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
    
      KinVtxFitter fitter(
        {Ds_ttracks->at(trk1_idx), Ds_ttracks->at(trk2_idx), Ds_ttracks->at(trk3_idx), pions_ttracks->at(pi_idx)},
        {K_MASS, K_MASS, PI_MASS, PI_MASS},
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
      used_pi_id.emplace_back(pi_idx);

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
      cand.addUserFloat("fitted_pi_pt"  , fitter.daughter_p4(3).pt()); 
      cand.addUserFloat("fitted_pi_eta" , fitter.daughter_p4(3).eta());
      cand.addUserFloat("fitted_pi_phi" , fitter.daughter_p4(3).phi());

      // second mass hypothesis
      auto trk1p4 = fitter.daughter_p4(0);
      auto trk2p4 = fitter.daughter_p4(1);
      auto trk3p4 = fitter.daughter_p4(2);
      trk1p4.SetM(K_MASS);
      trk2p4.SetM(PI_MASS);
      trk3p4.SetM(K_MASS);
      cand.addUserFloat("barMassDs_fullfit",(trk1p4+trk2p4+trk3p4).M());
      cand.addUserFloat("fitted_barMass",(trk1p4+trk2p4+trk3p4+fitter.daughter_p4(3)).M());     
   
      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float trk1_iso03 = 0;
      float trk1_iso04 = 0;
      float trk2_iso03 = 0;
      float trk2_iso04 = 0;
      float trk3_iso03 = 0;
      float trk3_iso04 = 0;
      float pi_iso03  = 0;
      float pi_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {

        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        //edm::Ptr<pat::CompositeCandidate> trk(iso_tracks, iTrk);
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
        if (pi_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two trks
        if (trk1_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk2_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk3_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        //if (pi_idx == iTrk) continue;
        //if (trk1_idx == iTrk) continue;
        //if (trk2_idx == iTrk) continue;
        //if (trk3_idx == iTrk) continue;

        // add to final particle iso if dR < cone
        float dr_to_trk1 = deltaR(cand.userFloat("fitted_trk1_eta"), cand.userFloat("fitted_trk1_phi"), trk.eta(), trk.phi());
        float dr_to_trk2 = deltaR(cand.userFloat("fitted_trk2_eta"), cand.userFloat("fitted_trk2_phi"), trk.eta(), trk.phi());
        float dr_to_trk3 = deltaR(cand.userFloat("fitted_trk3_eta"), cand.userFloat("fitted_trk3_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_pi_eta") , cand.userFloat("fitted_pi_phi") , trk.eta(), trk.phi());
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
          pi_iso04 += trk.pt();
          if (dr_to_k < 0.3) pi_iso03 += trk.pt();
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
      cand.addUserFloat("pi_iso03" , pi_iso03 );
      cand.addUserFloat("pi_iso04" , pi_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );

      ret_val->push_back(cand);
    } // for(size_t Ds_idx = 0; Ds_idx < ditrktons->size(); ++Ds_idx) {
  } // for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx)

  for (auto & cand: *ret_val){
    cand.addUserInt("n_pi_used", std::count(used_pi_id.begin(),used_pi_id.end(),cand.userInt("pi_idx")));
    cand.addUserInt("n_trk1_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk1_idx"))+std::count(used_trk2_id.begin(),used_trk2_id.end(),cand.userInt("trk1_idx")));
    cand.addUserInt("n_trk2_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk2_idx"))+std::count(used_trk2_id.begin(),used_trk2_id.end(),cand.userInt("trk2_idx")));
    cand.addUserInt("n_trk3_used", std::count(used_trk1_id.begin(),used_trk1_id.end(),cand.userInt("trk3_idx"))+std::count(used_trk3_id.begin(),used_trk3_id.end(),cand.userInt("trk3_idx")));

  }

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BsToPiDsBuilder);
