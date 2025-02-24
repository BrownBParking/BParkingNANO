///////////////////////// Code to produce K* candidates ////////////////////////


#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"




class DsBuilder : public edm::global::EDProducer<> {

  
public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit DsBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    trk3_selection_{cfg.getParameter<std::string>("trk3Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {

      //output
       produces<pat::CompositeCandidateCollection>();

    }

  ~DsBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk3_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; //input PF cands this is sorted in pT in previous step
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_; //input TTracks of PF cands
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};


void DsBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //std::cout<<"beginning of Ds selection"<<std::endl;

  //inputs  
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
 
  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> Ds_out(new pat::CompositeCandidateCollection());

  
  // main loop
  for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ){

    edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
    if(!trk1_selection_(*trk1_ptr)) continue; 
    
    for(size_t trk2_idx = 0; trk2_idx < pfcands->size(); ++trk2_idx) {

     if (trk1_idx == trk2_idx) continue;
     edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
     if (trk1_ptr->charge() == trk2_ptr->charge()) continue;
     if(!trk2_selection_(*trk2_ptr)) continue;


     for(size_t trk3_idx = 0; trk3_idx < pfcands->size(); ++trk3_idx) {

      if ((trk1_idx == trk3_idx) | (trk2_idx == trk3_idx)) continue;
      edm::Ptr<pat::CompositeCandidate> trk3_ptr( pfcands, trk3_idx );
      if (trk1_ptr->charge() == trk3_ptr->charge()) continue;
      if(!trk3_selection_(*trk3_ptr)) continue;

      //std::cout<<"got 2-track pair"<<std::endl; 
      // create a Ds candidate; add first quantities that can be used for pre fit selection
      pat::CompositeCandidate Ds_cand;
      auto trk1_p4=trk1_ptr->polarP4();
      auto trk2_p4=trk2_ptr->polarP4();
      auto trk3_p4=trk3_ptr->polarP4();

      trk1_p4.SetM(K_MASS);
      trk2_p4.SetM(K_MASS);
      trk3_p4.SetM(PI_MASS);

      //adding stuff for pre fit selection
      Ds_cand.setP4(trk1_p4 + trk2_p4 + trk3_p4);
      Ds_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));

      // save indices
      Ds_cand.addUserInt("trk1_idx", trk1_idx );
      Ds_cand.addUserInt("trk2_idx", trk2_idx );
      Ds_cand.addUserInt("trk3_idx", trk3_idx );

      // save cands      
      Ds_cand.addUserCand("trk1", trk1_ptr );
      Ds_cand.addUserCand("trk2", trk2_ptr );
      Ds_cand.addUserCand("trk3", trk3_ptr );

      auto leadTrkPt = std::max({trk1_p4.pt(), trk2_p4.pt(), trk3_p4.pt()});
      Ds_cand.addUserFloat("leadTrkPt", leadTrkPt);

      //second mass hypothesis
      trk1_p4.SetM(K_MASS);
      trk2_p4.SetM(PI_MASS);
      trk3_p4.SetM(K_MASS);

      Ds_cand.addUserFloat("barMass", (trk1_p4 + trk2_p4 + trk3_p4).M() );
 
      // selection before fit
      if( !pre_vtx_selection_(Ds_cand) ) continue;
      //std::cout<<"got pre selected candidiate"<<std::endl;    
      
      KinVtxFitter fitter(
        {ttracks->at(trk1_idx), ttracks->at(trk2_idx), ttracks->at(trk3_idx)},
        { K_MASS, K_MASS, PI_MASS },
        {K_SIGMA, K_SIGMA, 2.4e-7} //K and PI sigma equal...
         );
      // std::cout<<fitter.success()<<std::endl;
       if ( !fitter.success() ) continue;           

       auto fit_p4 = fitter.fitted_p4();

       // save quantities after fit
       Ds_cand.addUserFloat("sv_chi2", fitter.chi2());
       Ds_cand.addUserFloat("sv_ndof", fitter.dof()); 
       Ds_cand.addUserFloat("sv_prob", fitter.prob());    
       Ds_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
       Ds_cand.addUserFloat("fitted_pt", fitter.fitted_candidate().globalMomentum().perp() );
       Ds_cand.addUserFloat("fitted_eta", fitter.fitted_candidate().globalMomentum().eta() );
       Ds_cand.addUserFloat("fitted_phi", fitter.fitted_candidate().globalMomentum().phi() );

       Ds_cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, Ds_cand.p4())
        );
       Ds_cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
       auto lxy = l_xy(fitter, *beamspot);
       Ds_cand.addUserFloat("l_xy", lxy.value());
       Ds_cand.addUserFloat("l_xy_unc", lxy.error());

       // second mass hypothesis
       auto fitted_trk1= fitter.daughter_p4(0);
       auto fitted_trk2= fitter.daughter_p4(1);
       auto fitted_trk3= fitter.daughter_p4(2);

       fitted_trk1.SetM(K_MASS);
       fitted_trk2.SetM(PI_MASS);
       fitted_trk3.SetM(K_MASS);

       Ds_cand.addUserFloat("fitted_barMass", (fitted_trk1+fitted_trk2+fitted_trk3).M() );
                                         
       // after fit selection
       //std::cout<<fitter.prob()<<std::endl;
       if( !post_vtx_selection_(Ds_cand) ) continue;
       //std::cout<<"got post selected candidate"<<std::endl;
       Ds_out->emplace_back(Ds_cand);
      }
    }
  }
  
  evt.put(std::move(Ds_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DsBuilder);

