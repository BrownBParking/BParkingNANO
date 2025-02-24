import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

DToKPiPi = cms.EDProducer(
       'DBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization   
        trk2Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization
        trk3Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization 
        beamSpot = cms.InputTag("offlineBeamSpot"),
        preVtxSelection = cms.string( 
#        ''
#        ' pt()> 3.0 '
        ' ( (mass() > 1.82 && mass() < 1.9) )'
#        ' ( (mass() > 1.56 && mass() < 2.16) )'
        ),
        postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.01'
        '&& (  (userFloat("fitted_mass") > 1.82 && userFloat("fitted_mass") < 1.9)  )'
#        '&& (  (userFloat("fitted_mass") > 1.56 && userFloat("fitted_mass") < 2.16)  )'
        '&& userFloat("l_xy")/userFloat("l_xy_unc") > 5. '
)
)
B0ToKD = cms.EDProducer(
    'B0ToKDBuilder',
    D = cms.InputTag('DToKPiPi'),
    DTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    trgMuon    = cms.InputTag("muonTrgSelector:trgMuons"),
    kaonSelection = cms.string('pt > 1.0 && abs(eta)<2.4'),
    isoTracksSelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection = cms.string(
#        'pt > 5.'
        ' ( mass > 4.0 && mass < 7.0)'
#        '((mass > 5.18 && mass < 5.23) | (mass > 5.33 && mass < 5.38))'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.'
#        '&& userFloat("sv_dRtrgMu") > 0.4'
        '&& userFloat("fitted_cos_theta_2D") >= 0.999 '
        '&& (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.0)'
				'&& userFloat("l_xy")/userFloat("l_xy_unc") > 7. '
#        '&& ((userFloat("fitted_mass") > 5.18 && userFloat("fitted_mass") < 5.23) | (userFloat("fitted_mass") > 5.33 && userFloat("fitted_mass") < 5.38))'
    )
)

################################### Tables #####################################
KstarToKPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("KstarToKPi"),
    cut = cms.string(""),
    name = cms.string("Kstar"),
    doc = cms.string("Kstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      barMass = ufloat('barMass'),
      fitted_mass = ufloat('fitted_mass'),
      fitted_barMass = ufloat('fitted_barMass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      cos2D = ufloat('cos_theta_2D'),
      fitted_cos2D = ufloat('fitted_cos_theta_2D'),
      l_xy = ufloat('l_xy'),
      l_xy_unc = ufloat('l_xy_unc'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')
    )
)


DToKstarPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("DToKstarPi"),
    cut = cms.string(""),
    name = cms.string("D"),
    doc = cms.string("D Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      fitted_mass = ufloat('fitted_mass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      cos2D = ufloat('cos_theta_2D'),
      fitted_cos2D = ufloat('fitted_cos_theta_2D'),
      l_xy = ufloat('l_xy'),
      l_xy_unc = ufloat('l_xy_unc'),
      svprob = ufloat('sv_prob'),         
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx'),
      pi_idx = uint('pi_idx')
    )
)
DToKPiPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("DToKPiPi"),
    cut = cms.string(""),
    name = cms.string("D"),
    doc = cms.string("D Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      barMass = ufloat('barMass'),
      fitted_mass = ufloat('fitted_mass'),
      fitted_barMass = ufloat('fitted_barMass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      cos2D = ufloat('cos_theta_2D'),
      fitted_cos2D = ufloat('fitted_cos_theta_2D'),
      l_xy = ufloat('l_xy'),
      l_xy_unc = ufloat('l_xy_unc'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx'),
      trk3_idx = uint('trk3_idx')
    )
)
B0ToKDTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("B0ToKD"),
    cut = cms.string(""),
    name = cms.string("B0ToKD"),
    doc = cms.string("B0ToKD Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trk1Idx = uint('trk1_idx'),
        trk2Idx = uint('trk2_idx'),
        trk3Idx = uint('trk3_idx'),
        KIdx = uint('K_idx'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        dRtrgMu = ufloat('sv_dRtrgMu'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # D0 fitted in b0 vertex
        #fit_D0_mass = ufloat('fitted_D0_mass'),
        #fit_D0_pt = ufloat('fitted_D0_pt'),
        #fit_D0_eta = ufloat('fitted_D0_eta'),
        #fit_D0_phi = ufloat('fitted_D0_phi'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_trk1_pt = ufloat('fitted_trk1_pt'),
        fit_trk1_eta = ufloat('fitted_trk1_eta'),
        fit_trk1_phi = ufloat('fitted_trk1_phi'),
        fit_trk2_pt = ufloat('fitted_trk2_pt'),
        fit_trk2_eta = ufloat('fitted_trk2_eta'),
        fit_trk2_phi = ufloat('fitted_trk2_phi'),
        fit_K_pt = ufloat('fitted_K_pt'),
        fit_K_eta = ufloat('fitted_K_eta'),
        fit_K_phi = ufloat('fitted_K_phi'),
        fit_trk3_pt = ufloat('fitted_trk3_pt'),
        fit_trk3_eta = ufloat('fitted_trk3_eta'),
        fit_trk3_phi = ufloat('fitted_trk3_phi'),
        trk1_iso03 = ufloat('trk1_iso03'),
        trk1_iso04 = ufloat('trk1_iso04'),
        trk2_iso03 = ufloat('trk2_iso03'),
        trk2_iso04 = ufloat('trk2_iso04'),
        K_iso03 = ufloat('K_iso03'),
        K_iso04 = ufloat('K_iso04'),
        trk3_iso03  = ufloat('trk3_iso03'),
        trk3_iso04  = ufloat('trk3_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        n_trk3_used = uint('n_trk3_used'),
        n_trk1_used = uint('n_trk1_used'),
        n_trk2_used = uint('n_trk2_used'),
        n_K_used = uint('n_K_used'),

    )
)

CountB0ToKD = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("B0ToKD")
)    

########################### Sequencies  ############################

B0ToKDSequence = cms.Sequence(
    (DToKPiPi + B0ToKD)
)

