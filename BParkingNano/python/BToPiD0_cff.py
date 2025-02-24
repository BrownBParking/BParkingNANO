import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

D0ToKPi = cms.EDProducer(
       'D0Builder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 1. && abs(eta)<2.4'), #need optimization   
        trk2Selection = cms.string('pt > 1. && abs(eta)<2.4'), #need optimization
        beamSpot = cms.InputTag("offlineBeamSpot"),
        preVtxSelection = cms.string( 
#        ''
#        ' pt()> 2.0 '
        '( (mass() > 1.82 && mass() < 1.9) )'
#        '( (mass() > 1.56 && mass() < 2.16) )'
        ),
        postVtxSelection = cms.string(
        'userFloat("sv_prob") > 1.e-2'
        '&&(  (userFloat("fitted_mass")>1.82 && userFloat("fitted_mass")<1.9) )'
#        '&&(  (userFloat("fitted_mass")>1.56 && userFloat("fitted_mass")<2.16) )'
				'&& userFloat("l_xy")/userFloat("l_xy_unc") > 5. '
)
)

BToPiD0 = cms.EDProducer(
    'BToPiD0Builder',
    D0s = cms.InputTag('D0ToKPi'),
    D0sTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    pions = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    trgMuon    = cms.InputTag("muonTrgSelector:trgMuons"),
    piSelection = cms.string('pt > 1. && abs(eta)<2.4'),
    isoTracksSelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection = cms.string(
#        'pt > 5.'# && userFloat("min_dr") > 0.4 '
        ' (( mass > 4. && mass < 7.))'
#        ' ((mass > 4.93 && mass < 5.1) | (mass > 5.44 && mass < 5.61))'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.0 '
#        '&& userFloat("sv_dRtrgMu") > 0.4'
        '&& userFloat("fitted_cos_theta_2D") > 0.999 '
        '&& (( userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.))'
#        '&&((userFloat("fitted_mass") > 4.93 && userFloat("fitted_mass") < 5.1)|( userFloat("fitted_mass") > 5.44 && userFloat("fitted_mass") < 5.61))'
				'&& userFloat("l_xy")/userFloat("l_xy_unc") > 7. '
    )
)

################################### Tables #####################################

D0ToKPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("D0ToKPi"),
    cut = cms.string(""),
    name = cms.string("D0"),
    doc = cms.string("D0 Variables"),
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

BToPiD0Table = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToPiD0"),
    cut = cms.string(""),
    name = cms.string("BToPiD0"),
    doc = cms.string("BToPiD0 Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trk1Idx = uint('trk1_idx'),
        trk2Idx = uint('trk2_idx'),
        piIdx = uint('pi_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
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
        # D0 fitted in b+ vertex
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
        fit_pi_pt = ufloat('fitted_pi_pt'),
        fit_pi_eta = ufloat('fitted_pi_eta'),
        fit_pi_phi = ufloat('fitted_pi_phi'),
        # additional mass hypothesis
        barMass = ufloat ('barMass'),
        fit_barMass = ufloat('fitted_barMass'),
        fit_barD0_mass = ufloat('barMassD0_fullfit'),
        trk1_iso03 = ufloat('trk1_iso03'),
        trk1_iso04 = ufloat('trk1_iso04'),
        trk2_iso03 = ufloat('trk2_iso03'),
        trk2_iso04 = ufloat('trk2_iso04'),
        pi_iso03  = ufloat('pi_iso03'),
        pi_iso04  = ufloat('pi_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        n_pi_used = uint('n_pi_used'),
        n_trk1_used = uint('n_trk1_used'),
        n_trk2_used = uint('n_trk2_used'),
    )
)

CountBToPiD0 = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToPiD0")
)    

########################### Sequencies  ############################

BToPiD0Sequence = cms.Sequence(
    (D0ToKPi * BToPiD0)
)

