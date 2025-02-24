from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.BParkingNano.trgbits_cff import *



##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.tracksBPark_cff import *

## B collections
#from PhysicsTools.BParkingNano.BToKLL_cff import *
#from PhysicsTools.BParkingNano.BToKstarLL_cff import *
from PhysicsTools.BParkingNano.BToPiD0_cff import *
from PhysicsTools.BParkingNano.BsToPiDs_cff import *
from PhysicsTools.BParkingNano.B0ToKD_cff import *

nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +           
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              globalTablesMC + genWeightsTable + genParticleBParkTables + lheInfoTable) 



def nanoAOD_customizeMuonTriggerBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoBKeeSequence     = cms.Sequence( electronsBParkSequence + electronBParkTables)
    process.nanoBKstarEESequence = cms.Sequence( electronsBParkSequence + electronBParkTables)
    return process

def nanoAOD_customizeTriggerBitsBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

#def nanoAOD_customizeBToKLL(process):
#    process.nanoBKeeSequence   = cms.Sequence( process.nanoBKeeSequence + BToKEESequence    + BToKeeTable   )
#    process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
#    return process

#three possibilities for K*LL
#def nanoAOD_customizeBToKstarLL(process):
#    process.nanoBKstarLLSequence   = cms.Sequence( KstarToKPiSequence + BToKstarLLSequence + KstarToKPiTable + BToKstarLLTables )
#    return process

#def nanoAOD_customizeBToKstarEE(process):
#    process.nanoBKstarEESequence   = cms.Sequence( process.nanoBKstarEESequence + BToKstarEESequence + BToKstarEETable + KstarToKPiTable )
#    return process

#def nanoAOD_customizeBToKstarMuMu(process):
#    process.nanoBKstarMuMuSequence = cms.Sequence( BToKstarMuMuSequence + BToKstarMuMuTable + KstarToKPiTable )
#    return process

#BToPiD0
def nanoAOD_customizeBToPiD0(process):
    process.nanoBPiD0Sequence = cms.Sequence( BToPiD0Sequence + BToPiD0Table + D0ToKPiTable )
    return process
#BsToPiDs
def nanoAOD_customizeBsToPiDs(process):
    process.nanoBsPiDsSequence = cms.Sequence( BsToPiDsSequence + BsToPiDsTable + DsToKKPiTable )
    return process

#B0ToKD
def nanoAOD_customizeB0ToKD(process):
    process.nanoB0KDSequence = cms.Sequence( B0ToKDSequence + B0ToKDTable + DToKPiPiTable)
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'electronsForAnalysis:SelectedElectrons', 'selectedElectronsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonBParkSequence, process.muonBParkMC)
        path.replace(process.electronsBParkSequence, process.electronBParkMC)
        path.replace(process.tracksBParkSequence, process.tracksBParkMC)
