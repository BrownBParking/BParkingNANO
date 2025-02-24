from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.setDefault('maxEvents', -1)
options.setDefault('tag', '10215')
options.parseArguments()

globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
    options.inputFiles = ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'] if not options.isMC else \
                         [
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD*.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD5000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD10000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD15000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD20000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD25000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD30000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD35000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD40000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD45000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD50000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD55000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD60000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD65000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD70000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD75000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD80000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD85000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD90000.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKPiPi_29Sep2020_MINIAOD95000.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_0.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_1.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_2.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_3.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_4.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_5.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_6.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_7.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_8.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_9.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_10.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_11.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_12.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_13.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_14.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_15.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToKDToKKstar0_19Oct2020_MINIAOD_16.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BsToRhoDs_ToKKPi_06Jul20.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/B0ToPiDToKPiPi_30Jul2020_MINIAOD.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BsToPiDsToPiPhi_10Aug2020_MINIAOD.root'
                         #'file:/eos/uscms/store/user/tkwon/BsToPiDs_ToKPiPi_MuFilter_TuneCP5_13TeV-pythia8-evtgen_03Nov20/MC_generation_BsToPiDs_ToKPiPi_MuFilter_TuneCP5_13TeV-pythia8-evtgen_MINIAOD_10Nov20/201110_145913/0000/BPH-RunIIAutumn18MiniAOD-00170_1.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpipi_27Sep2020_MINIAOD.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_0.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_1.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_2.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_3.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_4.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_5.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_6.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_7.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_8.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_9.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_10.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_11.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_12.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_13.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_14.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_15.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_16.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_17.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_18.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_19.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_20.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_21.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_22.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_23.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_24.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_25.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_26.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_27.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuTopiD0_10Oct2020_MINIAOD_28.root'
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD1.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD2.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD3.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD4.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD5.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD7.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD8.root',
                         #'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKpiRho_19Jun2020_MINIAOD9.root'
                         'file:/eos/uscms/store/user/tkwon/BParking/MC_Production/BuToKKpi_18May2020_MINIAOD.root'
                         ]
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')
# this is for the LowPt energy regression
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))







from PhysicsTools.BParkingNano.nanoBPark_cff import *
process = nanoAOD_customizeMuonTriggerBPark(process)
process = nanoAOD_customizeElectronFilteredBPark(process)
process = nanoAOD_customizeTrackFilteredBPark(process)
#process = nanoAOD_customizeBToKLL(process)
#process = nanoAOD_customizeBToKstarEE(process)
#process = nanoAOD_customizeBToKstarMuMu(process)
process = nanoAOD_customizeTriggerBitsBPark(process)
process = nanoAOD_customizeBToPiD0(process)
process = nanoAOD_customizeBsToPiDs(process)
process = nanoAOD_customizeB0ToKD(process)

# Path and EndPath definitions
#process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu )
#process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence + process.nanoBKeeSequence   + CountBToKee   )
#process.nanoAOD_KstarMuMu_step = cms.Path(process.nanoSequence + process.KstarToKPiSequence + process.nanoBKstarMuMuSequence + CountBToKstarMuMu )
#process.nanoAOD_KstarEE_step  = cms.Path(process.nanoSequence + process.KstarToKPiSequence + process.nanoBKstarEESequence + CountBToKstarEE  )
process.nanoAOD_step  = cms.Path(process.nanoSequence)
process.nanoAOD_BToPiD0_step  = cms.Path(process.nanoSequence + process.nanoBPiD0Sequence + CountBToPiD0  )
process.nanoAOD_BsToPiDs_step  = cms.Path(process.nanoSequence + process.nanoBsPiDsSequence + CountBsToPiDs  )
process.nanoAOD_B0ToKD_step  = cms.Path(process.nanoSequence + process.nanoB0KDSequence + CountB0ToKD  )


# customisation of the process.
if options.isMC:
   from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
   nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                #process.nanoAOD_KMuMu_step,
                                #process.nanoAOD_Kee_step, 
                                #process.nanoAOD_KstarMuMu_step,
                                #process.nanoAOD_KstarEE_step,
                                #process.nanoAOD_step,
																#process.nanoAOD_BToPiD0_step,
																#process.nanoAOD_BsToPiDs_step,
																process.nanoAOD_B0ToKD_step,
                                process.endjob_step, 
                                process.NANOAODoutput_step
                               )
if options.wantFullRECO:
    process.schedule = cms.Schedule(
                                    #process.nanoAOD_KMuMu_step,
                                    #process.nanoAOD_Kee_step, 
                                    #process.nanoAOD_KstarMuMu_step,
                                    #process.nanoAOD_KstarEE_step,
																		#process.nanoAOD_step,
																		#process.nanoAOD_BToPiD0_step,
																		#process.nanoAOD_BsToPiDs_step,
																		process.nanoAOD_B0ToKD_step,
                                    process.endjob_step, 
                                    process.FEVTDEBUGHLToutput_step, 
                                    process.NANOAODoutput_step
                                    )
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   #'nanoAOD_step',
                                   #'nanoAOD_BToPiD0_step'
                                   #'nanoAOD_BsToPiDs_step'
                                   'nanoAOD_B0ToKD_step'
                                   #'nanoAOD_KMuMu_step', 
                                   #'nanoAOD_Kee_step',
                                   #'nanoAOD_KstarMuMu_step',
                                   #'nanoAOD_KstarEE_step'
																	 ))


### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
