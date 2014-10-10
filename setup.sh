#!/bin/bash
date

################################################
# Sunday, April 13, 2014: 12:48 hr             #                                                                                                                          
################################################

WORKDIR=/nfs/dust/cms/user/pooja/
ANALYSISPATH=$WORKDIR/scratch/plot-macro/tau-hadronic/h2tautau-analysis/

#INPUT FILE
COMMONFILEPATH=$WORKDIR/samples/
FILEPATH=$COMMONFILEPATH/syncEx_tautau_Christian_Oct2014/output_crab/
FILENAME=$FILEPATH/TTFH_MC_SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola

# LUMI FILE
LUMIPATH=$ANALYSISPATH/Lumi/new
LUMIFILE=LUMI_INFO_TTFH_MC_SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola

# CONFIG FILE
CONFIGPATH=$ANALYSISPATH/Config
CONFIGTYPE=$CONFIGPATH/analysis
CONFIGFILE=$CONFIGTYPE/2012.cfg

GREP[1]="2011.cfg"
GREP[2]="matching.cfg"
GREP[3]="notrigger.cfg"
GREP[4]="summer12_1jet.cfg"
GREP[5]="very_loose_ostau_isolation.cfg"
GREP[6]="2012.cfg"
GREP[7]="electron_iso.cfg"
GREP[8]="fall11_4jets.cfg"               
GREP[9]="highpt_selection.cfg"
GREP[10]="highpt_selection_notrigger_matching.cfg"
GREP[11]="new_fakerates.cfg"
GREP[12]="same-sign-tau-pairs.cfg"
GREP[13]="summer12_2jets.cfg"
GREP[14]="summer12_3jets.cfg"
GREP[15]="summer12_4jets.cfg"
GREP[16]="tau_down.cfg"
GREP[17]="tau_up.cfg"
GREP[18]="common.cfg"
GREP[19]="conversion_rejection.cfg"
GREP[20]="all_triplets.cfg"
GREP[21]="capped_fakerates.cfg"
GREP[22]="isoonly_fr.cfg"
GREP[23]="matching.cfg"


usage () {
    echo "@@@@@@ Higgs_tau_tau "
    echo "              ==Possible arguments=="
    echo "  makelib                                 - to build lib of MainAnalyzer"
    echo "  skim                                    - to execute skim.cc"
    echo "  main                                    - to execute main.cc"
    echo "  --help                                  - shows this help"
}


SetupEnv() {
#    source $SETUP
    source $LIB
}

MakeLib() {
    cd /nfs/dust/cms/user/pooja/scratch/plot-macro/AnalysisTool_WH
    source AnalysisToolUseThis
    ./configure --prefix=/nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool
    make
    make install
    echo '@@@@@@@@@@@@@ LIBERARIES ARE BUILD @@@@@@@@@@@@@@@@@@@@@@'
}

Skim() {
    make -B main_skim
# test prupose (single sample)
# ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root  /nfs/dust/cms/user/pooja/scratch/plot-macro/armin-file/TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6_*.root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg

# test purpose (huge dataset)
#    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root /nfs/dust/cms/user/pooja/samples/h2tautau-8TeV/TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6_28*root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
#    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV*root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6_4.root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
#    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root  TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6_29.root  /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
}


Main() {
    make -B main
#    ./main skimmed.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main finalskim.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_patch2/src/MyRootMaker/MyRootMaker/test/GluGluToHToTauTau_M-125.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_H2Tau/src/MyRootMaker/MyRootMaker/test/VBF_HToTauTau_M-125.root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6_10.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_patch2/src/MyRootMaker/MyRootMaker/test/GluGluToHToTauTau_M-125_full.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/VBF_HToTauTau_M-125_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/GluGluToHToTauTau_M-125_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/SUSYGluGluToHToTauTau_M-120_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg  
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/DYJetsToLL_M-50_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_DYJetsToLL_M-50_TuneZ2star_8TeV_madgraph.root ../Config/analysis/2012.cfg
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/VBF_HToTauTau_M-125_tauPolarOff_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/GluGluToHToTauTau_M-125_tauPolarOff_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_14/src/MyRootMaker/MyRootMaker/GluGluToHToTauTau_M-125_full.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg
#    ./main /nfs/dust/cms/user/pooja/samples/cpStudy_Aug2014/GluGluToHToTauTau_M-125_MC_v8_vtxWithBS_test/GluGluToHToTauTau_M-125*.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 

#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6_4.root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg
    ./main $FILEPATH/TTFH_MC_SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola*root   ../Config/analysis/2012.cfg  

}

#if [ $# -lt 1 ]; then
#    usage
#    exit 1
#fi

if [ "$1" = "--help" ]; then
    usage

elif [ "$1" = "makelib" ]; then
    SetupEnv
    MakeLib

#main_skim
elif [ "$1" = "skim" ]; then
    SetupEnv
    Skim

    
#main
elif [ "$1" = "main" ]; then
    SetupEnv
    Main
fi

