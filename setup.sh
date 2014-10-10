 #!/bin/bash
date

################################################
# Sunday, April 13, 2014: 12:48 hr             #                                                                                                                          
################################################

WORKDIR=/nfs/dust/cms/user/pooja/
TOOLDIR=$WORKDIR/scratch/plot-macro/
ANALYSISPATH=$WORKDIR/scratch/plot-macro/tau-hadronic/h2tautau-analysis/

# INPUT FILE
COMMONFILEPATH=$WORKDIR/samples/
FILEPATH=$COMMONFILEPATH/syncEx_tautau_Christian_Oct2014/output_crab/
FILENAME=$FILEPATH/TTFH_MC_SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola

# LUMI FILE
LUMIPATH=$TOOLDIR/Lumi/new
LUMIFILE=$LUMIPATH/LUMI_INFO_TTFH_MC_SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola

# CONFIG FILE
CONFIGPATH=$ANALYSISPATH/Config
CONFIGDEFAULT="YES"
CONFIGANALYSIS=$CONFIGPATH/analysis
CONFIGSKIM=$CONFIGPATH/skim
CONFIGDEFAULTFILE=$CONFIGANALYSIS/2012.cfg

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
GREP[12]="summer12_1jets.cfg"
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
GREP[23]="inclusive_fakerates.cfg"
GREP[24]="one_triplet.cfg"
GREP[25]="new_fakerates.cfg"
GREP[26]="highpt_selection_notrigger.cfg"
GREP[27]="summer12.cfg"
GREP[28]="same-sign-tau-pairs.cfg"

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
    cd $TOOLDIR/AnalysisTool_WH
    source AnalysisToolUseThis
    ./configure --prefix=$TOOLDIR/lib_AnalysisTool
    make
    make install
    echo '@@@@@@@@@@@@@ LIBERARIES ARE BUILD @@@@@@@@@@@@@@@@@@@@@@'
}

Skim() {
    make -B main_skim
    ./main_skim $LUMIFILE"*.root" $FILENAME"*.root" $CONFIGSKIM/$GREP[18]
}


Main() {
    make -B main

    ## detault config 
    if [ "$CONFIGDEFAULT" = "YES" ]; then
	echo '@@@@@@@@@@@@@@@@@@@    USING DEFAULT CONFIG ::' $CONFIGDEFAULTFILE
	./main $FILENAME"*root" $LUMIFILE"*.root"  $CONFIGDEFAULTFILE

    else
	## using all config files
	for ConfigFileIndx in "${GREP[@]}"
	  do
	  echo '@@@@@@@@@@@@@@@@@@@    USING FOLLOWING CONFIG ::' $ConfigFileIndx
           #echo "./main" $FILENAME"*root" $LUMIFILE"*.root" $CONFIGANALYSIS/$ConfigFileIndx
	  ./main $FILENAME"*root" $LUMIFILE"*.root"   $CONFIGANALYSIS/$ConfigFileIndx
	  
	  echo '@@@@@@@@@@@@@ Going to Sleep'
	  sleep 5
	  echo ""
	done
    fi
    echo 'EXISTING..'
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

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

