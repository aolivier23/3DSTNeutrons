#!/bin/bash

#Usage example with jobsub:
#
#jobsub_submit -N 10 -M --OS=SL6 --group=dune --resource-provides=usage_model=OPPORTUNISTIC file:///dune/app/users/aolivier/nd_neutrons/EDep_grid.sh 
#dune/persistent/users/aolivier/edep-sim/neutrons '(.*)-edep_'${PROCESS}'\.root' /pnfs/dune/persistent/users/aolivier/neutrons/PosRes/default/hits 
#hits_with_meta /pnfs/dune/persistent/users/aolivier/tardir --reco GridNeutronHits --ana FSNeutrons --E-min 2.0
#
#Variables that I need for this job.  If I want to generate this script, I just need to change these directories
INDIR=root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/$0
echo "INDIR is ${INDIR}\n"
REGEX= $1 #'(.*)-edep_'${PROCESS}'\.root'
echo "REGEX is ${REGEX}\n"
#TODO: Couldn't I stream my output files as well?  
OUTDIR=$2
echo "Sending output to OUTDIR: ${OUTDIR}\n"
OUTBASE=$3_File${PROCESS}
echo "Output files start with ${OUTBASE}\n"
#TODO: Wrapper script that figures out username and guesses pnfs directories appropriately
TARDIR=/pnfs/$4
echo "Getting tar files from ${TARDIR}\n"
shift 5 #eat arguments used for required parameters

#Set up the Fermilab products I need
#Copied from Mike Kordosky's working area on 10/24/2017 at his suggestion.
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

#### one version #################
setup nutools v2_15_00 -q e14:debug
setup cmake v3_9_0
setup dk2nu v01_05_00e -q e14:debug
setup genie_xsec   v2_12_6 -q DefaultPlusMECWithNC
setup genie_phyopt v2_12_6  -q dkcharmtau

# edep-sim needs to know where a certain GEANT .cmake file is...
export PATH=/cvmfs/fermilab.opensciencegrid.org/products/larsoft/geant4/v4_10_3_p01a/Linux64bit+2.6-2.12-e14-debug/bin/:$PATH
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

export GXMLPATH=$GXMLPATH:/dune/app/users/kordosky/nd_sim/event_gen

#Taken from Jose Palomino's edepsim grid job
#TODO: Is the following "the section you should not change"?  Not messing with it for now, but one day I want to understand it.
setup jobsub_client

TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0
cd $TMP
# End of the section you should not change.

echo "Scratch directory: $TMP"

#Note: Using xrootd for input files, so I shouldn't have to copy them.

#tar up my edepsim directory.  Should be able to get it from the EDEP environment variables it sets up. 
ifdh cp -D ${TARDIR}/edepsim.tar.gz . 
tar -zxvf edepsim.tar.gz 
cd edep-sim 
source setup.sh 
cd ../ 
echo "Set up edepsim successfully.\n"

#untar EDepNeutrons and util
#TODO: Get names of Utilities and EdepNeutrons tar files automatically from CMake
ifdh cp -D ${TARDIR}/Utilities-1.0.0.tar.gz .
tar -zxvf Utilities.tar.gz
echo "Set up Utilities successfully.\n"

ifdh cp -D ${TARDIR}/EdepNeutrons-2.0.0.tar.gz .
tar -zxvf EdepNeutrons.tar.gz
source EdepNeutrons/bin/setup.sh
echo "Set up EdepNeutrons successfully.\n"

#My EDepApp command
#Forward all command line arguments to EdepApp
/usr/bin/time -v ./EdepNeutrons/bin/EdepApp --path ${INDIR} --regex ${REGEX} --reco-file ${OUTBASE}.root $@ &> ${OUTBASE}_log.txt 

#Copy results back to pnfs
ifdh cp -D ${OUTBASE}.root ${OUTDIR}
ifdh cp -D ${OUTBASE}_log.txt ${OUTDIR}
ifdh cp -D histos_${OUTBASE}.root ${OUTDIR}
ifdh cp -D JobThatMade_${OUTBASE}.txt ${OUTDIR}
