#!/bin/bash

## IMPORTANT! Load local CCDB copy.  Make sure these are in your startup script.
source /etc/profile.d/modules.sh
module load sqlite/dev
module load xrootd



# Set job info
export MCINDEX=0 #NOTE: This will be set by `setup.sh`.
export PREFIX="out_rgh_NH3_" #NOTE: Make sure this matches what you have in `submit_clasdis.sh`.
export BASENAME="${PREFIX}${MCINDEX}"
export NEVENTS=100 #NOTE: This will be set by `setup.sh`.
export GCARD="$RGH_SIM_HOME/rgh_physics.gcard"
export YAML="$RGH_SIM_HOME/rgh_physics.yaml"
export OUTDIR="$RGH_SIM_VOL_DIR"
export CLAS12DIR="/work/clas12/users/mfmce/coatjava/coatjava/bin"

# Obsolete if using `submit_clasdis.sh`
# export EVPFILE=$NEVENTS
# export BEAM_ENERGY=10.6
# export TARGET_TYPE="proton"
# export XMIN=0.05 #NOTE: These are kinematic limits for the clasdis generator
# export XMAX=0.95

mkdir -p $OUTDIR

function check_task_status() {

    # Set arguments
    export TASK_NAME=$1
    export TASK_OUTPUT=$2
    export TASK_STATUS=$3
    export TASK_EXIT_CODE=$4

    # Print info
    echo "Task name:   $TASK_NAME"
    echo "Task output: $TASK_OUTPUT"
    echo "Task status: $TASK_STATUS"

    # Exit if status is nonzero
    if [ $TASK_STATUS -ne 0 ]; then
        echo "Task failed. Exit code $TASK_EXIT_CODE."
        exit $TASK_EXIT_CODE
    fi
    
}

# Run MC event Simulations (nh3 target 10.6 GeV)
export OUTDIR_LUND="${OUTDIR}/lund"
mkdir -p $OUTDIR_LUND
cd $OUTDIR_LUND #NOTE: Since clasdis does not like long input path names and just truncates them just cd to here and use basename
export LUNDFILE=$BASENAME #NOTE: This cannot be too long, otherwise clasdis will truncate it.
#NOTE: This step is run separately in `submit_clasdis.sh` now: clasdis --beam $BEAM_ENERGY --targ $TARGET_TYPE --trig $NEVENTS --nmax $EVPFILE --path $LUNDFILE --x $XMIN $XMAX
ls -lrth $OUTDIR_LUND/${BASENAME}*clasdis*.dat
export LUND_TASK_STATUS=$?
cd - #NOTE: cd back to wherever you were before clasdis
export LUNDFILE=`ls $OUTDIR_LUND/${BASENAME}clasdis*.dat`
check_task_status "clasdis" $LUNDFILE $LUND_TASK_STATUS 1

# Run GEMC detector simulation
export OUTDIR_GEMC="${OUTDIR}/gemc"
mkdir -p $OUTDIR_GEMC
export GEMCFILE=$OUTDIR_GEMC/$BASENAME.hipo
gemc $GCARD -SAVE_ALL_MOTHERS=1 -SKIPREJECTEDHITS=1 -NGENP=50 -INTEGRATEDRAW="*" -USE_GUI=0 -RUNNO=11 -INPUT_GEN_FILE="LUND, $LUNDFILE" -OUTPUT="hipo, $GEMCFILE" -N=$NEVENTS
export GEMC_TASK_STATUS=$?
check_task_status "gemc" $GEMCFILE $GEMC_TASK_STATUS 2

# Run CLAS12 reconstruction
export OUTDIR_REC="${OUTDIR}/cooked"
mkdir -p $OUTDIR_REC
export RECFILE=$OUTDIR_REC/$BASENAME.hipo
$CLAS12DIR/recon-util -i $GEMCFILE -o $RECFILE -n $NEVENTS -y $YAML
export REC_TASK_STATUS=$?
check_task_status "recon-util" $RECFILE $REC_TASK_STATUS 3

# Create DST
export OUTDIR_DST="${OUTDIR}/dst"
mkdir -p $OUTDIR_DST
export DSTFILE=$OUTDIR_DST/$BASENAME.hipo
export FILTER_BANKS='RUN::*,RAW::epics,RAW::scaler,HEL::flip,HEL::online,REC::*,RECFT::*,MC::RecMatch,MC::GenMatch,MC::Particle,MC::User,MC::Header,MC::Lund,MC::Event'
hipo-utils -filter -b $FILTER_BANKS -merge -o $DSTFILE $RECFILE
export DST_TASK_STATUS=$?
check_task_status "hipo-utils -filter" $DSTFILE $DST_TASK_STATUS 4

# Clean up intermediate files
#rm $LUNDFILE #NOTE: Keep these files because you need to know how much XS you generated and that is recorded in the file names.
rm $GEMCFILE
rm $RECFILE
