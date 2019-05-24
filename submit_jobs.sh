#!/bin/bash

# Define path for job directories
BASE_PATH=/path/to/job/directory
mkdir -p $BASE_PATH

# Set processes
PROCESSES=( \
    #SMHiggsToZZTo4L \
    #ZZTo2e2mu \
    #ZZTo4mu \
    #ZZTo4e \
    #GluGluToHToTauTau \
    #VBF_HToTauTau \
    #TTbar \
    #W1JetsToLNu \
    #W2JetsToLNu \
    #W3JetsToLNu \
    #DYJetsToLL \
    #Run2012B_TauPlusX\
    #Run2012C_TauPlusX\
    #DY2JetsToLL \
    #DY3JetsToLL \
    #DY4JetsToLL \
    #Run2012B_SingleMu\
    #Run2012C_SingleMu\
    #Run2012B_DoubleMuParked \
    #Run2012C_DoubleMuParked \
    #Run2012B_DoubleElectron \
    #Run2012C_DoubleElectron \
    )

# Create JDL files and job directories
for PROCESS in ${PROCESSES[@]}
do
    python create_job.py $PROCESS $BASE_PATH
done

# Submit jobs
THIS_PWD=$PWD
for PROCESS in ${PROCESSES[@]}
do
    cd $BASE_PATH/$PROCESS
    condor_submit job.jdl
    cd $THIS_PWD
done
