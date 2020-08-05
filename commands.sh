#!/bin/sh -l
# parameters: $1 number of events, $2 configuration file
# echo pwd `pwd`: /home/cmsusr/CMSSW_5_3_32/src
# echo $USER cmsusr
sudo chown $USER /mnt/vol

mkdir workspace
cd workspace
# For the plain github action with docker, the area would be available in /mnt/vol
git clone git://github.com/cms-opendata-analyses/AOD2NanoAODOutreachTool  AOD2NanoAOD
cd AOD2NanoAOD
scram b -j8

if [ -z "$1" ]; then nev=10000; else nev=$1; fi
if [ -z "$2" ]; then config=data_cfg.py; else config=$2; fi
eventline=$(grep maxEvents configs/$config)
sed -i "s/$eventline/process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32($nev) )/g" configs/$config
cmsRun configs/$config

cp output.root /mnt/vol/
echo  ls -l /mnt/vol
ls -l /mnt/vol



