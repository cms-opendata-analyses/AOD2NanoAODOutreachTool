#!/bin/sh -l

# echo pwd `pwd`: /home/cmsusr/CMSSW_5_3_32/src
# echo $USER cmsusr
sudo chown $USER /mnt/vol

mkdir workspace
cd workspace
# Clone done here for the compatibility with the minikube workflow 
# For the plain github action with docker, the are would be available in /mnt/vol
git clone git://github.com/katilp/AOD2NanoAODOutreachTool  AOD2NanoAOD
cd AOD2NanoAOD
scram b -j8

#cmsRun configs/simulation_cfg.py
cmsRun configs/data_cfg.py

cp output.root /mnt/vol/ 
echo  ls -l /mnt/vol 
ls -l /mnt/vol



