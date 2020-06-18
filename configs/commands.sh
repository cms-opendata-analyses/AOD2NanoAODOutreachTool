#!/bin/sh -l

mkdir workspace
cd workspace
git clone git://github.com/katilp/AOD2NanoAODOutreachTool -b v1.2 AOD2NanoAOD
cd AOD2NanoAOD
scram b -j8

#cmsRun configs/simulation_cfg.py
cmsRun configs/data_cfg.py

ls -l 

sudo chown -R cmsusr /mountedvolume
chmod 755 /mountedvolume
mkdir /mountedvolume/outputs
#cp *.pdf /mountedvolume/outputs
ls -l /mountedvolume/outputs
touch /mountedvolume/outputs/empty.out