# Tool to convert AOD to NanoAOD

Tool to convert AOD to NanoAOD file format for the purpose of education and outreach

## Description

The tool can be used to read events from CMS AOD files and convert them to a reduced NanoAOD data format. Note that the tool is published for the documentation of the related datasets below and may need significant experiment-specific knowledge to be used.

## Setup CMSSW

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_5_3_32
```

## Build module

```bash
cd CMSSW_5_3_32/src
cmsenv
mkdir workspace
cd workspace
git clone <THIS REPOSITORY>
cd AOD2NanoAOD
scram b -j8
```

## Test configuration locally

```bash
cmsRun configs/simulation_cfg.py
cmsRun configs/data_cfg.py
```

## Create jobs for lxplus batch system

```bash
./submit_jobs.sh /path/to/job/directory
```

## Merge job files

```bash
./merge_jobs.py /path/to/job/outputs
```
