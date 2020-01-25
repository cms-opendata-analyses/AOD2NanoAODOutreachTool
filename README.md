# Tool to convert AOD to NanoAOD

Tool to convert AOD to NanoAOD file format for the purpose of education and outreach

## Description

The tool can be used to read events from CMS AOD files and convert them to a reduced NanoAOD data format. Note that the tool is published for the documentation of the related datasets below and may need significant experiment-specific knowledge to be used.

## Setup CMSSW

KLP: In case running on lxplus, set first the slc6 environment (see [instructions](http://cms-sw.github.io/singularity.html))

```bash
cmssw-slc6
```

In case CMSSW is set up outside of the [CMS Open Data VM](http://opendata.cern.ch/record/252), source the following script.

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

And check-out the appropriate CMSSW release using following call.

```bash
cmsrel CMSSW_5_3_32
```

## Build module

```bash
cd CMSSW_5_3_32/src
cmsenv
mkdir workspace
cd workspace
git clone git://github.com/cms-opendata-analyses/AOD2NanoAODOutreachTool -b 2012 AOD2NanoAOD
cd AOD2NanoAOD
scram b -j8
```

## Test configuration locally

```bash
cmsRun configs/simulation_cfg.py
cmsRun configs/data_cfg.py
```

## Example scripts for batch system submission

You can use the following script to submit to any [HTCondor](https://research.cs.wisc.edu/htcondor/) batch system.

```bash
./submit_jobs.sh /path/to/job/directory
```

You can merge the job files with the following script.

```bash
./merge_jobs.py /path/to/job/outputs
```
