# Tool to convert AOD to NanoAOD

Tool to convert AOD to NanoAOD file format for the purpose of education and outreach

## Description

The tool can be used to read events from CMS AOD files and convert them to a reduced NanoAOD data format. Note that the tool is published for the documentation of the related datasets below and may need significant experiment-specific knowledge to be used.

## Setup CMSSW

In case CMSSW is set up outside of the [CMS Open Data VM](http://opendata.cern.ch/record/252) or of the [CMS Open Data docker container](http://opendata.cern.ch/docs/cms-guide-docker), source the following script.

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
git clone git://github.com/cms-opendata-analyses/AOD2NanoAODOutreachTool -b v1.2 AOD2NanoAOD
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

## Test workflows

This repository containes two github action workflows, which run the test workflow on the CMS open data container using github free resources.  

The workflow in [main.yml](.github/workflows/main.yml) runs a test job in a docker container. The run commands are passed in [commands.sh](commands.sh).

The workflow in [main_argo.yml](.github/workflow/main_argo.yml) sets up a minikube environment and runs a workflow defined with argo workflow engine. The workflow definition and run commands are in  [argo-workflow.yaml](argo-workflow.yaml).

The ouput is returned as a github artifact. The workflows are triggered by a pull request.
