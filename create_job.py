#!/usr/bin/env python


import os
import sys


jdl = """\
universe = vanilla
executable = ./{PROCESS}.sh
output = out/$(ProcId).$(ClusterID).out
error = err/$(ProcId).$(ClusterID).err
log = log/$(ProcId).$(ClusterID).log
requirements = (OpSysAndVer =?= "SLCern6")
max_retries = 3
RequestCpus = 1
+MaxRuntime = 1800
queue arguments from arguments.txt\
"""


def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


def parse_arguments():
    if not len(sys.argv) == 3:
        raise Exception("./create_job.py PROCESS PATH_TO_JOBDIR")
    return {"process": sys.argv[1], "jobdir": sys.argv[2]}


def main(args):
    process = args["process"]
    print("Process: %s" % process)

    # Build argument list
    print("Filelist:")
    arguments = []
    counter = 0
    for filename in os.listdir("data/"):
        if process in filename:
            print("    %s." % filename)
            for line in open("data/" + filename, "r").readlines():
                arguments.append("%u %s %s" % (counter, process, line))
                counter += 1
    print("Number of jobs: %u" % len(arguments))

    # Create jobdir and subdirectories
    jobdir = os.path.join(args["jobdir"], process)
    print("Jobdir: %s" % jobdir)
    mkdir(jobdir)
    mkdir(os.path.join(jobdir, "out"))
    mkdir(os.path.join(jobdir, "log"))
    mkdir(os.path.join(jobdir, "err"))

    # Write jdl file
    out = open(os.path.join(jobdir, "job.jdl"), "w")
    out.write(jdl.format(PROCESS=process))
    out.close()

    # Write argument list
    arglist = open(os.path.join(jobdir, "arguments.txt"), "w")
    for a in arguments:
        arglist.write(a)
    arglist.close()

    # Write job file
    jobfile = open("job.sh", "r").read()
    job = open(os.path.join(jobdir, "{PROCESS}.sh".format(PROCESS=process)), "w")
    job.write(jobfile)
    job.close()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
