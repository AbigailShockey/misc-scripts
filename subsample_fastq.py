#!/usr/bin/env python 3

# This script subsamples gzipped fastq files down to a specified number of reads
# Requires seqtk to be installed and in your PATH

import glob
import subprocess as sub
import shlex
import sys

files = glob.glob("*.fastq.gz")

nums = [1000,2500,5000,7500,10000,25000,50000,75000,100000,250000,500000,750000,1000000]

for num in nums:
    for file in files:
        handle = file.split("_")[0]
        extension = file.split("_")[1].replace(".gz","")
        subsampleFile = f"{handle}-{num}_{extension}"
        # run seqtk and subsample reads
        subsampleFile = open(f"{subsampleFile}","w")
        cmd = shlex.split(f"seqtk sample -s100 {file} {num}")
        sub.Popen(cmd, stdout=subsampleFile).wait()

files = glob.glob("*.fastq")
files = ' '.join(files)
cmd = shlex.split(f"gzip {files}")
sub.Popen(cmd).wait()
