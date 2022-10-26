#!/usr/bin/env python 3

# This script 'contaminates' a list of gzipped fastqs with a percentage of reads from gzipped fastqs from a list
# Requires seqtk to be installed and in your PATH, as well as the pyfastx python package

import glob
import pyfastx
import subprocess as sub
import shlex
import sys

fastqs = sys.argv[1]
contams = sys.argv[2]

fastq_files = []
with open(fastqs,"r") as fastqFile:
    for fastq in fastqFile:
        fastq = fastq.strip()
        fastq_files.append(fastq)

contam_files = []
with open(contams,"r") as contamFile:
    for contam in contamFile:
        contam = contam.strip()
        contam_files.append(contam)

# remove contamination file from list of reads
files = list(set(fastq_files)-set(contam_files))

nums = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

fwd = []
rev = []
contam_fwd = []
contam_rev = []

for num in nums:
    for file in files:
        # load fastq file
        fq = pyfastx.Fastq(file)
        # number of reads in fastq
        num_reads = len(fq)
        # number of reads to keep
        prop_reads = int(num_reads * num)
        # number of reads to toss
        num_rm = num_reads - prop_reads
        # percent of reads kept
        perc = str(int(num * 100)) + "percent"
        # percent of reads tossed
        perc_rm = str((100-int(num * 100))) + "percent"
        # set up subsample output file
        handle = file.split("_")[0]
        extension = file.split("_")[1].replace(".gz","")
        subsampleFile = f"{handle}-{perc}_{extension}"
        if "R1" in extension:
            fwd.append(subsampleFile)
        if "R2" in extension:
            rev.append(subsampleFile)
        # run seqtk and subsample reads
        subsampleFile = open(f"{subsampleFile}","w")
        cmd = shlex.split(f"seqtk sample -s100 {file} {prop_reads}")
        sub.Popen(cmd, stdout=subsampleFile).wait()
        for contam_file in contamination_files:
            contam_handle = contam_file.split("_")[0]
            contam_extension = contam_file.split("_")[1].replace(".gz","")
            contam_subsampleFile = f"{contam_handle}-{perc_rm}_{contam_extension}"
            if "R1" in contam_extension:
                contam_fwd.append(contam_subsampleFile)
            if "R2" in contam_extension:
                contam_rev.append(contam_subsampleFile)
            contam_subsampleFile = open(f"{contam_subsampleFile}","w")
            cmd = shlex.split(f"seqtk sample -s100 {contam_file} {num_rm}")
            sub.Popen(cmd, stdout=contam_subsampleFile).wait()

fwd = list(set(fwd))
rev = list(set(rev))
contam_fwd = list(set(contam_fwd))
contam_rev = list(set(contam_rev))

for reads in fwd:
    for contam_reads in contam_fwd:
        reads_handle = reads.split("_")[0]
        contam_reads_handle = contam_reads.split("_")[0]
        combinedFile = open(f"{reads_handle}_{contam_reads_handle}_R1.fastq","w")
        cmd = shlex.split(f"cat {reads} {contam_reads}")
        sub.Popen(cmd, stdout=combinedFile).wait()

for reads in rev:
    for contam_reads in contam_rev:
        reads_handle = reads.split("_")[0]
        contam_reads_handle = contam_reads.split("_")[0]
        combinedFile = open(f"{reads_handle}_{contam_reads_handle}_R2.fastq","w")
        cmd = shlex.split(f"cat {reads} {contam_reads}")
        sub.Popen(cmd, stdout=combinedFile).wait()

files = glob.glob("*.fastq")
files = ' '.join(files)
cmd = shlex.split(f"gzip {files}")
sub.Popen(cmd).wait()
