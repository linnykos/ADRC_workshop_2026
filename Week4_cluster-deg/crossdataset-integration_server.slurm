#!/bin/bash
#
# Week 4: Batch correction, part 2 — the SLURM job for the server half
# Kevin Z. Lin, 2026-08-20
#
# Submits crossdataset-integration_server.R as a batch job. From the folder
# holding both files, on the server:
#
#     sbatch crossdataset-integration_server.slurm
#
# That prints a job ID. Then:
#
#     squeue -u $USER                 # is it queued (PD) or running (R)?
#     scancel JOBID                   # stop it
#     tail -f crossdataset-integration_JOBID.out
#
# Every line below starting with #SBATCH is a directive to the scheduler, not a
# comment. They must all appear before the first real command.

# --- Things you MUST change ---------------------------------------------------

# Which allocation to bill and which set of machines to run on. There is no
# sensible default: these are specific to your cluster and your group, and the
# job will be rejected outright without them. On UW Hyak they look like
# --account=adrc and --partition=compute; on the SCRI HPC they are different
# again. Ask whoever gave you the account.
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION

# --- Resources ----------------------------------------------------------------

#SBATCH --job-name=crossdataset-integration
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# Seurat's integration is largely single-threaded unless you set up
# future::plan() yourself, which this script does not. A handful of cores helps
# the linear algebra underneath; asking for 16 buys you a longer wait in the
# queue and nothing else.
#SBATCH --cpus-per-task=4

# This is the directive most likely to need raising, and the failure is easy to
# misread. IntegrateData() holds several copies of a genes-by-cells matrix at
# once. If the job disappears with no R error in the .out file, look in the .err
# file for a scheduler line about the memory limit or an oom-kill event — that
# is the kernel stopping the job, not R, so there is no traceback to read.
#
# 64G is generous for the job as shipped, which integrates 10,000 nuclei against
# 9,354. The six-region variant of the same script (133,007 on the query side) is
# what wanted 128G. Raise this, never lower it, when you point the script at your
# own data.
#SBATCH --mem=64G

# Wall clock. The job is killed when this runs out, part-way, with nothing
# saved. Ask for more than you think you need; unused time is not charged. The
# job as shipped takes about 9 minutes, so this is roughly 13x headroom, which is
# the right way round to be wrong.
#SBATCH --time=2:00:00

# --- Logs ---------------------------------------------------------------------

# %j is the job ID, so resubmitting does not overwrite the previous logs.
# Everything the R script cat()s lands in .out; errors and warnings land in .err.
#SBATCH --output=crossdataset-integration_%j.out
#SBATCH --error=crossdataset-integration_%j.err

# Uncomment both to be emailed when the job ends or fails.
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=YOUR_NETID@uw.edu

# --- The job itself -----------------------------------------------------------

# How R is made available differs by cluster: a module, a container, or a conda
# environment. Uncomment whichever applies and check the version it gives you.
# module load R/4.4.1
# module load apptainer

Rscript crossdataset-integration_server.R
