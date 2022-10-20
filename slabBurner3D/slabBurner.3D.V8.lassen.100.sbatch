#!/bin/tcsh
# see https://hpc.llnl.gov/documentation/tutorials/using-lc-s-sierra-systems#BatchScripts

### LSF syntax
#BSUB -nnodes 1                   #number of nodes
#BSUB -W 5                      #walltime in minutes
#BSUB -G sunyb                   #account
#BSUB -e error.%J.txt             #stderr
#BSUB -o out.%J.txt             #stdout
#BSUB -J test                    #name of job
#BSUB -q pdebug                   #queue to use

### Shell scripting
date; hostname
echo -n 'JobID is '; echo $LSB_JOBID
cd /p/gpfs1/mcgurn4/


# not -n should be nnodes*2
### Launch parallel executable
jsrun -n2 -r2 -a20 -c20 ./a.out