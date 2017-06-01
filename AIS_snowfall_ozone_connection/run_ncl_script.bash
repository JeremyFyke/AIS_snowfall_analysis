#!/bin/tcsh
#
# LSF batch script to run an NCL script
#
#BSUB -P P93300301           # project code
#BSUB -W 0:30               # wall-clock time (hrs:mins)
#BSUB -n 1                   # number of tasks in job         
#BSUB -J ozone            # job name
#BSUB -R "rusage[mem=512000]"     # memory to reserve, in MB
#BSUB -o myjob.%J.out        # output file name in which %J is replaced by the job ID
#BSUB -e myjob.%J.err        # error file name in which %J is replaced by the job ID
#BSUB -q geyser              # queue


source /glade/u/apps/opt/lmod/4.2.1/init/tcsh
module load ncl
ncl plot_oz_sf_jan.ncl

