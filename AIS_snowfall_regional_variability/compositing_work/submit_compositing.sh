#!/bin/bash

#BSUB -P P93300601
#BSUB -n 1
#BSUB -W 12:00
#BSUB -q geyser
#BSUB -J generate_composites 
#BSUB -o generate_composites.%J.out       
#BSUB -e generate_composites.%J.err      
source /glade/u/apps/opt/lmod/4.2.1/init/bash
module load ncl
ncl generate_composites.ncl
