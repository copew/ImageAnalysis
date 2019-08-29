#!/bin/bash
#
#PBS -N Matlab
#PBS -m be
#PBS -k oe


/usr/local/Cluster-Apps/matlab/R2017b/bin/matlab -nodesktop -nosplash -nodisplay <<EOF

addpath('/rds/user/ww234/hpc-work/ImageAnalysis')
tum_lymph_cluster_hpc('/rds/user/ww234/hpc-work/itpt/593987.fits')
;exit
EOF
