#!/bin/bash

myscriptdir=/rds/user/ww234/hpc-work/ImageAnalysis
myfilesdir=/rds/user/ww234/hpc-work/itpt

submit=${myscriptdir}/cellcluster_submit.sh
prepare=${myscriptdir}/cellcluster_prepare.sh
func=${myscriptdir}/tum_lymph_cluster_hpc_diff_parameters.m
#!images=${myfilesdir}/repeatfiles.txt #! the repeats
#!images=${myfilesdir}/file_list1.txt #! For real
#!images=${myfilesdir}/total_list.txt #! total list
#!images=${myscriptdir}/testfiles.txt #! For debug
images=${myscriptdir}/largefiles.txt #! for the two large files that often get left out

mkdir outputfiles_tum_lymph7
cd outputfiles_tum_lymph7
while read line ; do sbatch ${submit} ${prepare} ${func} "${myfilesdir}/${line}.fits" ; done<${images}
while read line ; do echo ${submit} ${prepare} ${func} "${myfilesdir}/${line}.fits" ; done<${images} #!For debug

cd ..
