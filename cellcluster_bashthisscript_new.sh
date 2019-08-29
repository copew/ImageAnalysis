#!/bin/bash

myscriptdir=/rds/user/ww234/hpc-work/ImageAnalysis
myfilesdir=/rds/user/ww234/hpc-work/itpt

submit=${myscriptdir}/cellcluster_submit.sh
prepare=${myscriptdir}/cellcluster_prepare.sh
func=${myscriptdir}/tum_lymph_cluster_hpc.m
#! images=${myfilesdir}/file_list1.txt #! For real
images=${myscriptdir}/testfiles.txt #! For debug


mkdir outputfiles_tum_lymph
cd outputfiles_tum_lymph
while read line ; do sbatch ${submit} ${prepare} ${func} "${myfilesdir}/${line}.fits" ; done<${images}
while read line ; do echo ${submit} ${prepare} ${func} "${myfilesdir}/${line}.fits" ; done<${images} #!For debug

cd ..
