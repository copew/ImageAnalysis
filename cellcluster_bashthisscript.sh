myscriptdir=/rds-d4/user/ww234/hpc-work/ImageAnalysis
myfilesdir=/rds-d4/user/ww234/hpc-work/itpt

submit=${myscriptdir}/cellcluster_submit_himemlongtime.sh
prepare=${myscriptdir}/cellcluster_prepare.sh
func=${myscriptdir}/cluster_buffer_knn_hpc.m
images=${myfilesdir}/file_list2.txt

mkdir outputfiles
cd outputfiles
while read line ; do sbatch ${submit} ${prepare} ${func} ${line} ; done<${images}
#!this line submit all the jobs in the file list

#!while read line ; do echo ${submit} ${prepare} ${func} ${line} ; done<${images} #! for debug this will list all the sbatch submissions, to test, copy one line and paste after "sbatch "
cd ..

