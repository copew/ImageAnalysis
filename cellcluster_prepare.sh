#!/bin/bash
#
#PBS -N Matlab
#PBS -m be 
#PBS -k oe
 
date

func=$1
ref=$2

echo "working on file ${ref}"

#INCLUDE MATLAB CALL

#We may have to include -nojvm or there is a memory error
#-nodesktop -nosplash -nodisplay -nojvm together work
#Some Matlab functions like gzip require java so cannot
#use -nojvm option
/usr/local/Cluster-Apps/matlab/R2017b/bin/matlab -nodesktop -nosplash -nodisplay <<EOF

[pa,af,~]=fileparts('${func}');
addpath(pa);
disp(['Path is ' pa])
disp(['Function is ' af])
% disp([af,'(''','${ref}',''')']) %For submitting as string
dofunc=sprintf('%s(''%s'')',af,'${ref}')

% dofunc=sprintf('%s(%s)',af,'${ref}')

eval(dofunc)
;exit
EOF
