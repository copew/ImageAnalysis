function scp_get_files(thisnum)
eval(['!scp ww234@login-cpu.hpc.cam.ac.uk:/rds-d4/user/ww234/hpc-work/ImageAnalysis/outputfiles/' num2str(thisnum) '* .'])
