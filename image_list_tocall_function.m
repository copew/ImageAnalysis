% this file is for test running a function in order to run the tumour lymphocyte cluster
% on hpc.  this is to test if the function is working

%this is an abbreviated short list of images
image_list=[619857, 619872, 619905, 625951]; %603288, 593987, 
for this_image = 1:size(image_list, 2)
    image_filenumber = image_list(this_image);
    tum_lymph_cluster_hpc(image_filenumber)
end
