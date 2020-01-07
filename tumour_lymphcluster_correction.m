%this is for fixing a stupid mistake in the tum_lymph_cluster_hpc that caused the tumour lymph cluster count to be exactly the same as lymph_lymph cluster count

image_list = [619872, 619857, 619905];
cluster_size = [5];
lymph_cluster_size = [5];
buffer_size = 100;


for image = 1:size(image_list,2)
    image_filenumber = image_list(image);

% loading file

load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph5/' num2str(image_filenumber) '/'  num2str(image_filenumber) '_workspace_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat']);
image_list = [619872, 619857, 619905];

% calculating
for i = 1:size(lymphocyte_polygon, 2)
    if isempty(lymphocyte_polygon{i})
        tumour_lymphcluster_count{i} = [];
        continue
    end
    
    for l = 1:size(lymphocyte_polygon{i},2)
        if isempty(lymphocyte_polygon{i}{l})
            tumour_lymphcluster_count{i}{j} = [];
            continue
        end
        
        tumour_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
        tumour_lymphcluster_count{i}{l} = sum(tumour_lymphcluster{i}{l});
    end
    
end

%create an output file

tumour_lymphcluster_count_combined = extractmycells(tumour_lymphcluster_count);

%turn into csv files

csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/correct_tumour_lymphcluster/' num2str(image_filenumber) '_tumour_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumour_lymphcluster_count_combined);



end