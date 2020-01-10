%this is for fixing a stupid mistake in the tum_lymph_cluster_hpc that caused the tumour lymph cluster count to be exactly the same as lymph_lymph cluster count

image_list = [603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160,603288, 603298,603292, 619872, 619857, 619905];
cluster_size = [5];
lymph_cluster_size = [5];
buffer_size = 100;


for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    
    % loading file
    
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph5/' num2str(image_filenumber) '/'  num2str(image_filenumber) '_workspace_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat']);
    image_list = [603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160,603288, 603298,603292, 619872, 619857, 619905];
    
    % calculating
    for i = 1:size(lymphocyte_polygon, 2)
        if isempty(lymphocyte_polygon{i})
            lymph_lymphcluster_count{i} = [];
            tumour_lymphcluster_count{i} = [];
            continue
        end
        
        for l = 1:size(lymphocyte_polygon{i},2)
            if isempty(lymphocyte_polygon{i}{l})
                lymph_lymphcluster_count{i}{j} = [];
                tumour_lymphcluster_count{i}{j} = [];
                continue
            end
            lymph_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
            lymph_lymphcluster_count{i}{l} = sum(lymph_lymphcluster{i}{l});
            tumour_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
            tumour_lymphcluster_count{i}{l} = sum(tumour_lymphcluster{i}{l});
        end
        
    end
    
    %create an output file
    
    lymph_lymphcluster_count_combined = extractmycells(lymph_lymphcluster_count);
    tumour_lymphcluster_count_combined = extractmycells(tumour_lymphcluster_count);
    
    %turn into csv files
    csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/correct_tumour_lymphcluster/' num2str(image_filenumber) '_lymph_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], lymph_lymphcluster_count_combined);
    csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/correct_tumour_lymphcluster/' num2str(image_filenumber) '_tumour_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumour_lymphcluster_count_combined);
    
    
end



% % % % % % % for i = 1:size(lymphocyte_polygon, 2)
% % % % % % %     if isempty(lymphocyte_polygon{i})
% % % % % % %         tumour_lymphcluster_count{i} = [];
% % % % % % %         continue
% % % % % % %     end
% % % % % % %
% % % % % % %     for l = 1:size(lymphocyte_polygon{i},2)
% % % % % % %         if isempty(lymphocyte_polygon{i}{l})
% % % % % % %             tumour_lymphcluster_count{i}{l} = [];
% % % % % % %             continue
% % % % % % %         end
% % % % % % %
% % % % % % %         tumour_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
% % % % % % %         tumour_lymphcluster_count{i}{l} = sum(tumour_lymphcluster{i}{l});
% % % % % % %     end
% % % % % % %
% % % % % % % end
% % % % % % %
% % % % % % % % create an output file
% % % % % % %
% % % % % % % tumour_lymphcluster_count_combined = extractmycells(tumour_lymphcluster_count);
% % % % % % %
% % % % % % % % turn into csv files
% % % % % % %
% % % % % % % csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/correct_tumour_lymphcluster/' num2str(image_filenumber) '_tumour_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '2.csv'], tumour_lymphcluster_count_combined);
% % % % % % %
% % % % % % %
% % % % % % %
% % % % % % % end