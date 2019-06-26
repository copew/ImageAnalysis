%this is trying to combine all the cell arrays from the knn measurement

image_list = [603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160,603288, 603298,603292];

% %%
% %loading the data
% for i = 1:length(image_list)
% load_buffer=['buffer_' num2str(image_list(i)) '=load(''' num2str(image_list(i)) '_lymph_buffer_distances.mat'',''lymph_buffer_distances'');']
% eval(load_buffer)
% end
% 
% for i = 1:length(image_list)
% load_in=['in_' num2str(image_list(i)) '=load(''' num2str(image_list(i)) '_lymph_in_distances.mat'',''lymph_in_distances'');']
% eval(load_in)
% end
% 
%make a list of filenames
buffer_list=[];
in_list=[];
for i = 1:length(image_list)
    buffer_list{i} = ['buffer_' num2str(image_list(i))];
    in_list{i} = ['in_' num2str(image_list(i))];
end

%%

%for the buffer zones
for i = 1:length(image_list)
    load([num2str(image_list(i)) '_lymph_buffer_distances.mat'],'lymph_buffer_distances');
    lymph_buffer_dist_notempty = lymph_buffer_distances(~cellfun('isempty',lymph_buffer_distances));
    
    for j = 1:size(lymph_buffer_dist_notempty, 2)
        lymph_buffer_dist_notempty{j} = lymph_buffer_dist_notempty{j}(~cellfun('isempty',lymph_buffer_dist_notempty{j}));
    end
    
    final_array = [];
    for k = 1:size(lymph_buffer_dist_notempty, 2)
        for l = 1:size(lymph_buffer_dist_notempty{k}, 2)
            temp_array = nan(51,size(lymph_buffer_dist_notempty{k}{l},2));
            temp_array(1:size(lymph_buffer_dist_notempty{k}{l},1),1:size(lymph_buffer_dist_notempty{k}{l},2))=lymph_buffer_dist_notempty{k}{l};
            final_array = [final_array,temp_array];
        end
    end
    save_statement = [buffer_list{i} '= final_array;'];
    eval(save_statement);
end

%for the in zones
for i = 1:length(image_list)
    load([num2str(image_list(i)) '_lymph_in_distances.mat'],'lymph_in_distances');
    lymph_in_dist_notempty = lymph_in_distances(~cellfun('isempty',lymph_in_distances));
    
    for j = 1:size(lymph_in_dist_notempty, 2)
        lymph_in_dist_notempty{j} = lymph_in_dist_notempty{j}(~cellfun('isempty',lymph_in_dist_notempty{j}));
    end
    
    final_array = [];
    for k = 1:size(lymph_in_dist_notempty, 2)
        for l = 1:size(lymph_in_dist_notempty{k}, 2)
            temp_array = nan(51,size(lymph_in_dist_notempty{k}{l},2));
            temp_array(1:size(lymph_in_dist_notempty{k}{l},1),1:size(lymph_in_dist_notempty{k}{l},2))=lymph_in_dist_notempty{k}{l};
            final_array = [final_array,temp_array];
        end
    end
    save_statement = [in_list{i} '= final_array;'];
    eval(save_statement);

end
    

for i = 1:length(image_list)
    csvwrite([num2str(image_list(i)) '_knn_in_distance.csv'], eval(in_list{i}));
    csvwrite([num2str(image_list(i)) '_knn_buffer_distance.csv'], eval(buffer_list{i}));
end
  
    
    
    