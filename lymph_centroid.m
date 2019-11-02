% this is to calculate centroid value for each of the lymphocyte clusters for subsequent
% analysis

%original image list
% image_list = [593708,593960,593971,593978,593987,594000,594006,594017,594035,594044,594051,594060,594080,594088,594095,594105,594110,594116,594137,595806,597768,597777,597786,597794,597800,597812,597820,597830,597838,599797,599803,599813,599831,599850,602915,602921,602927,602935,602942,602945,602951,602952,602958,602962,602966,602971,602976,602983,602994,603006,603007,603055,603245,603253,603257,603262,603263,603269,603271,603279,603283,603288,603292,603298,603904,603911,603919,603922,603925,603928,603941,603944,603952,603956,603963,603966,603970,603976,603982,603987,604045,604047,604055,604057,604061,604066,604071,604077,604078,604083,604089,604094,604098,604105,607165,607166,619460,619465,619470,619476,619483,619489,619496,619503,619508,619515,619527,619533,619538,619539,619550,619556,619571,619579,619583,619678,619793,619799,619803,619809,619815,619832,619837,619842,619847,619852,619857,619859,619863,619868,619872,619877,619884,619896,619905,619911,619916,619922,619926,619932,619937,619942,619953,619954,625322,625333,625338,625865,625876,625887,625891,625895,625908,625911,625916,625923,625930,625936,625946,625951,625958,626018,626047,626102,626103,626160,626166,626171,626172];

image_list = [593708]; %,593960,593971,593978,593987,594000];

for image = 1:size(image_list,2)
    
    image_filenumber = image_list(image);
    
    %% loading files
    
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_polygon.mat'])

    %% working out centroid
    
    lymphcluster_centroid = [];
    
    for this_core = 1:size(lymphocyte_polygon, 2)
        if isempty(lymphocyte_polygon{this_core})
            continue
        end
        lymphcluster_centroid{this_core} = [];
        for this_polygon = 1:size(lymphocyte_polygon{this_core}, 2)
            [x,y] = centroid(lymphocyte_polygon{this_core}{this_polygon});
            lymphcluster_centroid{this_core} = [lymphcluster_centroid{this_core}; [x,y]];
            
        end
       % lymphcluster_centroid{this_core} = num2cell(lymphcluster_centroid{this_core});
    end
    
    
    %% saving files
    
    
end
