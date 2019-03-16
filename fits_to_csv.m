% This converts .fits files into .csv, and keeping only the X-coordinates, Y-coordinates
% and cell type info for processing in cluster

image_list = [603298,603292,603288,603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160];

for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
%     info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
%     image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];
    
    %create indexing
    X_ind = 3;
    Y_ind = 4;
    cell_ind = size(data,2);
    
    %trim data
    data_trimmed = data;
    for i = 1:size(data, 2)
        data_trimmed{i} = data_trimmed{i}(data{cell_ind}~=0);
    end

    % Now remove duplicate points from cell identification tiling
    full_coords = zeros(1,length(data_trimmed{Y_ind}));
    for i = 1:length(data_trimmed{Y_ind}) 
        full_coords(i) = [data_trimmed{X_ind}(i)*1000000000 + data_trimmed{Y_ind}(i)]; 
    end
    [~, unique_coords, ~] = unique(full_coords,'stable');
    for i = 1:length(data_trimmed)
        data_trimmed{i} = data_trimmed{i}(unique_coords);
    end
    
    % and combine tumour and normal
    data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;
    
    data_mini = [data_trimmed{X_ind} data_trimmed{Y_ind} data_trimmed{cell_ind}];
    
    %convert to csv file
    %csvwrite([num2str(image_filenumber) '_mini.csv'], data_mini)
    
    %save as .mat file
    save(['/Users/cope01/Documents/MATLAB/mat_file/' num2str(image_filenumber) '.mat'], 'data_mini');
end
    