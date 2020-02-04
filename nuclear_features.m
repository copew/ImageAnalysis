%this is a total list
image_list =[603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160, 603288, 603298,603292, 619872, 619857, 619905];

% this is for testing
%image_list = [603283];

for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    image_path_stem = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis/IT_PT_zone';
    
    % loading files
    data = fitsread([image_path_stem '/' num2str(image_filenumber) '.fits'],'binarytable');
    info = fitsinfo([image_path_stem '/' num2str(image_filenumber) '.fits']);
    
    % create indexing
    X_ind = 3;
    Y_ind = 4;
    overlap = 22;
    s2n=61;
    cell_ind = 62;
    
    
    %trim data
    data_trimmed = data;
    
    
    %remove non cells
    for i = 1:size(data, 2)
        data_trimmed{i} = data_trimmed{i}(data{cell_ind}~=0);
    end
    %remove overlaps
    data_tmp = data_trimmed;
    for i = 1:size(data_trimmed, 2)
        data_trimmed{i} = data_trimmed{i}(data_tmp{overlap} =='F');
    end
    %remove low signal to noise ratio (set threshold at 1.3 - discussed with A Dariush)
    data_tmp = data_trimmed;
    for i = 1:size(data_trimmed, 2)
        data_trimmed{i} = data_trimmed{i}(data_tmp{s2n} >= 1.3);
    end
    
    % and combine tumour and normal
    data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;
    
    
    
    % so now data_trimmed is the array that we need
    
    %now removing column 58 (slide name), 53 & 52 (pixel:micron ratio), 22 (overlap), 21
    %(convexity)
    for j = 1:size(data_trimmed, 2)
        if isa(data_trimmed{j}, 'double')
            continue
        elseif isa(data_trimmed{j}, 'char')
            data_trimmed{j} = double(data_trimmed{j});
        else
             data_trimmed{j} = cell2mat(data_trimmed{j});
             data_trimmed{j} = double(data_trimmed{j});
        end
        
    end
    
    %create subset
    data_trimmed_t = [];
    data_trimmed_l = [];
    data_trimmed_s = [];
    
    for i = 1:size(data_trimmed, 2)
        data_trimmed_t{i} = data_trimmed{i}(data_trimmed{cell_ind} == 1);
        data_trimmed_l{i} = data_trimmed{i}(data_trimmed{cell_ind} == 2);
        data_trimmed_s{i} = data_trimmed{i}(data_trimmed{cell_ind} == 3);
        
    end
    
    
    %rows will be as such:
    % %     mean = 1;
    % %     median=2;
    % %     min=3;
    % %     max=4;
    % %     stdev=5;
    % %     firstq=6;
    % %     thirdq=7;
    % %     iqr=8;
    
    
    summary_statistics_t = zeros(8, size(data_trimmed_t, 2));
    summary_statistics_l = zeros(8, size(data_trimmed_l, 2));
    summary_statistics_s = zeros(8, size(data_trimmed_s, 2));
    
    for k = 1:size(data_trimmed_t, 2)
        summary_statistics_t(1,k) = mean(data_trimmed_t{k});
        summary_statistics_t(2,k) = median(data_trimmed_t{k});
        summary_statistics_t(3,k) = min(data_trimmed_t{k});
        summary_statistics_t(4,k) = max(data_trimmed_t{k});
        summary_statistics_t(5,k) = std(data_trimmed_t{k});
        summary_statistics_t(6,k) = prctile(data_trimmed_t{k}, 25);
        summary_statistics_t(7,k) = prctile(data_trimmed_t{k}, 75);
        summary_statistics_t(8,k) = iqr(data_trimmed_t{k});
    end
    
    for k = 1:size(data_trimmed_l, 2)
        summary_statistics_l(1,k) = mean(data_trimmed_l{k});
        summary_statistics_l(2,k) = median(data_trimmed_l{k});
        summary_statistics_l(3,k) = min(data_trimmed_l{k});
        summary_statistics_l(4,k) = max(data_trimmed_l{k});
        summary_statistics_l(5,k) = std(data_trimmed_l{k});
        summary_statistics_l(6,k) = prctile(data_trimmed_l{k}, 25);
        summary_statistics_l(7,k) = prctile(data_trimmed_l{k}, 75);
        summary_statistics_l(8,k) = iqr(data_trimmed_l{k});
    end
    
    for k = 1:size(data_trimmed_s, 2)
        summary_statistics_s(1,k) = mean(data_trimmed_s{k});
        summary_statistics_s(2,k) = median(data_trimmed_s{k});
        summary_statistics_s(3,k) = min(data_trimmed_s{k});
        summary_statistics_s(4,k) = max(data_trimmed_s{k});
        summary_statistics_s(5,k) = std(data_trimmed_s{k});
        summary_statistics_s(6,k) = prctile(data_trimmed_s{k}, 25);
        summary_statistics_s(7,k) = prctile(data_trimmed_s{k}, 75);
        summary_statistics_s(8,k) = iqr(data_trimmed_s{k});
    end
 
  csvwrite(['./nuclear_features/' num2str(image_filenumber) '_ss_t.csv'], summary_statistics_t);
  csvwrite(['./nuclear_features/' num2str(image_filenumber) '_ss_l.csv'], summary_statistics_l);
  csvwrite(['./nuclear_features/' num2str(image_filenumber) '_ss_s.csv'], summary_statistics_s);
    
end
