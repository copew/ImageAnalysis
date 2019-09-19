% this is because in the original script i have forgotten to assess the overlap between
% the tumour buffer zone and the lymphocyte clusters

%ok I am stupid - the list needs to be run now i've changed file names.  the ones in the
%done list need to be completely re-run as the other outputs are messed up. the ones in
%the image_list now can be run with the new filename.  14th sep 11:45

%% make a list of files and directories, and set parameters
%image_list = [625951,626018,626047,626102,626103,626160,603288, 603298,603292];

%image_list = [603283,603279,603262,603269,603263,603257,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945]; %this is the list that was caused by bad naming strategy....

image_list = [602951, 619857, 619872, 619905]; %this is a list that is somehow missing from this script but is in the core list.... a bit weird


% done: 603283,603279,603262,603269,603263,603257,603253,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,,,

%done with correct file name: 602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,
% ones with issues: 603271,603253, 602994,604094, 619877, 626172, 625946

src_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph';
output_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap';

X_ind = 1;
Y_ind = 2;
cell_ind = 3;


%% now make the actual calculation

for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    
    load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_polygon.mat']);
    load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_buffer_in.mat']);
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis/mat_file_new/' num2str(image_filenumber) '.mat']);
    data_mini = num2cell(data_mini, 3);
    if size(lymphocyte_polygon, 2) ~= size(tumour_buffer_in, 2)
        warning(['the lymphocyte and tumour polygons have different number of cores ', num2str(image_filenumber)]);
    end
    
    tumour_core_total = size(tumour_buffer_in, 2);
    lymph_core_total = size(lymphocyte_polygon, 2);
    
    % if same sizes then dont' worry about it....
    if lymph_core_total == tumour_core_total
        %do nothing
    else
        diff = lymph_core_total - tumour_core_total;
        if diff > 0
            for d = 1:diff
                tumour_buffer_in{1, tumour_core_total+d} = [];
                
                
            end
        else
            for d = 1:abs(diff)
                lymphocyte_polygon{1, lymph_core_total+d} = [];
            end
        end
    end
    
    tbuffer_l_intersection = [];
    
    for i=1:size(tumour_buffer_in, 2) %looking at the cores with tumours
        tbuffer_l_intersection{i}=[];
        tbuffer_l_intersection_area{i}=[];
        tbuffer_l_intersection_lymph_count{i} = [];
        
        for j=1:size(tumour_buffer_in{i}, 2)
            if isempty(tumour_buffer_in{i}{j})
                continue
            end
            for k = 1:size(lymphocyte_polygon{i}, 2)
                if isempty(lymphocyte_polygon{i}{k})
                    continue
                end
                tmp = intersect(tumour_buffer_in{i}{j}, lymphocyte_polygon{i}{k});
                if tmp.NumRegions == 0
                    continue
                end
                tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
                tmp_lymph_count = sum(tmp_lymph);
                tbuffer_l_intersection{i} = [tbuffer_l_intersection{i} tmp];
                tbuffer_l_intersection_area{i} = [tbuffer_l_intersection_area{i} area(tmp)];
                tbuffer_l_intersection_lymph_count{i} = [tbuffer_l_intersection_lymph_count{i} tmp_lymph_count];
            end
        end
    end
    
    save([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tbuffer_l_intersection.mat'], 'tbuffer_l_intersection');
    csvwrite([output_dir '/count/' num2str(image_filenumber) '_lymphbuffer_intersection_count.csv'], tbuffer_l_intersection_lymph_count);
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphbuffer_intersection_count.csv'], tbuffer_l_intersection_lymph_count);%saving in both locations
    csvwrite([output_dir '/area/' num2str(image_filenumber) '_intersectionbuffer_area.csv'], tbuffer_l_intersection_area);
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersectionbuffer_area.csv'], tbuffer_l_intersection_area);

    
end

% %% This bit is to fix a problem caused by forgotten to rename the file....
% 
% % XXX still need to redo the image_list
% 
% image_list = [603283,603279,603262,603269,603263,603257,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945];
%  
% src_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph';
% output_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap';
% 
% X_ind = 1;
% Y_ind = 2;
% cell_ind = 3;
% 
% 
% for image = 1:size(image_list,2)
%     image_filenumber = image_list(image);
%     
%     load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_polygon.mat']);
%     load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_polygon_in.mat']);
%     load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis/mat_file/' num2str(image_filenumber) '.mat']);
%     data_mini = num2cell(data_mini, 3);
%     if size(lymphocyte_polygon, 2) ~= size(tumour_polygon_in, 2)
%         warning(['the lymphocyte and tumour polygons have different number of cores ', num2str(image_filenumber)]);
%     end
%     
%     tumour_core_total = size(tumour_polygon_in, 2);
%     lymph_core_total = size(lymphocyte_polygon, 2);
%     
%     % if same sizes then dont' worry about it....
%     if lymph_core_total == tumour_core_total
%         %do nothing
%     else
%         diff = lymph_core_total - tumour_core_total;
%         if diff > 0
%             for d = 1:diff
%                 tumour_polygon_in{1, tumour_core_total+d} = [];
%                 
%                 
%             end
%         else
%             for d = 1:abs(diff)
%                 lymphocyte_polygon{1, lymph_core_total+d} = [];
%             end
%         end
%     end
%     
%     t_l_intersection = [];
%     
%     for i=1:size(tumour_polygon_in, 2) %looking at the cores with tumours
%         t_l_intersection{i}=[];
%         t_l_intersection_area{i}=[];
%         t_l_intersection_lymph_count{i} = [];
%         
%         for j=1:size(tumour_polygon_in{i}, 2)
%             if isempty(tumour_polygon_in{i}{j})
%                 continue
%             end
%             for k = 1:size(lymphocyte_polygon{i}, 2)
%                 if isempty(lymphocyte_polygon{i}{k})
%                     continue
%                 end
%                 tmp = intersect(tumour_polygon_in{i}{j}, lymphocyte_polygon{i}{k});
%                 if tmp.NumRegions == 0
%                     continue
%                 end
%                 tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
%                 tmp_lymph_count = sum(tmp_lymph);
%                 t_l_intersection{i} = [t_l_intersection{i} tmp];
%                 t_l_intersection_area{i} = [t_l_intersection_area{i} area(tmp)];
%                 t_l_intersection_lymph_count{i} = [t_l_intersection_lymph_count{i} tmp_lymph_count];
%             end
%         end
%     end
%     
%     save([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_t_l_intersection.mat'], 't_l_intersection');
%     csvwrite([output_dir '/count/' num2str(image_filenumber) '_lymph_intersection_count.csv'], t_l_intersection_lymph_count);
%     csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_intersection_count.csv'], t_l_intersection_lymph_count);%saving in both locations
%     csvwrite([output_dir '/area/' num2str(image_filenumber) '_intersection_area.csv'], t_l_intersection_area);
%     csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersection_area.csv'], t_l_intersection_area);
% 
%     
% end