%% looking at spatial relationship between tumour clusters and lymphocyte clusters
% all the indexing should be based on the tumour clusters

%image_list=[603288, 593987, 619857, 619872, 619905, 625951];

% this is the full list
%image_list = [603283,603279,603271, 603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966];%,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160];

% this is the running list - to start
image_list = [625946,625951,626018,626047,626102,626103,626160];




%done
% 626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,
% 619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926]; %,
%593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105, 604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922]; %,603262,603269,603263,603257,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966, 602945,602951,602942,602935,602927,602921,602915,593960,593708,593978]; %,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,603283,603279, 603271, 603253, 602994, 604094, 603952, 603952, 594116,594110,595806,594137]; %





src_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph';
output_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap';

X_ind = 1;
Y_ind = 2;
cell_ind = 3;


%% now the loop through the images

for image = 1:size(image_list, 2)
    % loading data
    
    image_filenumber = image_list(image);
    
    
    load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_polygon.mat']);
    load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_buffer_in.mat']);
    load([src_dir '/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_polygon_in.mat']);
    %the data fiels are: tumour_polygon_in, tumour_buffer_in and lymphocyte_polygon
    
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis/mat_file_new/' num2str(image_filenumber) '.mat']);
    
    data_mini = num2cell(data_mini, [1 3]);
    
    % clear variables
    t_l_intersection = [];
    t_l_intersection_area = [];
    t_l_intersection_lymph_count = [];
    t_l_centroid = [];
    t_l_centroid_count=[];
    tbuffer_l_intersection=[];
    tbuffer_l_intersection_area=[];
    tbuffer_l_intersection_lymph_count = [];
    tbuffer_l_centroid = [];
    tbuffer_l_centroid_count=[];
    
    
    %% using tumour polygons as the index/base
    
    % first of all look at overlaps between tumour polygons and lymphocyte polygons
    
    for i=1:size(tumour_polygon_in, 2)
        if isempty(tumour_polygon_in{i})
            t_l_intersection{i}=[];
            t_l_intersection_area{i}=[];
            t_l_intersection_lymph_count{i} = [];
            continue
        end
        for j=1:size(tumour_polygon_in{i}, 2)
            
            if isempty(tumour_polygon_in{i}{j})
                t_l_intersection{i}{j} = [];
                t_l_intersection_area{i}{j} = [];
                t_l_intersection_lymph_count{i}{j} = [];
                continue
            end
            
            %             if isempty(lymphocyte_polygon{i})
            %                 t_l_intersection{i} = [];
            %                 t_l_intersection_area{i} = [];
            %                 t_l_intersection_lymph_count{i} = [0];
            %                 t_l_centroid_count{i} = [0];
            %                 continue
            %             end
            
            t_l_intersection{i}{j}=[];
            t_l_intersection_area{i}{j}=[0];
            t_l_intersection_lymph_count{i}{j} = [0];
            t_l_centroid{i}{j}=[];
            t_l_centroid_count{i}{j}=[0];
            
            
            
            %clear the centroid tmp measurements
            tmp_l_centroid_count = [];
            tmp_l_centroid = [];
            
            
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
            
            
            
            for k = 1:size(lymphocyte_polygon{i}, 2)
                if isempty(lymphocyte_polygon{i}{k})
                    continue
                end
                tmp = intersect(tumour_polygon_in{i}{j}, lymphocyte_polygon{i}{k});
                
                if tmp.NumRegions == 0
                    continue
                end
                tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
                tmp_lymph_count = sum(tmp_lymph);
                t_l_intersection{i}{j} = [t_l_intersection{i}{j} tmp];
                t_l_intersection_area{i}{j} = t_l_intersection_area{i}{j}+area(tmp);
                t_l_intersection_lymph_count{i}{j} = t_l_intersection_lymph_count{i}{j} + tmp_lymph_count;
                
                %whilst still looping through tumour polygons, calculate where the lymphoid
                %centroid are and see whether it fit
                [x,y] = centroid(lymphocyte_polygon{i}{k});
                tmp_l_centroid = [tmp_l_centroid; [x,y]];
                tmp_centroid_in = inpolygon(tmp_l_centroid(:,1), tmp_l_centroid(:,2), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
                tmp_l_centroid_count = sum(tmp_centroid_in);
                t_l_centroid_count{i}{j} = t_l_centroid_count{i}{j} + tmp_l_centroid_count;
            end
            
            
        end %end of lymphocyte polygon{i}
    end %end of tumour polygon{i}
    
    save([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_t_l_intersection2.mat'], 't_l_intersection');
    
    t_l_intersection_lymph_count_cell = extractmycells(t_l_intersection_lymph_count);
    t_l_intersection_area_cell = extractmycells(t_l_intersection_area);
    t_l_centroid_count_cell = extractmycells(t_l_centroid_count);
    
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_intersection_count2.csv'], t_l_intersection_lymph_count_cell);
    csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_lymph_intersection_count2.csv'], t_l_intersection_lymph_count_cell);
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersection_area2.csv'], t_l_intersection_area_cell);
    csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/area/' num2str(image_filenumber) '_intersection_area2.csv'], t_l_intersection_area_cell);
    
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_centroidbuffer_count.csv'], t_l_centroid_count_cell);
    csvwrite([ '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_centroidbuffer_count.csv'], t_l_centroid_count_cell);
    
    
    
    %% this is added to work out the overlap or intersection, between lymphocyte and buffer zones
    
    
    for i=1:size(tumour_buffer_in, 2)
        if isempty(tumour_buffer_in{i})
            tbuffer_l_intersection{i}=[];
            tbuffer_l_intersection_area{i}=[];
            tbuffer_l_intersection_lymph_count{i} = [];
            continue
        end
        
        for j=1:size(tumour_buffer_in{i}, 2)
            if isempty(tumour_buffer_in{i}{j})
                tbuffer_l_intersection{i}{j} = [];
                tbuffer_l_intersection_area{i}{j} = [];
                tbuffer_l_intersection_lymph_count{i}{j} = [];
                continue
            end
            
            tbuffer_l_intersection{i}{j}=[];
            tbuffer_l_intersection_area{i}{j}=[0];
            tbuffer_l_intersection_lymph_count{i}{j} = [0];
            tbuffer_l_centroid{i}{j}=[];
            tbuffer_l_centroid_count{i}{j}=[0];
            
            %clear the centroid tmp measurements
            tmp_l_centroid_count = [];
            tmp_l_centroid = [];
            
            for k = 1:size(lymphocyte_polygon{i}, 2)
                %                 if isempty(lymphocyte_polygon{i})
                %                     tbuffer_l_intersection{i} = [];
                %                     tbuffer_l_intersection_area{i} = [0];
                %                     tbuffer_l_intersection_lymph_count{i} = [0];
                %                     continue
                %                 end
                if isempty(lymphocyte_polygon{i}{k})
                    continue
                end
                
                tmp = intersect(tumour_buffer_in{i}{j}, lymphocyte_polygon{i}{k});
                
                if tmp.NumRegions == 0
                    continue
                end
                tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
                tmp_lymph_count = sum(tmp_lymph);
                tbuffer_l_intersection{i}{j} = [tbuffer_l_intersection{i}{j} tmp];
                tbuffer_l_intersection_area{i}{j} = tbuffer_l_intersection_area{i}{j}+area(tmp);
                tbuffer_l_intersection_lymph_count{i}{j} = tbuffer_l_intersection_lymph_count{i}{j} + tmp_lymph_count;
                
                [x,y] = centroid(lymphocyte_polygon{i}{k});
                tmp_l_centroid = [tmp_l_centroid; [x,y]];
                tmp_centroid_in = inpolygon(tmp_l_centroid(:,1), tmp_l_centroid(:,2), tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
                tmp_l_centroid_count = sum(tmp_centroid_in);
                tbuffer_l_centroid_count{i}{j} = tbuffer_l_centroid_count{i}{j} + tmp_l_centroid_count;
                
                
                
            end
        end
    end
    
    save([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_tbuffer_l_intersection2.mat'], 'tbuffer_l_intersection');
    
    tbuffer_l_intersection_lymph_count_cell = extractmycells(tbuffer_l_intersection_lymph_count);
    tbuffer_l_intersection_area_cell = extractmycells(tbuffer_l_intersection_area);
    tbuffer_l_centroid_count_cell = extractmycells(tbuffer_l_centroid_count);
    
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphbuffer_intersection_count2.csv'], tbuffer_l_intersection_lymph_count_cell);
    csvwrite([ '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_lymphbuffer_intersection_count2.csv'], tbuffer_l_intersection_lymph_count_cell);
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersectionbuffer_area2.csv'], tbuffer_l_intersection_area_cell);
    csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/area/' num2str(image_filenumber) '_intersectionbuffer_area2.csv'], tbuffer_l_intersection_area_cell);
    csvwrite([output_dir '/outputfiles_tum_lymph/' num2str(image_filenumber) '/' num2str(image_filenumber) '_centroidbuffer_count.csv'], tbuffer_l_centroid_count_cell);
    csvwrite([ '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_centroidbuffer_count.csv'], tbuffer_l_centroid_count_cell);
    
    
end %end of tumour polygon
