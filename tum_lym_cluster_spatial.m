%% looking at spatial relationship between tumour clusters and lymphocyte clusters
% all the indexing should be based on the tumour clusters

image_list=[603288, 593987, 619857, 619872, 619905, 625951];


src_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph';
output_dir = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap';

X_ind = 1;
Y_ind = 2;
cell_ind = 3;


%% now the loop through the images

for i = 1:size(image_list, 2)



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
tbuffer_l_intersection=[];
tbuffer_l_intersection_area=[];
tbuffer_l_intersection_lymph_count = [];


%% using tumour polygons as the index/base

% first of all look at overlaps between tumour polygons and lymphocyte polygons



for i=1:size(tumour_polygon_in, 2)
    
    for j=1:size(tumour_polygon_in{i}, 2)
        t_l_intersection{i}{j} = [];
        t_l_intersection_area{i}{j} = [0];
        t_l_intersection_lymph_count{i}{j}= [0];
%         tbuffer_l_intersection{i}{j}=[];
%         tbuffer_l_intersection_area{i}{j}=[0];
%         tbuffer_l_intersection_lymph_count{i}{j} = [0]
        
        
        if isempty(tumour_polygon_in{i})
            t_l_intersection{i} = [];
            t_l_intersection_area{i} = [];
            t_l_intersection_lymph_count{i} = [];
%             tbuffer_l_intersection{i}=[];
%             tbuffer_l_intersection_area{i}=[];
%             tbuffer_l_intersection_lymph_count{i} = [];
            continue
        end
        if isempty(tumour_polygon_in{i}{j})
            t_l_intersection{i}{j} = [];
            t_l_intersection_area{i}{j} = [];
            t_l_intersection_lymph_count{i}{j} = [];

            continue
        end
      
        for k = 1:size(lymphocyte_polygon{i}, 2)
            if isempty(lymphocyte_polygon{i})
                t_l_intersection{i} = [];
                t_l_intersection_area{i} = [];
                t_l_intersection_lymph_count{i} = [];

                continue
            end
            if isempty(lymphocyte_polygon{i}{k})
                continue
            end
            
            tmp = intersect(tumour_polygon_in{i}{j}, lymphocyte_polygon{i}{k});
            
            if tmp.NumRegions == 0
                continue
            end
            tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
            tmp_lymph_count = sum(tmp_lymph);
            t_l_intersection{i}{j} = [t_l_intersection{i}{j}, tmp];
            t_l_intersection_area{i}{j} = [t_l_intersection_area{i}{j}+area(tmp)];
            t_l_intersection_lymph_count{i}{j} = [t_l_intersection_lymph_count{i}{j} + tmp_lymph_count];
        end
    end
end

% save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) 'test_t_l_intersection.mat'], 't_l_intersection');
% csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) 'test_lymph_intersection_count.csv'], t_l_intersection_lymph_count);
% csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) 'test_intersection_area.csv'], t_l_intersection_area);

end
% 
% %% this is added to work out the overlap or intersection, between lymphocyte and buffer zones
%                 tbuffer_l_intersection{i}=[];
%                 tbuffer_l_intersection_area{i}=[];
%                 tbuffer_l_intersection_lymph_count{i} = [];
%             tbuffer_l_intersection{i}{j}=[];
%             tbuffer_l_intersection_area{i}{j}=[];
%             tbuffer_l_intersection_lymph_count{i}{j} = [];
% 
% 
% 
% tbuffer_l_intersection = [];
% 
% for i=1:size(overlap_list, 2) %looking at the cores with tumours
%     tbuffer_l_intersection{overlap_list(i)}=[];
%     tbuffer_l_intersection_area{overlap_list(i)}=[];
%     tbuffer_l_intersection_lymph_count{overlap_list(i)} = [];
%     
%     for j=1:size(tumour_buffer_in{overlap_list(i)}, 2)
%         if isempty(tumour_buffer_in{overlap_list(i)}{j})
%             continue
%         end
%         for k = 1:size(lymphocyte_polygon{overlap_list(i)}, 2)
%             if isempty(lymphocyte_polygon{overlap_list(i)}{k})
%                 continue
%             end
%             tmp = intersect(tumour_buffer_in{overlap_list(i)}{j}, lymphocyte_polygon{overlap_list(i)}{k});
%             if tmp.NumRegions == 0
%                 continue
%             end
%             tmp_lymph = inpolygon(in_core{overlap_list(i)}{X_ind}(in_core{overlap_list(i)}{cell_ind}==2), in_core{overlap_list(i)}{Y_ind}(in_core{overlap_list(i)}{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
%             tmp_lymph_count = sum(tmp_lymph);
%             tbuffer_l_intersection{overlap_list(i)} = [tbuffer_l_intersection{overlap_list(i)} tmp];
%             tbuffer_l_intersection_area{overlap_list(i)} = [tbuffer_l_intersection_area{overlap_list(i)} area(tmp)];
%             tbuffer_l_intersection_lymph_count{overlap_list(i)} = [tbuffer_l_intersection_lymph_count{overlap_list(i)} tmp_lymph_count];
%         end
%     end
% end
% 
% save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tbuffer_l_intersection.mat'], 'tbuffer_l_intersection');
% csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphbuffer_intersection_count.csv'], tbuffer_l_intersection_lymph_count);
% %csvwrite([ '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_lymphbuffer_intersection_count.csv'], tbuffer_l_intersection_lymph_count);
% csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersectionbuffer_area.csv'], tbuffer_l_intersection_area);
% %csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/area/' num2str(image_filenumber) '_intersectionbuffer_area.csv'], tbuffer_l_intersection_area);
% 
% 
% 
% %% this is to work out the relationship between the centroid of the lymphoid clusters and the tumour clusters
