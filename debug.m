save(['./' num2str(image_filenumber) '/tumour_polygon_in.mat'], 'tumour_polygon_in');
save(['./' num2str(image_filenumber) '/tumour_buffer_in.mat'], 'tumour_buffer_in');
save(['./' num2str(image_filenumber) '/lymphocyte_polygon.mat'], 'lymphocyte_polygon');

%% looking at spatial relationship between tumour clusters and lymphocyte clusters


% first of all look at overlaps
% need to select the cores with both tumour polygons and lymphocyte polygons
overlap_list = intersect(tumour_core_list, lymphocyte_core_list);
t_l_intersection = [];

%%
for i=1:size(overlap_list, 2) %looking at the cores with tumours
    t_l_intersection{overlap_list(i)}=[];
    t_l_intersection_area{overlap_list(i)}=[];
    t_l_intersection_lymph_count{overlap_list(i)} = [];
    
    for j=1:size(tumour_polygon_in{overlap_list(i)}, 2)
        if isempty(tumour_polygon_in{overlap_list(i)}{j})
            continue
        end
        for k = 1:size(lymphocyte_polygon{overlap_list(i)}, 2)
            if isempty(lymphocyte_polygon{overlap_list(i)}{k})
                continue
            end
            tmp = intersect(tumour_polygon_in{overlap_list(i)}{j}, lymphocyte_polygon{overlap_list(i)}{k});
            if tmp.NumRegions == 0
                continue
            end
            tmp_lymph_count = sum(tmp_lymph);
            tmp_lymph = inpolygon(in_core{overlap_list(i)}{X_ind}(in_core{overlap_list(i)}{cell_ind}==2), in_core{overlap_list(i)}{Y_ind}(in_core{overlap_list(i)}{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
            t_l_intersection{overlap_list(i)} = [t_l_intersection{overlap_list(i)} tmp];
            t_l_intersection_area{overlap_list(i)} = [t_l_intersection_area{overlap_list(i)} area(tmp)];
            t_l_intersection_lymph_count{overlap_list(i)} = [t_l_intersection_lymph_count{overlap_list(i)} tmp_lymph_count];
        end
    end
end

save(['./' num2str(image_filenumber) '/t_l_intersection.mat'], 't_l_intersection');
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_intersection_count.csv'], t_l_intersection_lymph_count);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersection_area.csv'], t_l_intersection_area);


%%
%calculate area
for i = 1:size(core_list, 2)
    %core area
    core_area{core_list(i)} = area(core_polygon{i});
    %tumour polygon area
end

for t=1:size(tumour_polygon_in, 2)
    for j = 1:size(tumour_polygon_in{t}, 2)
        tumour_polygon_area{t}{j} = area(tumour_polygon_in{t}{j});
    end
    %tumour buffer zone area
    for k = 1:size(tumour_buffer_in{t},2)
        tumour_buffer_area{t}{k} = area(tumour_buffer_in{t}{k});
    end
end

    %lymphocyte cluster area
    
for x = 1:size(lymphocyte_polygon, 2)   
    for l = 1:size(lymphocyte_polygon{x}, 2)
        lymphocyte_polygon_area{x}{l} = area(lymphocyte_polygon{x}{l});
    end 
end

%create an output file
core_area_combined = horzcat(core_area{:});
%core_area_combined(cellfun(@isempty, core_area_combined))=[];
tumour_in_area_combined=horzcat(tumour_polygon_area{:});
tumour_in_area_combined(cellfun(@isempty, tumour_in_area_combined))=[];
tumour_buffer_area_combined=horzcat(tumour_buffer_area{:});
tumour_buffer_area_combined(cellfun(@isempty, tumour_buffer_area_combined))=[];
lymphocyte_polygon_area_combined=horzcat(lymphocyte_polygon_area{:});
lymphocyte_polygon_area_combined(cellfun(@isempty, lymphocyte_polygon_area_combined))=[];

%turn into csv files
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_area_in.csv'], tumour_in_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_area_buffer.csv'], tumour_buffer_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_core_area.csv'], core_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_area.csv'], lymphocyte_polygon_area_combined);

%% 
%get the lymphocyte counts as well as tumour cell count
for i = 1:size(core_list, 2)
    for j = 1:size(tumour_polygon_in{i}, 2)
        if isempty(tumour_polygon_in{i}{j})
            continue
        end
        lymph_in{i}{j}=inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==2), data_trimmed{Y_ind}(data_trimmed{cell_ind}==2), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        lymph_in_count{i}{j} = sum(lymph_in{i}{j});
        tumourcell_in{i}{j}=inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==1), data_trimmed{Y_ind}(data_trimmed{cell_ind}==1), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        tumourcell_in_count{i}{j} = sum(tumourcell_in{i}{j});
    end
    
    for k = 1:size(tumour_buffer_in{i},2)
        if isempty(tumour_buffer_in{i}{k})
            continue
        end
        lymph_buffer{i}{k} = inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==2), data_trimmed{Y_ind}(data_trimmed{cell_ind}==2), tumour_buffer_in{i}{k}.Vertices(:,1), tumour_buffer_in{i}{k}.Vertices(:,2));
        lymph_buffer_count{i}{k} = sum(lymph_buffer{i}{k});
        tumourcell_buffer{i}{j}=inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==1), data_trimmed{Y_ind}(data_trimmed{cell_ind}==1), tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
        tumourcell_buffer_count{i}{j} = sum(tumourcell_buffer{i}{j});
    end
    
    XXXXXXXXX
        if exist(lymphocyte_polygon{i})==1
            for l = 1:size(lymphocyte_polygon{i},2)
                if isempty(lymphocyte_polygon{i}{l})
                    continue
                end
                lymph_lymphcluster{i}{l} = inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==2), data_trimmed{Y_ind}(data_trimmed{cell_ind}==2), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
                lymph_lymphcluster_count{i}{l} = sum(lymph_lymphcluster{i}{l});
                
            end
     
    
    end
    
    
    
    
end

%create an output file
lymph_in_count_combined=horzcat(lymph_in_count{:});
lymph_buffer_count_combined=horzcat(lymph_buffer_count{:});
tumourcell_in_count_combined=horzcat(tumourcell_in_count{:});
tumourcell_buffer_count_combined=horzcat(tumourcell_buffer_count{:});
lymph_lymphcluster_count_combined = horzcat(lymph_lymphcluster_count{:});

%turn into csv files
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_in_count.csv'], lymph_in_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_buffer_count.csv'], lymph_buffer_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumourcell_in_count.csv'], tumourcell_in_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumourcell_buffer_count.csv'], tumourcell_buffer_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_lymphcluster_count.csv'], lymph_lymphcluster_count_combined);
