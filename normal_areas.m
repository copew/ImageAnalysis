% load the matlab workspace


image_list= [603288];
for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/outputfiles_tum_lymph7/' num2str(image_filenumber) '/' num2str(image_filenumber) '_workspace_l20_b100.mat']);

    for i = 1:size(core_list, 2)
        this_core = core_list(i); 
        normal_polygon{this_core} = substract(core_polygon{this_core}, tumour_polygon_in{this_core});
        normal_polygon{this_core} = substract(normal_polygon{this_core}, tumour_buffer_in{this_core});
        lymph_normal_logi{this_core} = inpolygon(in_core{this_core}(:,1), in_core{this_core}(:,2), normal_polygon{this_core}.Vertices(:,1), normal_polygon{this_core}.Vertices(:,2));
        lymph_normal_count{this_core} = sum(lymph_normal_logi{this_core});
    
    
    
    
end













%% looking at spatial relationship between tumour clusters, buffer and lymphocyte clusters

% indexing using tumour polygons

t_l_intersection = {};
t_l_intersection_area={};
t_l_intersection_lymph_count_logi={};
t_l_intersection_lymph_count= {};
t_l_centroid_count = {};

tbuffer_l_intersection = {};
tbuffer_l_intersection_area={};
tbuffer_l_intersection_lymph_count_logi={};
tbuffer_l_intersection_lymph_count= {};
tbuffer_l_centroid_count = {};

lymph_centroid_in = {};
lymph_centroid_buffer = {};
lymph_centroid_in_logi = {};
lymph_centroid_buffer_logi = {};

for i=1:size(tumour_polygon_in2, 2) %selecting cores with tumour polygons
    if isempty(tumour_polygon_in2{i}) %in case some cores dont' have tumour
        continue
    end
    %         t_l_intersection{i} = [];
    %         t_l_intersection_area{i}={};
    %         t_l_intersection_lymph_count{i}= {};
    %         t_l_centroid_count{i} = {};
    
    if isempty(lymphocyte_polygon{i}) %if there is no lymphocyte in this core then all value for this becomes empty
        t_l_intersection{i} = cell(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        zeromat = zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        zerocell = mat2cell(zeromat, 1, ones(1, size(zeromat, 2)));
        t_l_intersection_area{i} = zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        t_l_intersection_lymph_count{i}= zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        t_l_centroid_count{i} = zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        tbuffer_l_intersection{i} = cell(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        tbuffer_l_intersection_area{i} = zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        tbuffer_l_intersection_lymph_count{i}= zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        tbuffer_l_centroid_count{i} = zerocell; %zeros(size(tumour_polygon_in2{i}, 1), size(tumour_polygon_in2{i}, 2));
        %this has made the corresponding cell for polygons, area and counts as
        %empty or zero respectively
        continue
    end
    
    for j = 1:size(tumour_polygon_in2{i}, 2) %indexing with tumour polygons
        t_l_intersection{i}{j} = polyshape(); %make sure it's empty
        t_l_intersection_area{i}{j} = [];
        t_l_intersection_lymph_count_logi{i}{j}= [];
        t_l_intersection_lymph_count{i}{j}  = [];
        t_l_centroid_count{i}{j} = [];
        tbuffer_l_intersection{i}{j} = polyshape(); %make sure it's empty
        tbuffer_l_intersection_area{i}{j} = [];
        tbuffer_l_intersection_lymph_count_logi{i}{j}=[];
        tbuffer_l_intersection_lymph_count{i}{j}  = [];
        tbuffer_l_centroid_count{i}{j} = [];
        lymph_centroid_in{i}{j} = [];
        lymph_centroid_buffer{i}{j} = [];
        lymph_centroid_in_logi{i}{j} = [];
        lymph_centroid_buffer_logi{i}{j} = [];
        
        if isempty(tumour_polygon_in2{i}{j})
            continue
        end
        
        
        for k = 1:size(lymphocyte_polygon{i}, 2)
            %now looking at the overlap between the particular lymphocyte polygon and
            %tumour polygon
            tmp = intersect(tumour_polygon_in2{i}{j}, lymphocyte_polygon{i}{k});
            %                 if tmp.NumRegions == 0
            %                     %if no actual overlap at all then move on
            %                     continue
            %                 end
            t_l_intersection{i}{j} = union(t_l_intersection{i}{j}, tmp); %a collection of polygons of t_l overlap
            [x,y] = centroid(lymphocyte_polygon{i}{k}); %a collection of centroid of the lymphoid polygons that have overlaps with this tumour cluster
            lymph_centroid_in{i}{j} = [lymph_centroid_in{i}{j}; [x,y]];
        end
        
        if isempty(t_l_intersection{i}{j})
            t_l_intersection_lymph_count{i}{j} = 0;
            t_l_intersection_area{i}{j} = 0;
            continue
        end
        
        t_l_intersection_lymph_count_logi{i}{j} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), t_l_intersection{i}{j}.Vertices(:,1), t_l_intersection{i}{j}.Vertices(:,2));
        t_l_intersection_lymph_count{i}{j} = sum(t_l_intersection_lymph_count_logi{i}{j});
        t_l_intersection_area{i}{j} = area(t_l_intersection{i}{j});
        
        if isempty(lymph_centroid_in{i}{j})
            lymph_centroid_in_logi{i}{j} = 0;
            t_l_centroid_count{i}{j} = 0;
            continue
        end
        
        lymph_centroid_in_logi{i}{j} = inpolygon(lymph_centroid_in{i}{j}(:,1), lymph_centroid_in{i}{j}(:,2), tumour_polygon_in2{i}{j}.Vertices(:,1), tumour_polygon_in2{i}{j}.Vertices(:,2));
        t_l_centroid_count{i}{j} = sum(lymph_centroid_in_logi{i}{j});
        
        
        for l = 1:size(lymphocyte_polygon{i}, 2)
            %now looking at the overlap between the particular lymphocyte polygon and
            %tumour buffer
            tmpb = intersect(tumour_buffer_in2{i}{j}, lymphocyte_polygon{i}{l});
            %                 if tmpb.NumRegions == 0
            %                     %if no actual overlap at all then move on
            %                     continue
            %                 end
            tbuffer_l_intersection{i}{j} = union(tbuffer_l_intersection{i}{j}, tmpb); %a collection of polygons of t_l overlap
            [x,y] = centroid(lymphocyte_polygon{i}{l}); %a collection of centroid of the lymphoid polygons that have overlaps with this tumour cluster
            lymph_centroid_buffer{i}{j} = [lymph_centroid_buffer{i}{j}; [x,y]];
        end
        
        if isempty(tbuffer_l_intersection{i}{j})
            tbuffer_l_intersection_lymph_count{i}{j} = 0;
            tbuffer_l_intersection_area{i}{j} = 0;
            continue
        end
        
        tbuffer_l_intersection_lymph_count_logi{i}{j} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), tbuffer_l_intersection{i}{j}.Vertices(:,1), tbuffer_l_intersection{i}{j}.Vertices(:,2));
        tbuffer_l_intersection_lymph_count{i}{j} = sum(tbuffer_l_intersection_lymph_count_logi{i}{j});
        tbuffer_l_intersection_area{i}{j} = area(tbuffer_l_intersection{i}{j});
        
        if isempty(lymph_centroid_buffer{i}{j})
            lymph_centroid_buffer_logi{i}{j} = 0;
            tbuffer_l_centroid_count{i}{j} = 0;
            continue
        end
        lymph_centroid_buffer_logi{i}{j} = inpolygon(lymph_centroid_buffer{i}{j}(:,1), lymph_centroid_buffer{i}{j}(:,2), tumour_buffer_in2{i}{j}.Vertices(:,1), tumour_buffer_in2{i}{j}.Vertices(:,2));
        tbuffer_l_centroid_count{i}{j} = sum(lymph_centroid_buffer_logi{i}{j});
        
    end
end

%converting this to csv compatible format
t_l_intersection_lymph_count_cell = extractmycells(t_l_intersection_lymph_count);
t_l_intersection_area_cell = extractmycells(t_l_intersection_area);
t_l_centroid_count_cell = extractmycells(t_l_centroid_count);

save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_t_l_intersection_' num2str(lymph_cluster_size) '_' num2str(buffer_size) '.mat'], 't_l_intersection');
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_intersection_count' '_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], t_l_intersection_lymph_count_cell);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersection_area_' num2str(lymph_cluster_size) '_' num2str(buffer_size) '.csv'], t_l_intersection_area_cell);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_centroid_in_count' '_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], t_l_centroid_count_cell);


%converting this to csv compatible format
tbuffer_l_intersection_lymph_count_cell = extractmycells(tbuffer_l_intersection_lymph_count);
tbuffer_l_intersection_area_cell = extractmycells(tbuffer_l_intersection_area);
tbuffer_l_centroid_count_cell = extractmycells(tbuffer_l_centroid_count);

save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tbuffer_l_intersection_' num2str(lymph_cluster_size) '_' num2str(buffer_size) '.mat'], 'tbuffer_l_intersection');
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphbuffer_intersection_count' '_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tbuffer_l_intersection_lymph_count_cell);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersectionbuffer_area_' num2str(lymph_cluster_size) '_' num2str(buffer_size) '.csv'], tbuffer_l_intersection_area_cell);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_centroidbuffer_in_count' '_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tbuffer_l_centroid_count_cell);



%%end % for debug



%%
%calculate area

%core area
for i = 1:size(core_list, 2)
    core_area{i} = area(core_polygon{i});
end


for i = 1:size(tumour_polygon_in2, 2)
    %tumour polygon area
    if isempty(tumour_polygon_in2{i})
        tumour_polygon_area{i} = [];
        tumour_buffer_area{i} = [];
        continue
    end
    for j = 1:size(tumour_polygon_in2{i}, 2)
        if isempty(tumour_polygon_in2{i}{j})
            tumour_polygon_area{i}{j} = [];
            tumour_buffer_area{i}{j} = [];
            continue
        end
        tumour_polygon_area{i}{j} = area(tumour_polygon_in2{i}{j});
        tumour_buffer_area{i}{j} = area(tumour_buffer_in2{i}{j});
    end
end

%lymphocyte cluster area
for i = 1:size(lymphocyte_polygon, 2)
    if isempty(lymphocyte_polygon{i})
        lymphocyte_polygon_area{i} = [];
        continue
    end
    for l = 1:size(lymphocyte_polygon{i}, 2)
        lymphocyte_polygon_area{i}{l} = area(lymphocyte_polygon{i}{l});
    end
end

%create an output file
core_area_combined = extractmycells(core_area);
tumour_in_area_combined=extractmycells(tumour_polygon_area);
tumour_buffer_area_combined=extractmycells(tumour_buffer_area);
lymphocyte_polygon_area_combined=extractmycells(lymphocyte_polygon_area);


%turn into csv files
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_area_in_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumour_in_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_area_buffer_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumour_buffer_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_core_area_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], core_area_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_area_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], lymphocyte_polygon_area_combined);

%%
%get the lymphocyte counts as well as tumour cell count
for i = 1:size(tumour_polygon_in2, 2)
    if isempty(tumour_polygon_in2{i})
        tumourcell_in_count{i} = [];
        tumourcell_buffer_count{i}=[];
        continue
    end
    
    for j = 1:size(tumour_polygon_in2{i}, 2)
        if isempty(tumour_polygon_in2{i}{j})
            tumourcell_in_count{i}{j} = [];
            tumourcell_buffer_count{i}{j}=[];
            lymph_in_count{i}{j} = [];
            continue
        end
        %now we calculate the number of lymphocyte within tumour polygon
        lymph_in{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), tumour_polygon_in2{i}{j}.Vertices(:,1), tumour_polygon_in2{i}{j}.Vertices(:,2));
        lymph_in_count{i}{j} = sum(lymph_in{i}{j});
        %now we calculate the number of tumour cells within tumour polygon
        tumourcell_in{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), tumour_polygon_in2{i}{j}.Vertices(:,1), tumour_polygon_in2{i}{j}.Vertices(:,2));
        tumourcell_in_count{i}{j} = sum(tumourcell_in{i}{j});
        % number of lymphocyte in buffer
        lymph_buffer{i}{j} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), tumour_buffer_in2{i}{j}.Vertices(:,1), tumour_buffer_in2{i}{j}.Vertices(:,2));
        lymph_buffer_count{i}{j} = sum(lymph_buffer{i}{j});
        %number of tumour cells in buffer
        tumourcell_buffer{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), tumour_buffer_in2{i}{j}.Vertices(:,1), tumour_buffer_in2{i}{j}.Vertices(:,2));
        tumourcell_buffer_count{i}{j} = sum(tumourcell_buffer{i}{j});
    end
end

for i = 1:size(lymphocyte_polygon, 2)
    if isempty(lymphocyte_polygon{i})
        lymph_lymphcluster_count{i} = [];
        tumour_lymphcluster_count{i} = [];
        continue
    end
    
    for l = 1:size(lymphocyte_polygon{i},2)
        if isempty(lymphocyte_polygon{i}{l})
            lymph_lymphcluster_count{i}{l} = [];
            tumour_lymphcluster_count{i}{l} = [];
            continue
        end
        lymph_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
        lymph_lymphcluster_count{i}{l} = sum(lymph_lymphcluster{i}{l});
        tumour_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
        tumour_lymphcluster_count{i}{l} = sum(tumour_lymphcluster{i}{l});
    end
    
end

%create an output file
lymph_in_count_combined=extractmycells(lymph_in_count);
lymph_buffer_count_combined=extractmycells(lymph_buffer_count);
tumourcell_in_count_combined=extractmycells(tumourcell_in_count);
tumourcell_buffer_count_combined=extractmycells(tumourcell_buffer_count);
lymph_lymphcluster_count_combined = extractmycells(lymph_lymphcluster_count);
tumour_lymphcluster_count_combined = extractmycells(tumour_lymphcluster_count);

%turn into csv files
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_in_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], lymph_in_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_buffer_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], lymph_buffer_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumourcell_in_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumourcell_in_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumourcell_buffer_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumourcell_buffer_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], lymph_lymphcluster_count_combined);
csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_lymphcluster_count_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.csv'], tumour_lymphcluster_count_combined);


%% clear workspace before next image when using the for loop

% % save('image_list.mat', 'image_list')
% % clear all
% % load('image_list')


