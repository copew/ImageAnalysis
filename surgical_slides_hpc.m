% This script works on the slides with invasive tumour

% an cellularity calculation is attempted
% also interaction between tumour cells, lymphocytes and normal cells

image_list = [600227]; 
% done: 602185,599871, 600227, 600245, 600268,600281 
% roo big: 605745, 
image_path_stem = '/Volumes/Wei/Neoadjuvant/trimmed_matlab_data/';


% function tum_lymph_cluster_hpc(image_filenumber_fullpath)
%
% %loading data
% disp(image_filenumber_fullpath)
% disp(class(image_filenumber_fullpath))
% if ischar(image_filenumber_fullpath)
% [image_path_stem,image_filenumber,~] = fileparts(image_filenumber_fullpath)
% image_filenumber = str2num(image_filenumber);
% if isempty(image_path_stem)
% image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt'
% warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
% end
% elseif isnumeric(image_filenumber_fullpath)
% image_filenumber = image_filenumber_fullpath;
% image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt'
% warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
% end

for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    
    % set parameters
    cluster_size = [5];
    lymph_cluster_size = [20];
    lymph_cutoff_size = [50];
    
    % make a folder to save files
    mkdir(num2str(image_filenumber));
    
    data= load([image_path_stem '/' num2str(image_filenumber) '.mat']);
    image_path = ['/Volumes/Wei/Neoadjuvant/NeoAdjuvant_Images/TransNeoSlides/' num2str(image_filenumber) '.svs'];
    
    %create indexing
    X_ind = 1;
    Y_ind = 2;
    cell_ind = 3;
    
    % get data into correct format
    data_trimmed = data.data_mini;
    data_trimmed = num2cell(data_trimmed, 1); %so that the cell index works later on
    
    
    
    %% this bit is for debug for plotting
    % First get the thumbnail
    image_info=imfinfo(image_path);
    thumbnail_height_scale_factor = image_info(1).Height/image_info(2).Height;
    thumbnail_width_scale_factor = image_info(1).Width/image_info(2).Width;
    thumbnail_overall_scale_factor = mean([thumbnail_height_scale_factor,thumbnail_width_scale_factor]);
    low_res_layer = length(image_info)-2;
    high_res_layer = 1;
    thumbnail_layer = 2;
    %By convention:
    % First level	Full resolution image
    % Second level	Thumbnail
    % Third level to N-2 Level	A reduction by a power of 2 (4:1 ratio, 16:1 ratio, 32:1 ratio, etc)
    % N-1 Level	Slide Label
    % N Level	Entire Slide with cropped region delineated in green
    %image_io=imread(image_path,'Index',low_res_layer);
    
    thumbnail_io=imread(image_path,'Index',thumbnail_layer);
    large_thumbnail_io = imresize(thumbnail_io,thumbnail_overall_scale_factor);
    
    %% identify tumour clusters
    
    base_cells = data_trimmed{cell_ind}==1;
    neighbour_cells = data_trimmed{cell_ind}==1;

    
    [all_multi_real_distances, all_indexes] = pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
    
    i = 0;
    if size(all_multi_real_distances,1)<max(cluster_size)+1
        all_multi_real_distances(size(all_multi_real_distances,1)+1:max(cluster_size)+1,:)=NaN;
        this_tumour_cluster_boundary= [];
        continue;
    end
    
    %now use DBScan
    for this_clustsize = cluster_size
        i = i+1;
        C=0;
        
        n=sum(neighbour_cells);
        cluster_membership{i}=zeros(n,1);
        
        visited=false(n,1);
        isnoise=false(n,1);
        
        %Calculate epsilon from the knee point of a
        %1000-bin histogram to the k-1th nearest neighbour
        %Remember that first index is self, so no need
        %to subtract one.
        [hist_y, hist_x] = hist(all_multi_real_distances(this_clustsize,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        
        for this_cell=1:n
            if ~visited(this_cell)
                visited(this_cell)=true;
                Neighbours=setdiff(all_indexes(all_multi_real_distances(:,this_cell)<=epsilon/2,this_cell),this_cell)';
                if numel(Neighbours)<this_clustsize
                    % X(i,:) is NOISE
                    isnoise(this_cell)=true;
                else
                    C=C+1;
                    cluster_membership{i}(this_cell)=C;
                    k = 1;
                    while true
                        this_neighbour_cell = Neighbours(k);
                        
                        if ~visited(this_neighbour_cell)
                            visited(this_neighbour_cell)=true;
                            Neighbours2=setdiff(all_indexes(all_multi_real_distances(:,this_neighbour_cell)<=epsilon/2,this_neighbour_cell),this_neighbour_cell)';
                            if numel(Neighbours2)>=this_clustsize
                                Neighbours=[Neighbours Neighbours2];   %#ok
                            end
                        end
                        if cluster_membership{i}(this_neighbour_cell)==0
                            cluster_membership{i}(this_neighbour_cell)=C;
                        end
                        
                        k = k + 1;
                        if k > numel(Neighbours)
                            break;
                        end
                    end
                end
            end
        end
        
        x_data_here = data_trimmed{X_ind}(neighbour_cells);
        y_data_here = data_trimmed{Y_ind}(neighbour_cells);
        
        for this_cluster_number = 1:max(cluster_membership{i})
            this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
            this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
            this_tumour_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
            this_tumour_cluster_boundary{i}{this_cluster_number} = [this_x_subset(this_tumour_cluster_boundary_numbers),this_y_subset(this_tumour_cluster_boundary_numbers)];
        end
    end
    % the output that feed into the next step is this_tumour_cluster_boundary
  
    %% lymphocyte clustering
  
    % lymphocyte clusters
    
    base_cells = data_trimmed{cell_ind}==2;
    neighbour_cells = data_trimmed{cell_ind}==2;
    
    %ignore the really small chunks that has very few tumour cells
%     if size(in_core{this_core}{cell_ind}==2,1) < 10
%         continue
%     end
%     
    [all_multi_real_distances, all_indexes] = pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest', lymph_cluster_size+1);
    
    i = 0;
    if size(all_multi_real_distances,1)<max(lymph_cluster_size)+1
        all_multi_real_distances(size(all_multi_real_distances,1)+1:max(lymph_cluster_size)+1,:)=NaN;
        this_lymphocyte_cluster_boundary = [];
        continue;
    end
    
    %now use DBScan
    for this_clustsize = lymph_cluster_size
        i = i+1;
        C=0;
        
        n=sum(neighbour_cells);
        cluster_membership{i}=zeros(n,1);
        
        visited=false(n,1);
        isnoise=false(n,1);
        
        %Calculate epsilon from the knee point of a
        %1000-bin histogram to the k-1th nearest neighbour
        %Remember that first index is self, so no need
        %to subtract one.
        %[hist_y, hist_x] = hist(all_multi_real_distances{this_core}(this_clustsize,:),1000);
        [hist_y, hist_x] = hist(all_multi_real_distances(3,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        
        for this_cell=1:n
            if ~visited(this_cell)
                visited(this_cell)=true;
                Neighbours=setdiff(all_indexes(all_multi_real_distances(:,this_cell)<=epsilon/2,this_cell),this_cell)';
                if numel(Neighbours)<this_clustsize
                    % X(i,:) is NOISE
                    isnoise(this_cell)=true;
                else
                    C=C+1;
                    cluster_membership{i}(this_cell)=C;
                    k = 1;
                    while true
                        this_neighbour_cell = Neighbours(k);
                        
                        if ~visited(this_neighbour_cell)
                            visited(this_neighbour_cell)=true;
                            Neighbours2=setdiff(all_indexes(all_multi_real_distances(:,this_neighbour_cell)<=epsilon/2,this_neighbour_cell),this_neighbour_cell)';
                            if numel(Neighbours2)>=this_clustsize
                                Neighbours=[Neighbours Neighbours2];   %#ok
                            end
                        end
                        if cluster_membership{i}(this_neighbour_cell)==0
                            cluster_membership{i}(this_neighbour_cell)=C;
                        end
                        
                        k = k + 1;
                        if k > numel(Neighbours)
                            break;
                        end
                    end
                end
            end
        end
        
        x_data_here = data_trimmed{X_ind}(neighbour_cells);
        y_data_here = data_trimmed{Y_ind}(neighbour_cells);
        
        for this_cluster_number = 1:max(cluster_membership{i})
            this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
            this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
            this_lymphocyte_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
            this_lymphocyte_cluster_boundary{i}{this_cluster_number} = [this_x_subset(this_lymphocyte_cluster_boundary_numbers),this_y_subset(this_lymphocyte_cluster_boundary_numbers)];
        end
    end

    



% %% this bit is to match the dimension of the cluster boundaries with the core_list, as we
% % may need empty cells at the end of the list for cores without either tumour clusters or
% % lymphocyte clusters
% 
% % firstly work out the number of cores in total (max of core_total), and how many are in
% % each list. then work out the difference and append.
% tumour_core_total = size(this_tumour_cluster_boundary, 2);
% lymph_core_total = size(this_lymphocyte_cluster_boundary, 2);
% core_total = max(core_list);
% 
% % if same sizes then dont' worry about it....
% if core_total == tumour_core_total
%     %do nothing
% else
%     tumour_diff = core_total - tumour_core_total;
%     for diff = 1:tumour_diff
%         this_tumour_cluster_boundary{1, tumour_core_total + diff} = [];
%     end
% end
% 
% if core_total == lymph_core_total
%     %do nothing
% else
%     lymph_diff = core_total - lymph_core_total;
%     for diff = 1:lymph_diff
%         this_lymphocyte_cluster_boundary{1, lymph_core_total + diff} = [];
%     end
% end
% 
% 
% save(['./' num2str(image_filenumber) '/workspace_clusters.mat']);
% 
%% this is to check if the clusters are correct
figure
ax=gca();
hold on;
imshow(large_thumbnail_io);
hold(ax, 'on');
% draw the tumour cluster first
    for i = 1:size(this_tumour_cluster_boundary, 2)

        if isempty(this_tumour_cluster_boundary{i}) == 1
            continue
        end

        % this_tumour_cluster_boundary{i}{1} =
        % this_tumour_cluster_boundary{i}{1}(~cellfun('isempty',
        % this_tumour_cluster_boundary{i}{1})); % this may be introducing extra layers
        % unnecessarily
        for j = 1:size(this_tumour_cluster_boundary{i}, 2)
            hold on;
            if isempty(this_tumour_cluster_boundary{i}{j}) == 1
                continue
            end
            plot(this_tumour_cluster_boundary{i}{j}(:,1), this_tumour_cluster_boundary{i}{j}(:,2), 'k-');
        end
        hold on;
    end

% now draw the lymphocyte cluster
for i = 1:size(this_lymphocyte_cluster_boundary, 2)
    if isempty(this_lymphocyte_cluster_boundary{i}) == 1
        continue
    end

    %this_lymphocyte_cluster_boundary{i}{1} =
    %this_lymphocyte_cluster_boundary{i}{1}(~cellfun('isempty',
    %this_lymphocyte_cluster_boundary{i}{1})); % this may be introducing extralayers
    %unnecesarily
    for j = 1:size(this_lymphocyte_cluster_boundary{i}, 2)
        hold on;
        if isempty(this_lymphocyte_cluster_boundary{i}{j}) == 1
            continue
        end
        plot(this_lymphocyte_cluster_boundary{i}{j}(:,1), this_lymphocyte_cluster_boundary{i}{j}(:,2), 'r-');
    end
    hold on;
end
hold(ax, 'off')

end


%% now convert each into a polygon - old code
% %the this_tumour/lymphocyte_cluster_boundary is organised in such a way that it has core number
% %first, then some times it has extra layers before reaches the collection of number for
% %the clusters. so will unnest within each core
for i =1:size(core_list, 2) %First un-nest the cells
    if isempty(this_tumour_cluster_boundary{core_list(i)}) == 1;
        continue
    end
    this_tumour_cluster_boundary{core_list(i)} = extractmycells(this_tumour_cluster_boundary{core_list(i)});
end

for i =1:size(core_list, 2)
    if isempty(this_lymphocyte_cluster_boundary{core_list(i)}) ==1;
        continue
    end
    this_lymphocyte_cluster_boundary{core_list(i)} = extractmycells(this_lymphocyte_cluster_boundary{core_list(i)});
end


for i =1:size(core_list, 2)
    if isempty(this_tumour_cluster_boundary{core_list(i)}) == 1;
        continue
    end
    
    for j = 1:size(this_tumour_cluster_boundary{core_list(i)}, 2)
        if isempty(this_tumour_cluster_boundary{core_list(i)}{j}) == 1
            continue
        end
        tumour_polygon{core_list(i)}{j} = polyshape (this_tumour_cluster_boundary{core_list(i)}{j});
        %plot(tumour_polygon{core_list(i)}{j});
    end
end

for i =1:size(core_list, 2)
    if isempty(this_lymphocyte_cluster_boundary{core_list(i)}) == 1;
        continue
    end
    for k = 1:size(this_lymphocyte_cluster_boundary{core_list(i)}, 2)
        if isempty(this_lymphocyte_cluster_boundary{core_list(i)}{k}) == 1
            continue
        end
        lymphocyte_polygon{core_list(i)}{k} = polyshape (this_lymphocyte_cluster_boundary{core_list(i)}{k});
        %plot(tumour_polygon{core_list(i)}{j});
    end
end


% % this bit of the code is
% % %the this_tumour/lymphocyte_cluster_boundary is organised in such a way that it has core number
% % %first, then some times it has extra layers before reaches the collection of number for
% % %the clusters. so will unnest within each core
% for i =1:size(core_list, 2) %First un-nest the cells
%     this_tumour_cluster_boundary{core_list(i)} = extractmycells(this_tumour_cluster_boundary{core_list(i)});
%     this_lymphocyte_cluster_boundary{core_list(i)} = extractmycells(this_lymphocyte_cluster_boundary{core_list(i)});
% end
% % %Then get rid of empty cells from the array. Note that this then makes the
% % %core_list not suitable for indexing through. Can't remove items from the
% % %core list because there might be some tumours empty with lymphocytes
% % %present etc.,
total_core_polygon = 1:size(this_tumour_cluster_boundary, 2);

tumour_core_list = intersect(core_list, total_core_polygon(~cellfun('isempty',this_tumour_cluster_boundary))); %In case you need this later
% %this_tumour_cluster_boundary = this_tumour_cluster_boundary(~cellfun('isempty',this_tumour_cluster_boundary));
lymphocyte_core_list = intersect(core_list, total_core_polygon(~cellfun('isempty',this_lymphocyte_cluster_boundary))); %In case you need this later
% %this_lymphocyte_cluster_boundary = this_lymphocyte_cluster_boundary(~cellfun('isempty',this_lymphocyte_cluster_boundary));
%
% for i =1:size(this_tumour_cluster_boundary,2)
%     for j = 1:size(this_tumour_cluster_boundary{i}, 2)
%         try
%             warning('off','all')
%             tumour_polygon{i}{j} = polyshape(this_tumour_cluster_boundary{i}{j});
%             warning('on','all')
%         catch
%             warning('on','all')
%             warning(['Something went wrong with tumour polygons for core ' num2str(tumour_core_list(i)) ' moving on, the array is:'])
%             disp(this_tumour_cluster_boundary{i}{j})
%             continue
%         end
%         %plot(tumour_polygon{i}{j});
%     end
% end
%
% for i =1:size(this_lymphocyte_cluster_boundary,2)
%     for j = 1:size(this_lymphocyte_cluster_boundary{i}, 2)
%         try
%             warning('off','all')
%             lymphocyte_polygon{i}{j} = polyshape(this_lymphocyte_cluster_boundary{i}{j});
%             warning('on','all')
%         catch
%             warning('on','all')
%             warning(['Something went wrong with lymphocyte polygons for core ' num2str(lymphocyte_core_list(i)) ' moving on, the array is:'])
%             disp(this_lymphocyte_cluster_boundary{i}{j})
%             continue
%         end
%         %plot(lymphocyte_polygon{i}{j});
%     end
% end

%the end result is tumour_polygon, and lymphocyte_polygon

%% sort out overlaps


%need to work out tumour polygon if it overlaps.
%but the problem occurs at what to do with the old ones...
overlap = cell(0);
for i =1:size(tumour_polygon, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        if isempty(tumour_polygon{i}{j})
            continue
        end
        for k = (j+1):size(tumour_polygon{i}, 2)
            if isempty(tumour_polygon{i}{k})
                continue
            end
            overlap{i}{j}{k} = inpolygon(tumour_polygon{i}{k}.Vertices(:,1), tumour_polygon{i}{k}.Vertices(:,2), tumour_polygon{i}{j}.Vertices(:,1), tumour_polygon{i}{j}.Vertices(:,2));
            if sum(overlap{i}{j}{k})>0
                tumour_polygon{i}{j} = union(tumour_polygon{i}{j}, tumour_polygon{i}{k});
                tumour_polygon{i}{k} = [];
            else
                continue;
            end
        end
    end
    
end

%now remove empty polygons

% none of these work.... but will just add continue if empty for later on

%     find(~cellfun('isempty', x))
%
%     for i =1:size(core_list, 2)
%         for j = 1:size(tumour_polygon{core_list(i)}, 2)
%             these_empty{i}{j} = find(~cellfun('isempty', tumour_polygon{core_list(i)}{j}(:)));
%         end
%     end
%
%     tumour_polygon_non_empty= tumour_polygon;
%     tumour_polygon_non_empty{core_list(i)}(these_empty) = [];


%this is an atteempt with while loop, but not very sensible.
%     tumour_polygon_non_empty= tumour_polygon;
%     for i =1:size(core_list, 2)
%         j=0;
%         while j <size(tumour_polygon{core_list(i)}, 2)
%             j=j+1;
%             if isempty(tumour_polygon{core_list(i)}{j})
%                 tumour_polygon_non_empty{core_list(i)}(j) = [];
%                 j=0;
%             end
%         end
%     end



%%
%this bit then works out the buffer bit

tumour_buffer = cell(0);
for i = 1:size(tumour_polygon, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        if isempty(tumour_polygon{i}{j})
            continue
        end
        tumour_buffer{i}{j} = polybuffer(tumour_polygon{i}{j},100);
        %tumour_buffer{core_list(i)}{j} = polyshape(tumour_buffer{core_list(i)}{j});
        %this is not necessary as it is already a polyshape
    end
end

%sanity check :)
%     figure
%     ax=gca();
%     hold on;
%     imshow(large_thumbnail_io);
%     hold(ax, 'on');
%     for i = 1:size(core_list, 2)
%         %this_cluster_boundary{core_list(i)}{1} = this_cluster_boundary{core_list(i)}{1}(~cellfun('isempty', this_cluster_boundary{core_list(i)}{1}));
%         for j = 1:size(tumour_buffer{core_list(i)}, 2)
%             plot(tumour_polygon{core_list(i)}{j});
%             hold on;
%             plot(tumour_buffer{core_list(i)}{j})
%         end
%         hold on;
%     end
%     hold(ax, 'off')

%% exclude the area that falls outside the cores
%intersecting regions so only within core area is calculated
%
for i = 1:size(tumour_core_list, 2)
    %     for j = 1:size(tumour_polygon{i}, 2)
    %         if isempty(tumour_polygon{i}{j})
    %             continue
    %         end
    %         tumour_polygon_in{i}{j} = intersect(tumour_polygon{i}{j}, core_polygon{tumour_core_list(i)});
    %     end
    % the above is unnecessary as only the buffer goes outside
    if isempty(tumour_buffer{tumour_core_list(i)})
        continue
    end
    for k = 1:size(tumour_buffer{tumour_core_list(i)},2)
        if isempty(tumour_buffer{tumour_core_list(i)}{k})
            continue
        end
        tumour_buffer_in{tumour_core_list(i)}{k} = intersect(tumour_buffer{tumour_core_list(i)}{k}, core_polygon{tumour_core_list(i)});
    end
end

%thi sis just for the rest of the code
tumour_polygon_in = tumour_polygon;

%%
%exclude the overlapped regions between the tumour and buffer
%if tumour and buffer overlap: treat as tumour
%if buffer and buffer overlap: just carry on as normal


for i = 1:size(tumour_polygon, 2)
    for j = 1:size(tumour_polygon_in{i}, 2)
        if isempty(tumour_polygon_in{i}{j})
            continue
        end
        for k = 1:size(tumour_buffer{i},2)
            if isempty(tumour_buffer{i}{k})
                continue
            end
            tumour_buffer_in{i}{k} = subtract(tumour_buffer_in{i}{k}, tumour_polygon_in{i}{j});
        end
    end
end

%% now save the files
save(['./' num2str(image_filenumber) '/tumour_polygon_in.mat'], 'tumour_polygon_in');
save(['./' num2str(image_filenumber) '/tumour_buffer_in.mat'], 'tumour_buffer_in');
save(['./' num2str(image_filenumber) '/lymphocyte_polygon.mat'], 'lymphocyte_polygon');

%% looking at spatial relationship between tumour clusters and lymphocyte clusters

% first of all look at overlaps
% need to select the cores with both tumour polygons and lymphocyte polygons
overlap_list = intersect(tumour_core_list, lymphocyte_core_list);
t_l_intersection = [];

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
            tmp_lymph = inpolygon(in_core{overlap_list(i)}{X_ind}(in_core{overlap_list(i)}{cell_ind}==2), in_core{overlap_list(i)}{Y_ind}(in_core{overlap_list(i)}{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
            tmp_lymph_count = sum(tmp_lymph);
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
for i = 1:size(tumour_core_list, 2)
    %core area
    core_area{tumour_core_list(i)} = area(core_polygon{i});
    %tumour polygon area
    for j = 1:size(tumour_polygon_in{tumour_core_list(i)}, 2)
        tumour_polygon_area{tumour_core_list(i)}{j} = area(tumour_polygon_in{tumour_core_list(i)}{j});
    end
    %tumour buffer zone area
    for k = 1:size(tumour_buffer_in{tumour_core_list(i)},2)
        tumour_buffer_area{tumour_core_list(i)}{k} = area(tumour_buffer_in{tumour_core_list(i)}{k});
    end
end

%lymphocyte cluster area
for i = 1:size(lymphocyte_core_list, 2)
    for l = 1:size(lymphocyte_polygon{lymphocyte_core_list(i)}, 2)
        lymphocyte_polygon_area{lymphocyte_core_list(i)}{l} = area(lymphocyte_polygon{lymphocyte_core_list(i)}{l});
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
for i = 1:size(tumour_core_list, 2)
    for j = 1:size(tumour_polygon_in{tumour_core_list(i)}, 2)
        if isempty(tumour_polygon_in{tumour_core_list(i)}{j})
            continue
        end
        %now we calculate the number of lymphocyte within tumour polygon
        lymph_in{tumour_core_list(i)}{j}=inpolygon(in_core{tumour_core_list(i)}{X_ind}(in_core{tumour_core_list(i)}{cell_ind}==2), in_core{tumour_core_list(i)}{Y_ind}(in_core{tumour_core_list(i)}{cell_ind}==2), tumour_polygon_in{tumour_core_list(i)}{j}.Vertices(:,1), tumour_polygon_in{tumour_core_list(i)}{j}.Vertices(:,2));
        lymph_in_count{tumour_core_list(i)}{j} = sum(lymph_in{tumour_core_list(i)}{j});
        %now we calculate the number of tumour cells within tumour polygon
        tumourcell_in{tumour_core_list(i)}{j}=inpolygon(in_core{tumour_core_list(i)}{X_ind}(in_core{tumour_core_list(i)}{cell_ind}==1), in_core{tumour_core_list(i)}{Y_ind}(in_core{tumour_core_list(i)}{cell_ind}==1), tumour_polygon_in{tumour_core_list(i)}{j}.Vertices(:,1), tumour_polygon_in{tumour_core_list(i)}{j}.Vertices(:,2));
        tumourcell_in_count{tumour_core_list(i)}{j} = sum(tumourcell_in{tumour_core_list(i)}{j});
    end
    
    for k = 1:size(tumour_buffer_in{tumour_core_list(i)},2)
        if isempty(tumour_buffer_in{tumour_core_list(i)}{k})
            continue
        end
        lymph_buffer{tumour_core_list(i)}{k} = inpolygon(in_core{tumour_core_list(i)}{X_ind}(in_core{tumour_core_list(i)}{cell_ind}==2), in_core{tumour_core_list(i)}{Y_ind}(in_core{tumour_core_list(i)}{cell_ind}==2), tumour_buffer_in{tumour_core_list(i)}{k}.Vertices(:,1), tumour_buffer_in{tumour_core_list(i)}{k}.Vertices(:,2));
        lymph_buffer_count{tumour_core_list(i)}{k} = sum(lymph_buffer{tumour_core_list(i)}{k});
        tumourcell_buffer{tumour_core_list(i)}{k}=inpolygon(in_core{tumour_core_list(i)}{X_ind}(in_core{tumour_core_list(i)}{cell_ind}==1), in_core{tumour_core_list(i)}{Y_ind}(in_core{tumour_core_list(i)}{cell_ind}==1), tumour_buffer_in{tumour_core_list(i)}{k}.Vertices(:,1), tumour_buffer_in{tumour_core_list(i)}{k}.Vertices(:,2));
        tumourcell_buffer_count{tumour_core_list(i)}{k} = sum(tumourcell_buffer{tumour_core_list(i)}{k});
    end
end

for i = 1:size(lymphocyte_core_list, 2)
    for l = 1:size(lymphocyte_polygon{lymphocyte_core_list(i)},2)
        if isempty(lymphocyte_polygon{lymphocyte_core_list(i)}{l})
            continue
        end
        lymph_lymphcluster{lymphocyte_core_list(i)}{l} = inpolygon(in_core{lymphocyte_core_list(i)}{X_ind}(in_core{lymphocyte_core_list(i)}{cell_ind}==2), in_core{lymphocyte_core_list(i)}{Y_ind}(in_core{lymphocyte_core_list(i)}{cell_ind}==2), lymphocyte_polygon{lymphocyte_core_list(i)}{l}.Vertices(:,1), lymphocyte_polygon{lymphocyte_core_list(i)}{l}.Vertices(:,2));
        lymph_lymphcluster_count{lymphocyte_core_list(i)}{l} = sum(lymph_lymphcluster{lymphocyte_core_list(i)}{l});
        tumour_lymphcluster{lymphocyte_core_list(i)}{l} = inpolygon(in_core{lymphocyte_core_list(i)}{X_ind}(in_core{lymphocyte_core_list(i)}{cell_ind}==1), in_core{lymphocyte_core_list(i)}{Y_ind}(in_core{lymphocyte_core_list(i)}{cell_ind}==1), lymphocyte_polygon{lymphocyte_core_list(i)}{l}.Vertices(:,1), lymphocyte_polygon{lymphocyte_core_list(i)}{l}.Vertices(:,2));
        tumour_lymphcluster_count{lymphocyte_core_list(i)}{l} = sum(lymph_lymphcluster{lymphocyte_core_list(i)}{l});
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

%% clear workspace before next image when using the for loop

% save(image_list, 'image_list')
% clear all
% load('image_list')

