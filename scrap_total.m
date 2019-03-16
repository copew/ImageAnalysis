% This script requires an image filek (.svs) and its matching cell
% segmentation file (.fits).

image_filenumber = 619857;

cluster_size = [10];


% loading files and trim data
data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];

X_ind = 3;
Y_ind = 4;
cell_ind = size(data,2);

data_trimmed = data;
for i = 1:size(data, 2)
    data_trimmed{i} = data_trimmed{i}(data{cell_ind}~=0);
end

% and combine tumour and normal
data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;

num_total = size(data_trimmed{cell_ind},1);
num_tum_cells = sum(data_trimmed{cell_ind}==1);
num_ly_cells = sum(data_trimmed{cell_ind}==2);
num_str_cells = sum(data_trimmed{cell_ind}==3);

% Then compute the proportion of each cell type
prop_tum_cells = num_tum_cells/num_total;
prop_ly_cells = num_ly_cells/num_total;
prop_str_cells = num_str_cells/num_total;


%% Create a boundary around each core
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

%now try to convert to bw and then create boundary
this_image = imbinarize(large_thumbnail_io); % Binarize the image on all three colour layers
this_image_2 = any(~this_image,3); % Select positive pixels in any colour layer
this_image_expanded = bwdist(this_image_2) <= 100; % Expand the image to allow 'almost connected' cells
cc = bwconncomp(this_image_expanded,4); % Now work out connected clusters
big_cluster = false(1,cc.NumObjects); % This bit rejects 'crud' outside of large areas
for i = 1:cc.NumObjects
    big_cluster(i) = size(cc.PixelIdxList{i},1) > 1000000; % Absolute value for quick first pass, need to improve this for applicability XXXFIXME
end
grain = cell(0);
for i = 1:cc.NumObjects
    if big_cluster(i)
        grain{end+1} = false(size(this_image_2));
        grain{end}(cc.PixelIdxList{i}) = true; % Now create logical images of each core
    end
end
this_boundary = cell(size(grain));
for i = 1:size(grain,2)
    clear row col % Now for each core work out a starting point for boundary determination
    for j = 1:size(grain{i},2)
        row = min(find(grain{i}(:,j)));
        if row
            break
        end
    end
    col = j;
    this_boundary{i} = bwtraceboundary(grain{i},[row col],'S'); % Trace the boundary
end
% figure % Plot the results for sanity
% imshow(large_thumbnail_io)
% hold on;
% for i = 1:size(grain,2)
%     plot(this_boundary{i}(:,2),this_boundary{i}(:,1),'g','LineWidth',1);
% end

% Now convert each into a polygon
core_polygon= cell(0);
for i = 1:size(this_boundary, 2)
    core_polygon{i} = polyshape (this_boundary{i}(:,2),this_boundary{i}(:,1)); %note the doordinate is this way round
    %plot(core_polygon{i}) %in case the coordinates flipped
end



% for i = 1:8
%     plot(core_polygon{i});
%     plot(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), 'k-');
%     pause(3);
% end

% the useful output of this chunck is the core_polygon cell array

%%
% now subset coordinates of cells
in_polygon=cell(0)
in_core=cell(0)

for i = 1:size(core_polygon, 2)
    in_polygon{i} = inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, core_polygon{i}.Vertices(:,1), core_polygon{i}.Vertices(:,2));
    
    for j = 1:size(data_trimmed, 2)
        in_core{i}{j} = data_trimmed{j}(in_polygon{i}, :);
    end
    %plot(in_core{i}{X_ind}, in_core{i}{Y_ind}, 'bd') %again in case the
    
    % cells are in the wrong place
    
end


%%
%Now we can perform cluster analysis and neighbourhood analysis core by
%core, the cells are in each in_core and is same format as data_trimmed


% Now compute the tumour clusters

for this_core = 1:size(core_polygon,2)

    base_cells = in_core{this_core}{cell_ind}==1;
    neighbour_cells = in_core{this_core}{cell_ind}==1;

    [all_multi_real_distances{this_core}, all_indexes{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);

    i = 0;
    for this_clustsize = cluster_size
        i = i+1;

        %Now calculate the clusters using DBScan if same cell type

        C=0;

        n=sum(neighbour_cells);
        cluster_membership{i}=zeros(n,1);

        visited=false(n,1);
        isnoise=false(n,1);

        %Calculate epsilon from the knee point of a
        %1000-bin histogram to the k-1th nearest neighbour
        %Remember that first index is self, so no need
        %to subtract one.
        [hist_y, hist_x] = hist(all_multi_real_distances{this_core}(this_clustsize,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        rand_cell_order = randperm(n); %Randomise cell visitation order
        for this_cell=rand_cell_order
            if ~visited(this_cell)
                visited(this_cell)=true;
                Neighbours=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_cell)<=epsilon,this_cell),this_cell)';
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
                            Neighbours2=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_neighbour_cell)<=epsilon,this_neighbour_cell),this_neighbour_cell)';
                            if numel(Neighbours2)>=this_clustsize
                                Neighbours=[Neighbours Neighbours2];   %#ok
                            end
                        end
                        temp_cluster_membership{i}(this_neighbour_cell)=C; %New step for competitive cluster size comparison
                        k = k + 1;
                        if k > numel(Neighbours)
                            break;
                        end
                    end
                    for k2 = 1:(k-1) %Competitive cluster size comparison step
                        this_neighbour_cell = Neighbours(k2);
                        if cluster_membership{i}(this_neighbour_cell)==0
                            cluster_membership{i}(this_neighbour_cell)=C;
                        elseif sum(cluster_membership{1}(:) == cluster_membership{i}(this_neighbour_cell)) < k
                            cluster_membership{i}(this_neighbour_cell)=C;
                        end
                    end
                end

            end
        end
        x_data_here = in_core{this_core}{X_ind}(neighbour_cells);
        y_data_here = in_core{this_core}{Y_ind}(neighbour_cells);

        for this_cluster_number = 1:max(cluster_membership{i})
            this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
            this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
            this_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset, 1);
            this_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_cluster_boundary_numbers),this_y_subset(this_cluster_boundary_numbers)];
        end
    end
end

%this is tom's test fix.... sadly it didn't work.  so the above commented section is what
%it used to be.....
% for this_core = 1:size(core_polygon,2)
%     
%     base_cells = in_core{this_core}{cell_ind}==1;
%     neighbour_cells = in_core{this_core}{cell_ind}==1;
%     
%     [all_multi_real_distances{this_core}, all_indexes{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
%     
%     i = 0;
%     for this_clustsize = cluster_size
%         i = i+1;
%         
%         %Now calculate the clusters using DBScan if same cell type
%         
%         C=0;
%         
%         n=sum(neighbour_cells);
%         cluster_membership{i}=zeros(n,1);
%         
%         visited=false(n,1);
%         isnoise=false(n,1);
%         
%         %Calculate epsilon from the knee point of a
%         %1000-bin histogram to the k-1th nearest neighbour
%         %Remember that first index is self, so no need
%         %to subtract one.
%         [hist_y, hist_x] = hist(all_multi_real_distances{this_core}(this_clustsize,:),1000);
%         [epsilon, ~] = knee_pt(hist_y,hist_x,true);
%         %epsilon = epsilon*10; % FOR DEBUG ONLY
%         rand_cell_order = randperm(n);
%         for this_cell=rand_cell_order
%             if ~visited(this_cell)
%                 visited(this_cell)=true;
%                 Neighbours=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_cell)<=epsilon,this_cell),this_cell)';
%                 if numel(Neighbours)<this_clustsize
%                     % X(i,:) is NOISE
%                     isnoise(this_cell)=true;
%                 else
%                     Neighbours = Neighbours(randperm(length(Neighbours)));
%                     C=C+1;
%                     cluster_membership{i}(this_cell)=C;
%                     k = 1;
%                     while true
%                         this_neighbour_cell = Neighbours(k);
%                         
%                         if ~visited(this_neighbour_cell) % would be better to comment this out, but then takes unfeasibly long
%                             visited(this_neighbour_cell)=true;
%                             Neighbours2=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_neighbour_cell)<=epsilon,this_neighbour_cell),this_neighbour_cell)';
%                             %if numel(Neighbours2)>=this_clustsize
%                             Neighbours2 = setdiff(Neighbours2, Neighbours);
%                             Neighbours2 = Neighbours2(randperm(length(Neighbours2)));
%                             Neighbours=[Neighbours Neighbours2];   %#ok
%                             %end
%                         end
%                         temp_cluster_membership{i}(this_neighbour_cell)=C; %New step for competitive cluster size comparison
%                         k = k + 1;
%                         if k > numel(Neighbours)
%                             break;
%                         end
%                     end
%                     for k2 = 1:(k-1) %Competitive cluster size comparison step
%                         this_neighbour_cell = Neighbours(k2);
%                         if cluster_membership{i}(this_neighbour_cell)==0
%                             cluster_membership{i}(this_neighbour_cell)=C;
%                         elseif sum(cluster_membership{1}(:) == cluster_membership{i}(this_neighbour_cell)) < k
%                             cluster_membership{i}(this_neighbour_cell)=C;
%                         end
%                     end
%                 end
%                 
%             end
%         end
%         x_data_here = in_core{this_core}{X_ind}(neighbour_cells);
%         y_data_here = in_core{this_core}{Y_ind}(neighbour_cells);
%         
%         for this_cluster_number = 1:max(cluster_membership{i})
%             this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
%             this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
%             this_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset, 1);
%             this_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_cluster_boundary_numbers),this_y_subset(this_cluster_boundary_numbers)];
%         end
%     end
% end

%this is to check if the clusters are correct - and it is... phew
% figure
% ax=gca();
% hold on;
% imshow(large_thumbnail_io);
% hold(ax, 'on');
% for i = 1:size(this_cluster_boundary, 2)
%     this_cluster_boundary{i}{1} = this_cluster_boundary{i}{1}(~cellfun('isempty', this_cluster_boundary{i}{1}));
%     for j = 1:size(this_cluster_boundary{i}{1}, 2)
%      hold on;
%      plot(this_cluster_boundary{i}{1}{j}(:,1), this_cluster_boundary{i}{1}{j}(:,2), 'k-');
%     end
%     hold on
% end
%  hold(ax, 'off')


%now convert each into a polygon
for i = 1:size(this_cluster_boundary, 2)
    this_cluster_boundary{i}{1} = this_cluster_boundary{i}{1}(~cellfun('isempty', this_cluster_boundary{i}{1}));
    for j = 1:size(this_cluster_boundary{i}{1}, 2)
        tumour_polygon{i}{j} = polyshape (this_cluster_boundary{i}{1}{j});
        %plot(tumour_polygon{i}{j});
    end
end

%the end result is tumour_polygon

%%
%this bit then works out the buffer bit

tumour_buffer = cell(0);

for i = 1:size(this_cluster_boundary, 2)
    for j = 1:size(this_cluster_boundary{i}{1}, 2)
        tumour_buffer{i}{j} = polybuffer(tumour_polygon{i}{j},100);
    end
end

%sanity check :)
% figure
% ax=gca();
% hold on;
% imshow(large_thumbnail_io);
% hold(ax, 'on');
% for i = 1:size(this_cluster_boundary, 2)
%     this_cluster_boundary{i}{1} = this_cluster_boundary{i}{1}(~cellfun('isempty', this_cluster_boundary{i}{1}));
%     for j = 1:size(this_cluster_boundary{i}{1}, 2)
%         plot(tumour_polygon{i}{j});
%         hold on;
%         plot(tumour_buffer{i}{j})
%     end
%     hold on;
% end
% hold(ax, 'off')


%%
%intersecting regions

for i = 1:size(in_core, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        tumour_polygon_in{i}{j} = intersect(tumour_polygon{i}{j}, core_polygon{i});
        tumour_buffer_in{i}{j} = intersect(tumour_buffer{i}{j}, core_polygon{i});
    end
end


%%
%calculate area
for i = 1:size(in_core, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        tumour_polygon_area{i}{j} = area(tumour_polygon_in{i}{j});
        tumour_buffer_area{i}{j} = area(tumour_buffer_in{i}{j});
    end
end


%get the lymphocyte counts
for i = 1:size(in_core, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        lymph_in{i}{j}=inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        lymph_buffer{i}{j} = inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
    end
end


%lymphocyte counts divided by area
for i = 1:size(in_core, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        density_in{i}{j} = sum(lymph_in{i}{j})/tumour_polygon_area{i}{j};
        density_buffer{i}{j}=sum(lymph_buffer{i}{j})/tumour_buffer_area{i}{j};
    end
end

%combining density into a single vector for in, and another single vector for buffer
 density_in_combined = cell(0);
density_buffer_combined = cell(0);
for i = 1:size(in_core, 2)
        density_in_combined = [density_in_combined density_in{i}];
        density_buffer_combined=[density_buffer_combined density_buffer{i}];
end

csvwrite('619857_in.csv', density_in_combined)
csvwrite('619857_buffer.csv', density_buffer_combined)


