% This script requires an image filek (.svs) and its matching cell
% segmentation file (.fits).
%
%this is a total list
%image_list = [603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160];
%,603288, 603298,603292


%9th Dec 2019 it has been noticed that the older way of indexing are excluding tumour
%clusters without any lymphocytes overlap and hence not a good idea. MOdification today is
%to re-index overlap areas

%% A  list of images
% these are distinct ones with different IT/PT distribution
% image_list=[603288, 593987, 619857, 619872, 619905, 625951];


% % % % % %% this is for debugging purpose
% % % % % image_list= [593708];
% % % % %
% % % % %
% % % % % %[593971, 594006, 602915, 602942, 602976, 602994, 603253, 603269, 603271, 603283, 603298];
% % % % %
% % % % % % repeat: 593960,
% % % % % % [603922, 626172,597786,594110,  ];

% % % % % for image = 1:size(image_list,2)
% % % % %     image_filenumber = image_list(image);
% % % % %
% % % % %     image_path_stem = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis';
% % % % %
% % % % %     data = load([image_path_stem '/mat_file_new/' num2str(image_filenumber) '.mat']);
% % % % %     image_path = [image_path_stem '/IT_PT_zone/' num2str(image_filenumber) '.svs'];
% % % % %
% % % % %

%% load images and fits files

function tum_lymph_cluster_hpc(image_filenumber_fullpath)
disp(image_filenumber_fullpath)
disp(class(image_filenumber_fullpath))
if ischar(image_filenumber_fullpath)
    [image_path_stem,image_filenumber,~] = fileparts(image_filenumber_fullpath)
    image_filenumber = str2num(image_filenumber);
    if isempty(image_path_stem)
        image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt'
        warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
    end
elseif isnumeric(image_filenumber_fullpath)
    image_filenumber = image_filenumber_fullpath;
    image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt'
    warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
end

%%%loading files
% % % data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
% % % info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
% % % image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];
% %
    data= load([image_path_stem '/' num2str(image_filenumber) '.mat']);
    image_path = [image_path_stem '/' num2str(image_filenumber) '.svs'];

%%%%%%%%%%%%%%%% comment from function to here for debugging

% make a folder to save files
mkdir(num2str(image_filenumber))

%% set parameters and index

% set parameters
cluster_size = [5];
lymph_cluster_size = [5];
buffer_size = 100;

% and index
X_ind = 1;
Y_ind = 2;
cell_ind = 3;
% overlap = 22;
% s2n=61;

%% tidying up data

%trim data
data_trimmed = data.data_mini;
data_trimmed = num2cell(data_trimmed, 1); %so that the cell index works later on

%removed non cells, overlaps, low signal to noise ratio (set threshold at 1.3 - discussed with A Dariush) and combine tumour and normal

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

%     figure % Plot the results for sanity
%     imshow(large_thumbnail_io)
%     hold on;
%
%     %for i = 1:size(grain,2)
%         plot(this_boundary{i}(:,2),this_boundary{i}(:,1),'g','LineWidth',1);
%     %end

% Now convert each into a polygon
core_polygon= cell(0);
for i = 1:size(this_boundary, 2)
    core_polygon{i} = polyshape (this_boundary{i}(:,2),this_boundary{i}(:,1)); %note the doordinate is this way round
    %plot(core_polygon{i}); %in case the coordinates flipped
end


% the useful output of this chunck is the core_polygon cell array

%% now subset coordinates of cells
in_polygon=cell(0);
in_core=cell(0);

for i = 1:size(core_polygon, 2)
    in_polygon{i} = inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, core_polygon{i}.Vertices(:,1), core_polygon{i}.Vertices(:,2));
    for j = 1:size(data_trimmed, 2)
        in_core{i}{j} = data_trimmed{j}(in_polygon{i}, :);
    end
end

save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_workspace1.mat']);

%% clear and reload workspace for subsequent work

% for final code, change reload workspace to =0 to skip this step (true or 1, false or 0)
reload_workspace = 0;
if reload_workspace
    save('image_filenumber','image_filenumber');
    clear variables;
    load('image_filenumber');
    load(['./' num2str(image_filenumber) '/workspace.mat']);
end


%% Now compute the tumour and lymphocyte clusters with subsetting

%24th july 2019: changed epsilon to epsilon/2, now works well.  epsilon*3/4 is not very
%tight, and epsilon/3 is splitting up lymphoid aggregates into very silly small bits.
%settle with epsilon/2

%core_list = [];
core_list=[1:1:size(core_polygon,2)];

for this_core = 1:size(core_polygon,2)
    
    %tumour clusters
    
    base_cells = in_core{this_core}{cell_ind}==1;
    neighbour_cells = in_core{this_core}{cell_ind}==1;
    
    %ignore the really small chunks that has very few tumour cells
    %             if size(in_core{this_core}{cell_ind}==1,1) < 10
    %                 continue
    %             end
    
    [all_multi_real_distances{this_core}, all_indexes{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
    
    i = 0;
    if size(all_multi_real_distances{this_core},1)<max(cluster_size)+1
        all_multi_real_distances{this_core}(size(all_multi_real_distances{this_core},1)+1:max(cluster_size)+1,:)=NaN;
        this_tumour_cluster_boundary{this_core} = [];
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
        [hist_y, hist_x] = hist(all_multi_real_distances{this_core}(this_clustsize,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        
        for this_cell=1:n
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
        
        x_data_here = in_core{this_core}{X_ind}(neighbour_cells);
        y_data_here = in_core{this_core}{Y_ind}(neighbour_cells);
        
        for this_cluster_number = 1:max(cluster_membership{i})
            this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
            this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
            this_tumour_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
            this_tumour_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_tumour_cluster_boundary_numbers),this_y_subset(this_tumour_cluster_boundary_numbers)];
        end
    end
    % the output that feed into the next step is this_tumour_cluster_boundary
    
    
    
    % lymphocyte clusters
    
    base_cells = in_core{this_core}{cell_ind}==2;
    neighbour_cells = in_core{this_core}{cell_ind}==2;
    
    %ignore the really small chunks that has very few lymph cells
    %     if size(in_core{this_core}{cell_ind}==2,1) < 10
    %         continue
    %     end
    
    [all_multi_real_distances{this_core}, all_indexes{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', lymph_cluster_size+1);
    
    i = 0;
    if size(all_multi_real_distances{this_core},1)<max(lymph_cluster_size)+1
        all_multi_real_distances{this_core}(size(all_multi_real_distances{this_core},1)+1:max(lymph_cluster_size)+1,:)=NaN;
        this_lymphocyte_cluster_boundary{this_core} = [];
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
        [hist_y, hist_x] = hist(all_multi_real_distances{this_core}(3,:),1000);
        [epsilon, ~] = knee_pt(hist_y,hist_x,true);
        
        for this_cell=1:n
            if ~visited(this_cell)
                visited(this_cell)=true;
                Neighbours=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_cell)<=epsilon/2,this_cell),this_cell)';
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
                            Neighbours2=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_neighbour_cell)<=epsilon/2,this_neighbour_cell),this_neighbour_cell)';
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
        
        x_data_here = in_core{this_core}{X_ind}(neighbour_cells);
        y_data_here = in_core{this_core}{Y_ind}(neighbour_cells);
        
        for this_cluster_number = 1:max(cluster_membership{i})
            this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
            this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
            this_lymphocyte_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
            this_lymphocyte_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_lymphocyte_cluster_boundary_numbers),this_y_subset(this_lymphocyte_cluster_boundary_numbers)];
        end
    end
    
    %core_list = [core_list this_core];
    
end


%% this bit is to match the dimension of the cluster boundaries with the core_list, as we
% may need empty cells at the end of the list for cores without either tumour clusters or
% lymphocyte clusters

% firstly work out the number of cores in total (max of core_total), and how many are in
% each list. then work out the difference and append.
tumour_core_total = size(this_tumour_cluster_boundary, 2);
lymph_core_total = size(this_lymphocyte_cluster_boundary, 2);
core_total = max(core_list);

% if same sizes then dont' worry about it....
if core_total == tumour_core_total
    %do nothing
else
    tumour_diff = core_total - tumour_core_total;
    for diff = 1:tumour_diff
        this_tumour_cluster_boundary{1, tumour_core_total + diff} = [];
    end
end

if core_total == lymph_core_total
    %do nothing
else
    lymph_diff = core_total - lymph_core_total;
    for diff = 1:lymph_diff
        this_lymphocyte_cluster_boundary{1, lymph_core_total + diff} = [];
    end
end


% save(['./' num2str(image_filenumber) '/'  num2str(image_filenumber) '_workspace_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat']);


%% now convert each into a polygon - old code
% %the this_tumour/lymphocyte_cluster_boundary is organised in such a way that it has core number
% %first, then some times it has extra layers before reaches the collection of number for
% %the clusters. so will unnest within each core
for i =1:size(core_list, 2) %First un-nest the cells
    if isempty(this_tumour_cluster_boundary{core_list(i)}) == 1
        this_tumour_cluster_boundary{core_list(i)} = [];
        continue
    end
    this_tumour_cluster_boundary{core_list(i)} = extractmycells(this_tumour_cluster_boundary{core_list(i)});
end

for i =1:size(core_list, 2)
    if isempty(this_lymphocyte_cluster_boundary{core_list(i)}) ==1
        %this_lymphocyte_cluster_boundary{core_list(i)} = [];
        continue
    end
    this_lymphocyte_cluster_boundary{core_list(i)} = extractmycells(this_lymphocyte_cluster_boundary{core_list(i)});
end


for i =1:size(core_list, 2)
    if isempty(this_tumour_cluster_boundary{core_list(i)}) == 1
        tumour_polygon{core_list(i)} = [];
        continue
    end
    
    for j = 1:size(this_tumour_cluster_boundary{core_list(i)}, 2)
        if isempty(this_tumour_cluster_boundary{core_list(i)}{j}) == 1
            tumour_polygon{core_list(i)}{j} = [];
            continue
        end
        tumour_polygon{core_list(i)}{j} = polyshape (this_tumour_cluster_boundary{core_list(i)}{j});
        %plot(tumour_polygon{core_list(i)}{j});
    end
end

for i =1:size(core_list, 2)
    if isempty(this_lymphocyte_cluster_boundary{core_list(i)}) == 1
        lymphocyte_polygon{core_list(i)} = [];
        continue
    end
    for k = 1:size(this_lymphocyte_cluster_boundary{core_list(i)}, 2)
        if isempty(this_lymphocyte_cluster_boundary{core_list(i)}{k}) == 1
            lymphocyte_polygon{core_list(i)}{k} = [];
            continue
        end
        lymphocyte_polygon{core_list(i)}{k} = polyshape(this_lymphocyte_cluster_boundary{core_list(i)}{k});
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
%     % % %present etc.,
%     total_core_polygon = 1:size(this_tumour_cluster_boundary, 2);
%
%     tumour_core_list = intersect(core_list, total_core_polygon(~cellfun('isempty',this_tumour_cluster_boundary))); %In case you need this later
%     % %this_tumour_cluster_boundary = this_tumour_cluster_boundary(~cellfun('isempty',this_tumour_cluster_boundary));
%     lymphocyte_core_list = intersect(core_list, total_core_polygon(~cellfun('isempty',this_lymphocyte_cluster_boundary))); %In case you need this later
%

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

%%
%this bit then works out the buffer bit

tumour_buffer = cell(0);
for i = 1:size(tumour_polygon, 2)
    for j = 1:size(tumour_polygon{i}, 2)
        if isempty(tumour_polygon{i}{j})
            continue
        end
        tumour_buffer{i}{j} = polybuffer(tumour_polygon{i}{j},buffer_size);
        %tumour_buffer{core_list(i)}{j} = polyshape(tumour_buffer{core_list(i)}{j});
        %this is not necessary as it is already a polyshape
    end
end


%% exclude the area that falls outside the cores
%intersecting regions so only within core area is calculated
%
for i = 1:size(core_list, 2)
    %     for j = 1:size(tumour_polygon{i}, 2)
    %         if isempty(tumour_polygon{i}{j})
    %             continue
    %         end
    %         tumour_polygon_in{i}{j} = intersect(tumour_polygon{i}{j}, core_polygon{tumour_core_list(i)});
    %     end
    % the above is unnecessary as only the buffer goes outside
    if isempty(tumour_buffer{core_list(i)})
        continue
    end
    for k = 1:size(tumour_buffer{core_list(i)},2)
        if isempty(tumour_buffer{core_list(i)}{k})
            continue
        end
        tumour_buffer_in{core_list(i)}{k} = intersect(tumour_buffer{core_list(i)}{k}, core_polygon{core_list(i)});
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
save(['./' num2str(image_filenumber) '/'  num2str(image_filenumber) '_workspace_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat']);

save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_polygon_in_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat'], 'tumour_polygon_in');
save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tumour_buffer_in_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat'], 'tumour_buffer_in');
save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_polygon_l' num2str(lymph_cluster_size) '_b' num2str(buffer_size) '.mat'], 'lymphocyte_polygon');

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

for i=1:size(tumour_polygon_in, 2) %selecting cores with tumour polygons
    if isempty(tumour_polygon_in{i}) %in case some cores dont' have tumour
        continue
    end
    %         t_l_intersection{i} = [];
    %         t_l_intersection_area{i}={};
    %         t_l_intersection_lymph_count{i}= {};
    %         t_l_centroid_count{i} = {};
    
    if isempty(lymphocyte_polygon{i}) %if there is no lymphocyte in this core then all value for this becomes empty
        t_l_intersection{i} = cell(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        t_l_intersection_area{i} = zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        t_l_intersection_lymph_count{i}= zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        t_l_centroid_count{i} = zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        tbuffer_l_intersection{i} = cell(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        tbuffer_l_intersection_area{i} = zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        tbuffer_l_intersection_lymph_count{i}= zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        tbuffer_l_centroid_count{i} = zeros(size(tumour_polygon_in{i}, 1), size(tumour_polygon_in{i}, 2));
        %this has made the corresponding cell for polygons, area and counts as
        %empty or zero respectively
        continue
    end
    
    for j = 1:size(tumour_polygon_in{i}, 2) %indexing with tumour polygons
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
        
        if isempty(tumour_polygon_in{i}{j})
            continue
        end
        
        
        for k = 1:size(lymphocyte_polygon{i}, 2)
            %now looking at the overlap between the particular lymphocyte polygon and
            %tumour polygon
            tmp = intersect(tumour_polygon_in{i}{j}, lymphocyte_polygon{i}{k});
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
        
        lymph_centroid_in_logi{i}{j} = inpolygon(lymph_centroid_in{i}{j}(:,1), lymph_centroid_in{i}{j}(:,2), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        t_l_centroid_count{i}{j} = sum(lymph_centroid_in_logi{i}{j});
        
        
        for l = 1:size(lymphocyte_polygon{i}, 2)
            %now looking at the overlap between the particular lymphocyte polygon and
            %tumour buffer
            tmpb = intersect(tumour_buffer_in{i}{j}, lymphocyte_polygon{i}{l});
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
        lymph_centroid_buffer_logi{i}{j} = inpolygon(lymph_centroid_buffer{i}{j}(:,1), lymph_centroid_buffer{i}{j}(:,2), tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
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


for i = 1:size(tumour_polygon_in, 2)
    %tumour polygon area
    if isempty(tumour_polygon_in{i})
        tumour_polygon_area{i} = [];
        tumour_buffer_area{i} = [];
        continue
    end
    for j = 1:size(tumour_polygon_in{i}, 2)
        if isempty(tumour_polygon_in{i}{j})
            tumour_polygon_area{i}{j} = [];
            tumour_buffer_area{i}{j} = [];
            continue
        end
        tumour_polygon_area{i}{j} = area(tumour_polygon_in{i}{j});
        tumour_buffer_area{i}{j} = area(tumour_buffer_in{i}{j});
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
for i = 1:size(tumour_polygon_in, 2)
    if isempty(tumour_polygon_in{i})
        tumourcell_in_count{i} = [];
        tumourcell_buffer_count{i}=[];
        continue
    end
    
    for j = 1:size(tumour_polygon_in{i}, 2)
        if isempty(tumour_polygon_in{i}{j})
            tumourcell_in_count{i}{j} = [];
            tumourcell_buffer_count{i}{j}=[];
            lymph_in_count{i}{j} = [];
            continue
        end
        %now we calculate the number of lymphocyte within tumour polygon
        lymph_in{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        lymph_in_count{i}{j} = sum(lymph_in{i}{j});
        %now we calculate the number of tumour cells within tumour polygon
        tumourcell_in{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), tumour_polygon_in{i}{j}.Vertices(:,1), tumour_polygon_in{i}{j}.Vertices(:,2));
        tumourcell_in_count{i}{j} = sum(tumourcell_in{i}{j});
        % number of lymphocyte in buffer
        lymph_buffer{i}{j} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
        lymph_buffer_count{i}{j} = sum(lymph_buffer{i}{j});
        %number of tumour cells in buffer
        tumourcell_buffer{i}{j}=inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), tumour_buffer_in{i}{j}.Vertices(:,1), tumour_buffer_in{i}{j}.Vertices(:,2));
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
            lymph_lymphcluster_count{i}{j} = [];
            tumour_lymphcluster_count{i}{j} = [];
            continue
        end
        lymph_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==2), in_core{i}{Y_ind}(in_core{i}{cell_ind}==2), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
        lymph_lymphcluster_count{i}{l} = sum(lymph_lymphcluster{i}{l});
        tumour_lymphcluster{i}{l} = inpolygon(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{lymphocyte_core_list(i)}{Y_ind}(in_core{i}{cell_ind}==1), lymphocyte_polygon{i}{l}.Vertices(:,1), lymphocyte_polygon{i}{l}.Vertices(:,2));
        tumour_lymphcluster_count{i}{l} = sum(lymph_lymphcluster{i}{l});
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

end
