%% this is for plotting different tumour cluster size

%% A test list of images
% these are distinct ones with different IT/PT distribution
image_list=[619571]; %, 593987, 619857, 619872, 619905, 625951];
%image_list=[603288];

%%image_list = [597786];
% next. not sure if output is in the right format though

% need to run: ,

% done
%603283, 619799, 619857, 619872, 593960,593971, 594006, 602915, 602942, 602976, 602994, 603253, 603269, 603271,  603298, 603956, 619905

% no lymphocyte cluster at 20:
% 603956,625333,619538,603919,607166,619837,619942,625333,625887,625936,625951,626047, 626172


%[626047]; %[619877]; %, 626172, 625946]; %this is a list that has had issues on hpc

%done: 603271,603253, 602994,604094,597786,

%% load images and fits files

%function tum_lymph_cluster_hpc(image_filenumber)
for image = 1:size(image_list,2)
    
    image_filenumber = image_list(image);
    
    % set parameters
    % % % % % % % %     cluster_size = [5];
    % % % % % % % %     lymph_cluster_size = [20];
    % % % % % % % %     lymph_cutoff_size = [50];
    % % % % % % % %
    % % % % % % % %     % make a folder to save files
    % % % % % % % %     mkdir(num2str(image_filenumber));
    % % % % % % % %
    % % % % % % % %     % loading files
    % % % % % % % %     data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
    % % % % % % % %     info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
    % % % % % % % %     image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];
    % % % % % % % %
    % % % % % % % %     %     data= load(['/rds-d4/user/ww234/hpc-work/itpt' num2str(image_filenumber) '.mat']);
    % % % % % % % %     %     image_path = ["/rds-d4/user/ww234/hpc-work/itpt" num2str(image_filenumber) '.svs'];
    % % % % % % % %
    % % % % % % % %     %create indexing
    % % % % % % % %     X_ind = 3;
    % % % % % % % %     Y_ind = 4;
    % % % % % % % %     overlap = 22;
    % % % % % % % %     s2n=61;
    % % % % % % % %     cell_ind = 62;
    % % % % % % % %
    % % % % % % % %     %% tidying up data
    % % % % % % % %
    % % % % % % % %     %trim data
    % % % % % % % %     data_trimmed = data;
    % % % % % % % %
    % % % % % % % %     %remove non cells
    % % % % % % % %     for i = 1:size(data, 2)
    % % % % % % % %         data_trimmed{i} = data_trimmed{i}(data{cell_ind}~=0);
    % % % % % % % %     end
    % % % % % % % %     %remove overlaps
    % % % % % % % %     data_tmp = data_trimmed;
    % % % % % % % %     for i = 1:size(data_trimmed, 2)
    % % % % % % % %         data_trimmed{i} = data_trimmed{i}(data_tmp{overlap} =='F');
    % % % % % % % %     end
    % % % % % % % %     %remove low signal to noise ratio (set threshold at 1.3 - discussed with A Dariush)
    % % % % % % % %     data_tmp = data_trimmed;
    % % % % % % % %     for i = 1:size(data_trimmed, 2)
    % % % % % % % %         data_trimmed{i} = data_trimmed{i}(data_tmp{s2n} >= 1.3);
    % % % % % % % %     end
    % % % % % % % %
    % % % % % % % %     % and combine tumour and normal
    % % % % % % % %     data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;
    % % % % % % % %
    % % % % % % % %
        %% Create a boundary around each core
        %xxxx run this next
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
    
            figure % Plot the results for sanity
            imshow(large_thumbnail_io)
            hold on;
    
            for i = 1:size(grain,2)
                plot(this_boundary{i}(:,2),this_boundary{i}(:,1),'g','LineWidth',1);
            end
% %     
    % % % % % % % %
    % % % % % % % %
    % % % % % % % %     % Now convert each into a polygon
    % % % % % % % %     core_polygon= cell(0);
    % % % % % % % %     for i = 1:size(this_boundary, 2)
    % % % % % % % %         core_polygon{i} = polyshape (this_boundary{i}(:,2),this_boundary{i}(:,1)); %note the doordinate is this way round
    % % % % % % % %         %plot(core_polygon{i}); %in case the coordinates flipped
    % % % % % % % %     end
    % % % % % % % %
    % % % % % % % %
    % % % % % % % %     % the useful output of this chunck is the core_polygon cell array
    % % % % % % % %
    % % % % % % % %     %% now subset coordinates of cells
    % % % % % % % %     in_polygon=cell(0);
    % % % % % % % %     in_core=cell(0);
    % % % % % % % %
    % % % % % % % %     for i = 1:size(core_polygon, 2)
    % % % % % % % %         in_polygon{i} = inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, core_polygon{i}.Vertices(:,1), core_polygon{i}.Vertices(:,2));
    % % % % % % % %         for j = 1:size(data_trimmed, 2)
    % % % % % % % %             in_core{i}{j} = data_trimmed{j}(in_polygon{i}, :);
    % % % % % % % %         end
    % % % % % % % %     end
    % % % % % % % %
    
    
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/workspace1/'  num2str(image_filenumber) '_workspace1.mat']);
    cluster_size = [5];
    lymph_cluster_size = [5];
    lymph_cutoff_size = [50];
    data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
    info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
    image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];
    
    thumbnail_io=imread(image_path,'Index',thumbnail_layer);
    large_thumbnail_io = imresize(thumbnail_io,thumbnail_overall_scale_factor);
% %     
% %     figure % Plot the results for sanity
% %     imshow(large_thumbnail_io)
% %     hold on;
% %     for i = 1:size(core_polygon, 2)
% %         plot(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), 'b.')
% %           %in case the cells are in the wrong place
% %     end
    
    %save(['./' num2str(image_filenumber) '/workspace1.mat']);
    
    
    %% Now compute the tumour and lymphocyte clusters with subsetting
    
    %24th july 2019: changed epsilon to epsilon/2, now works well.  epsilon*3/4 is not very
    %tight, and epsilon/3 is splitting up lymphoid aggregates into very silly small bits.
    %settle with epsilon/2
    
    core_list = [];
    
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
        
% % % % % % % %         base_cells = in_core{this_core}{cell_ind}==2;
% % % % % % % %         neighbour_cells = in_core{this_core}{cell_ind}==2;
% % % % % % % %         
% % % % % % % %         %ignore the really small chunks that has very few tumour cells
% % % % % % % %         if size(in_core{this_core}{cell_ind}==2,1) < 10
% % % % % % % %             continue
% % % % % % % %         end
% % % % % % % %         
% % % % % % % %         [all_multi_real_distances{this_core}, all_indexes{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', lymph_cluster_size+1);
% % % % % % % %         
% % % % % % % %         i = 0;
% % % % % % % %         if size(all_multi_real_distances{this_core},1)<max(lymph_cluster_size)+1
% % % % % % % %             all_multi_real_distances{this_core}(size(all_multi_real_distances{this_core},1)+1:max(lymph_cluster_size)+1,:)=NaN;
% % % % % % % %             continue;
% % % % % % % %         end
% % % % % % % %         
% % % % % % % %         %now use DBScan
% % % % % % % %         for this_clustsize = lymph_cluster_size
% % % % % % % %             i = i+1;
% % % % % % % %             C=0;
% % % % % % % %             
% % % % % % % %             n=sum(neighbour_cells);
% % % % % % % %             cluster_membership{i}=zeros(n,1);
% % % % % % % %             
% % % % % % % %             visited=false(n,1);
% % % % % % % %             isnoise=false(n,1);
% % % % % % % %             
% % % % % % % %             %Calculate epsilon from the knee point of a
% % % % % % % %             %1000-bin histogram to the k-1th nearest neighbour
% % % % % % % %             %Remember that first index is self, so no need
% % % % % % % %             %to subtract one.
% % % % % % % %             %[hist_y, hist_x] = hist(all_multi_real_distances{this_core}(this_clustsize,:),1000);
% % % % % % % %             [hist_y, hist_x] = hist(all_multi_real_distances{this_core}(3,:),1000);
% % % % % % % %             [epsilon, ~] = knee_pt(hist_y,hist_x,true);
% % % % % % % %             
% % % % % % % %             for this_cell=1:n
% % % % % % % %                 if ~visited(this_cell)
% % % % % % % %                     visited(this_cell)=true;
% % % % % % % %                     Neighbours=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_cell)<=epsilon/2,this_cell),this_cell)';
% % % % % % % %                     if numel(Neighbours)<this_clustsize
% % % % % % % %                         % X(i,:) is NOISE
% % % % % % % %                         isnoise(this_cell)=true;
% % % % % % % %                     else
% % % % % % % %                         C=C+1;
% % % % % % % %                         cluster_membership{i}(this_cell)=C;
% % % % % % % %                         k = 1;
% % % % % % % %                         while true
% % % % % % % %                             this_neighbour_cell = Neighbours(k);
% % % % % % % %                             
% % % % % % % %                             if ~visited(this_neighbour_cell)
% % % % % % % %                                 visited(this_neighbour_cell)=true;
% % % % % % % %                                 Neighbours2=setdiff(all_indexes{this_core}(all_multi_real_distances{this_core}(:,this_neighbour_cell)<=epsilon/2,this_neighbour_cell),this_neighbour_cell)';
% % % % % % % %                                 if numel(Neighbours2)>=this_clustsize
% % % % % % % %                                     Neighbours=[Neighbours Neighbours2];   %#ok
% % % % % % % %                                 end
% % % % % % % %                             end
% % % % % % % %                             if cluster_membership{i}(this_neighbour_cell)==0
% % % % % % % %                                 cluster_membership{i}(this_neighbour_cell)=C;
% % % % % % % %                             end
% % % % % % % %                             
% % % % % % % %                             k = k + 1;
% % % % % % % %                             if k > numel(Neighbours)
% % % % % % % %                                 break;
% % % % % % % %                             end
% % % % % % % %                         end
% % % % % % % %                     end
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %             
% % % % % % % %             x_data_here = in_core{this_core}{X_ind}(neighbour_cells);
% % % % % % % %             y_data_here = in_core{this_core}{Y_ind}(neighbour_cells);
% % % % % % % %             
% % % % % % % %             for this_cluster_number = 1:max(cluster_membership{i})
% % % % % % % %                 this_x_subset = x_data_here(cluster_membership{i}==this_cluster_number);
% % % % % % % %                 this_y_subset = y_data_here(cluster_membership{i}==this_cluster_number);
% % % % % % % %                 this_lymphocyte_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
% % % % % % % %                 this_lymphocyte_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_lymphocyte_cluster_boundary_numbers),this_y_subset(this_lymphocyte_cluster_boundary_numbers)];
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %         
        core_list = [core_list this_core];
        
    end
    
    %save(['./' num2str(image_filenumber) '/workspace_clusters.mat']);
    
    %% this bit is to match the dimension of the cluster boundaries with the core_list, as we
    % may need empty cells at the end of the list for cores without either tumour clusters or
    % % % % % %     % lymphocyte clusters
    % % % % % %
    % % % % % %     % add a line to say if there is no lymphocyte cluster then everything is zero
    % % % % % %     if exist('this_lymphocyte_cluster_boundary', 'var') == 0
    % % % % % %         save(['./' num2str(image_filenumber) '/workspace_clusters20.mat']);
    % % % % % %         lymphocyte_polygon = [];
    % % % % % %         t_l_intersection = [];
    % % % % % %         tbuffer_l_intersection = [];
    % % % % % %
    % % % % % %         t_l_intersection_lymph_count=0;
    % % % % % %         t_l_intersection_area =0;
    % % % % % %         tbuffer_l_intersection_lymph_count=0;
    % % % % % %         tbuffer_l_intersection_area=0;
    % % % % % %         lymphocyte_polygon_area_combined=0;
    % % % % % %         lymph_lymphcluster_count_combined=0;
    % % % % % %
    % % % % % %         save(['./' num2str(image_filenumber) '/t_l_intersection20.mat'], 't_l_intersection');
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_intersection_count20.csv'], t_l_intersection_lymph_count);
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersection_area20.csv'], t_l_intersection_area);
    % % % % % %         save(['./' num2str(image_filenumber) '/lymphocyte_polygon20.mat'], 'lymphocyte_polygon');
    % % % % % %         save(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_tbuffer_l_intersection20.mat'], 'tbuffer_l_intersection');
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphbuffer_intersection_count20.csv'], tbuffer_l_intersection_lymph_count);
    % % % % % %         csvwrite([ '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/count/' num2str(image_filenumber) '_lymphbuffer_intersection_count20.csv'], tbuffer_l_intersection_lymph_count);
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_intersectionbuffer_area20.csv'], tbuffer_l_intersection_area);
    % % % % % %         csvwrite(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/area/' num2str(image_filenumber) '_intersectionbuffer_area20.csv'], tbuffer_l_intersection_area);
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_area20.csv'], lymphocyte_polygon_area_combined);
    % % % % % %         csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymph_lymphcluster_count20.csv'], lymph_lymphcluster_count_combined);
    % % % % % %         continue
    % % % % % %     end
    % % % % % %
    % % % % % %
    % % % % % %     % firstly work out the number of cores in total (max of core_total), and how many are in
    % % % % % %     % each list. then work out the difference and append.
    % % % % % %     tumour_core_total = size(this_tumour_cluster_boundary, 2);
    % % % % % %     lymph_core_total = size(this_lymphocyte_cluster_boundary, 2);
    % % % % % %     core_total = max(core_list);
    % % % % % %
    % % % % % %     % if same sizes then dont' worry about it....
    % % % % % %     if core_total == tumour_core_total
    % % % % % %         %do nothing
    % % % % % %     else
    % % % % % %         tumour_diff = core_total - tumour_core_total;
    % % % % % %         for diff = 1:tumour_diff
    % % % % % %             this_tumour_cluster_boundary{1, tumour_core_total + diff} = [];
    % % % % % %         end
    % % % % % %     end
    % % % % % %
    % % % % % %     if core_total == lymph_core_total
    % % % % % %         %do nothing
    % % % % % %     else
    % % % % % %         lymph_diff = core_total - lymph_core_total;
    % % % % % %         for diff = 1:lymph_diff
    % % % % % %             this_lymphocyte_cluster_boundary{1, lymph_core_total + diff} = [];
    % % % % % %         end
    % % % % % %     end
    % % % % % %
    % % % % % %    save(['./' num2str(image_filenumber) '/workspace_clusters20.mat']);
    
    % this is to check if the clusters are correct
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
        if isempty(this_tumour_cluster_boundary{i}{1}) == 1
            continue
        end
        
        % this_tumour_cluster_boundary{i}{1} =
        % this_tumour_cluster_boundary{i}{1}(~cellfun('isempty',
        % this_tumour_cluster_boundary{i}{1})); % this may be introducing extra layers
        % unnecessarily
        for j = 1:size(this_tumour_cluster_boundary{i}{1}, 2)
            hold on;
            if isempty(this_tumour_cluster_boundary{i}{1}{j}) == 1
                continue
            end
            plot(this_tumour_cluster_boundary{i}{1}{j}(:,1), this_tumour_cluster_boundary{i}{1}{j}(:,2), 'k-', 'LineWidth',1);
        end
        hold on;
    end
    %
    % % now draw the lymphocyte cluster
% % % % %     for i = 1:size(this_lymphocyte_cluster_boundary, 2)
% % % % %         if isempty(this_lymphocyte_cluster_boundary{i}) == 1
% % % % %             continue
% % % % %         end
% % % % %         if isempty(this_lymphocyte_cluster_boundary{i}{1}) == 1
% % % % %             continue
% % % % %         end
% % % % %         %this_lymphocyte_cluster_boundary{i}{1} =
% % % % %         %this_lymphocyte_cluster_boundary{i}{1}(~cellfun('isempty',
% % % % %         %this_lymphocyte_cluster_boundary{i}{1})); % this may be introducing extralayers
% % % % %         %unnecesarily
% % % % %         for j = 1:size(this_lymphocyte_cluster_boundary{i}{1}, 2)
% % % % %             hold on;
% % % % %             if isempty(this_lymphocyte_cluster_boundary{i}{1}{j}) == 1
% % % % %                 continue
% % % % %             end
% % % % %             plot(this_lymphocyte_cluster_boundary{i}{1}{j}(:,1), this_lymphocyte_cluster_boundary{i}{1}{j}(:,2), 'r-');
% % % % %         end
% % % % %         hold on;
    
    hold(ax, 'off')
end