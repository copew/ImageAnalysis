% This script requires an image filek (.svs) and its matching cell
% segmentation file (.fits).
%
%this is a total list
%image_list = [603298,603292,603288,603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160];
%
% buffer zone width 0.5mm

    %image_list = [603298,603292,603288,603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971, 602962];

% image_list=[603288] %, 593987, 619857] %, 619872, 619905, 625951];

%%
function cluster_buffer_mat(image_filenumber)
%for image = 1:size(image_list,2)
    %,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160];


    %image_filenumber = image_list(image);
    cluster_size = [5];
   
    % loading files
    %data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
    %info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
    data=load(['/rds-d4/user/ww234/hpc-work/itpt/' num2str(image_filenumber) '.mat']);
    image_path = ['/rds-d4/user/ww234/hpc-work/itpt/' num2str(image_filenumber) '.svs'];
    
    %create indexing
    X_ind = 1;
    Y_ind = 2;
    cell_ind = 3;
    
    %trim data
    
    for i = 1:size(data.data_mini, 2)
	data_trimmed{i} = data.data_mini(:,i);
        data_trimmed{i} = data_trimmed{i}(data.data_mini(:,cell_ind)~=0);
    end

    % Now remove duplicate points from cell identification tiling
    full_coords = zeros(1,length(data_trimmed{Y_ind}));
    for i = 1:length(data_trimmed{Y_ind}) 
        full_coords(i) = [data_trimmed{X_ind}(i)*1000000000 + data_trimmed{Y_ind}(i)]; 
    end
    [~, unique_coords, ~] = unique(full_coords,'stable');
    for i = 1:length(data_trimmed)
        data_trimmed{i} = data_trimmed{i}(unique_coords);
    end
    
    % and combine tumour and normal
    data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;
    
    
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
%     for i = 1:size(grain,2)
%         plot(this_boundary{i}(:,2),this_boundary{i}(:,1),'g','LineWidth',1);
%     end
    
    % Now convert each into a polygon
    core_polygon= cell(0);
    for i = 1:size(this_boundary, 2)
        core_polygon{i} = polyshape (this_boundary{i}(:,2),this_boundary{i}(:,1)); %note the doordinate is this way round
        %plot(core_polygon{i}); %in case the coordinates flipped
    end
    
    
    % the useful output of this chunck is the core_polygon cell array
    
    %%
    % now subset coordinates of cells
    in_polygon=cell(0);
    in_core=cell(0);

    for i = 1:size(core_polygon, 2)
        in_polygon{i} = inpolygon(data_trimmed{X_ind}, data_trimmed{Y_ind}, core_polygon{i}.Vertices(:,1), core_polygon{i}.Vertices(:,2));
        for j = 1:size(data_trimmed, 2)
            in_core{i}{j} = data_trimmed{j}(in_polygon{i}, :);
        end
    end

    % figure % Plot the results for sanity
    % imshow(large_thumbnail_io)
    % hold on;
    % for i = 1:size(core_polygon, 2)
    %     plot(in_core{i}{X_ind}(in_core{i}{cell_ind}==1), in_core{i}{Y_ind}(in_core{i}{cell_ind}==1), 'b.') 
    %       % in case the cells are in the wrong place
    % end
    
    
    
    
%% 
 % Now compute the tumour clusters with subsetting
   
    core_list = [];
    
    for this_core = 1:size(core_polygon,2)
        base_cells = in_core{this_core}{cell_ind}==1;
        neighbour_cells = in_core{this_core}{cell_ind}==1;
        
        %ignore the really small chunks that has very few tumour cells
        if size(in_core{this_core}{cell_ind}==1,1) < 10
            continue
        end
    
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
                        C=C+1
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
                this_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset);
                this_cluster_boundary{this_core}{i}{this_cluster_number} = [this_x_subset(this_cluster_boundary_numbers),this_y_subset(this_cluster_boundary_numbers)];
            end
        end
    core_list = [core_list this_core];    
    end
    
 
    %%
    %this is to check if the clusters are correct - and it is... phew
%     figure
%     ax=gca();
%     hold on;
%     imshow(large_thumbnail_io);
%     hold(ax, 'on');
%     for i = 1:size(this_cluster_boundary, 2)
%         if isempty(this_cluster_boundary{i}) == 1;
%             continue
%         end
%         this_cluster_boundary{i}{1} = this_cluster_boundary{i}{1}(~cellfun('isempty', this_cluster_boundary{i}{1}));
%         for j = 1:size(this_cluster_boundary{i}{1}, 2)
%             hold on;
%             plot(this_cluster_boundary{i}{1}{j}(:,1), this_cluster_boundary{i}{1}{j}(:,2), 'k-');
%         end
%         hold on
%     end
%     hold(ax, 'off')
%   
  %now convert each into a polygon
empty_cores = [];
for i =1:size(core_list, 2)
        if isempty(this_cluster_boundary{i}) == 1;
 empty_cores(end+1) = i; %Make a list of empty cores        
    continue
        end
        this_cluster_boundary{core_list(i)}{1} = this_cluster_boundary{core_list(i)}{1}(~cellfun('isempty', this_cluster_boundary{core_list(i)}{1}));
        for j = 1:size(this_cluster_boundary{core_list(i)}{1}, 2)
            tumour_polygon{core_list(i)}{j} = polyshape (this_cluster_boundary{core_list(i)}{1}{j});
%             plot(tumour_polygon{core_list(i)}{j});
        end
    end
core_list(empty_cores) = []; %Remove the empty core from the processing list    

    %the end result is tumour_polygon
    
    %%
    %this bit then works out the buffer bit
    
    tumour_buffer = cell(0);
    
    for i = 1:size(core_list, 2)
        for j = 1:size(this_cluster_boundary{core_list(i)}{1}, 2)
            tumour_buffer{core_list(i)}{j} = polybuffer(tumour_polygon{core_list(i)}{j},100);
        end
    end
    
    %sanity check :)
%     figure
%     ax=gca();
%     hold on;
%     imshow(large_thumbnail_io);
%     hold(ax, 'on');
%     for i = 1:size(core_list, 2)
%         this_cluster_boundary{core_list(i)}{1} = this_cluster_boundary{core_list(i)}{1}(~cellfun('isempty', this_cluster_boundary{core_list(i)}{1}));
%         for j = 1:size(this_cluster_boundary{core_list(i)}{1}, 2)
%             plot(tumour_polygon{core_list(i)}{j});
%             hold on;
%             plot(tumour_buffer{core_list(i)}{j})
%         end
%         hold on;
%     end
%     hold(ax, 'off')
    
    
    %%
    %intersecting regions
    
    for i = 1:size(core_list, 2)
        for j = 1:size(tumour_polygon{core_list(i)}, 2)
            tumour_polygon_in{core_list(i)}{j} = intersect(tumour_polygon{core_list(i)}{j}, core_polygon{core_list(i)});
        end
        for k = 1:size(tumour_buffer{core_list(i)},2)
            tumour_buffer_in{core_list(i)}{k} = intersect(tumour_buffer{core_list(i)}{k}, core_polygon{core_list(i)});
        end
    end
    
    
    %%
    %calculate area
    for i = 1:size(core_list, 2)
        for j = 1:size(tumour_polygon{core_list(i)}, 2)
            tumour_polygon_area{core_list(i)}{j} = area(tumour_polygon_in{core_list(i)}{j});
        end
        for k = 1:size(tumour_buffer{core_list(i)},2)
            tumour_buffer_area{core_list(i)}{k} = area(tumour_buffer_in{core_list(i)}{k});
        end
    end
    
    %%
    %get the lymphocyte counts
    for i = 1:size(core_list, 2)
        for j = 1:size(tumour_polygon{core_list(i)}, 2)
            lymph_in{core_list(i)}{j}=inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==2), data_trimmed{Y_ind}(data_trimmed{cell_ind}==2), tumour_polygon_in{core_list(i)}{j}.Vertices(:,1), tumour_polygon_in{core_list(i)}{j}.Vertices(:,2));
        end
        for k = 1:size(tumour_buffer{core_list(i)},2)
            lymph_buffer{core_list(i)}{k} = inpolygon(data_trimmed{X_ind}(data_trimmed{cell_ind}==2), data_trimmed{Y_ind}(data_trimmed{cell_ind}==2), tumour_buffer_in{core_list(i)}{k}.Vertices(:,1), tumour_buffer_in{core_list(i)}{k}.Vertices(:,2));
        end
    end
    
    %%
    %lymphocyte counts divided by area
    for i = 1:size(core_list, 2)
        for j = 1:size(tumour_polygon{core_list(i)}, 2)
            density_in{core_list(i)}{j} = sum(lymph_in{core_list(i)}{j})/tumour_polygon_area{core_list(i)}{j};
        end
        for k = 1:size(tumour_buffer{core_list(i)},2)
            density_buffer{core_list(i)}{k}=sum(lymph_buffer{core_list(i)}{k})/tumour_buffer_area{core_list(i)}{k};
        end
    end
    %%
    %combining density into a single vector for in, and another single vector for buffer
    density_in_combined = cell(0);
    density_buffer_combined = cell(0);
    for i = 1:size(core_list, 2)
        density_in_combined = [density_in_combined density_in{core_list(i)}];
        density_buffer_combined=[density_buffer_combined density_buffer{core_list(i)}];
    end
    
%    csvwrite([num2str(image_filenumber) '_in.csv'], density_in_combined);
%    csvwrite([num2str(image_filenumber) '_buffer.csv'], density_buffer_combined)
    dlmwrite([num2str(image_filenumber) '_in.csv'], density_in_combined, 'precision', 10);
    dlmwrite([num2str(image_filenumber) '_buffer.csv'], density_buffer_combined, 'precision', 10);

    
%end
end


