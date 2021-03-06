% A master script for drawing boundaries around FITS data and checking the
% fit
% this_pool = cbupool(128);
% this_pool.ResourceTemplate = '-l nodes=^N^,mem=256GB,walltime=168:00:00';
% parpool(this_pool,this_pool.NumWorkers)
%dbstop if error

cluster_size = [5000]; % Define the size of the cluster of interest for distance in terms of number of cells

visualise = 0; %Whether or not to check relative scaling of FITS and SVS
image_filenumber = 625951

data = fitsread(['./IT_PT_zone/625951.fits'],'binarytable');
info = fitsinfo(['./IT_PT_zone/625951.fits']);
image_path = ['./IT_PT_zone/625951.svs'];

%%

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

%%
% %Create output file
% outfile = ['./clustering_data_multi_distance_fouth.csv'];
% detail_dir = ['./detail_newlist/'];
% if ~exist(detail_dir,'dir')
% mkdir(detail_dir)
% end
% outfile_detail_stem = [detail_dir 'clustering_detail_'];
% fileID = fopen(outfile,'w');

%Write output file header
% all_combinations = combvec(1:4,1:4); % 0:4 includes rubbish, 1:4 excludes
all_combinations = combvec(1:3, 1:3);% 0:4 includes rubbish, 1:4 excludes
all_cells=combvec(1,1)

key{1} = 'rubbish';
key{2} = 'tumour';
key{3} = 'lymphocyte';
key{4} = 'stroma';
key{5} = 'normal';
header_string = [];
for this_comb_outside = 1:size(all_combinations,2)
    header_string = [header_string ',Av_Mean_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',Av_Bootstrap_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',iqr_Mean_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1} ',iqr_Bootstrap_Distance_' key{all_combinations(1,this_comb_outside)+1} '_to_' key{all_combinations(2,this_comb_outside)+1}];
end

% fprintf(fileID,['Slide_ID,Cluster_Size,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal' header_string '\n']);
% fclose(fileID);
%%
%parfor thisfile = 1:size(filenames,1) %could parallelise here by subject
%for thisfile = 1:size(filenames,1) %could parallelise here by subject
this_comb_inside = 0
    
% Assume that X_global is always the third field, Y_global the fourth
% field, and cell type the final field.

X_ind = 3;
Y_ind = 4;
cell_ind = size(data,2);


data_trimmed = data;
data_trimmed{cell_ind} = data_trimmed{cell_ind}(data{cell_ind}~=0);
data_trimmed{X_ind} = data_trimmed{X_ind}(data{cell_ind}~=0);
data_trimmed{Y_ind} = data_trimmed{Y_ind}(data{cell_ind}~=0);
data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;
data_trimmed{end+1} = ones(size(data_trimmed{cell_ind}, 1),1);


    %assert(max(data{cell_ind})==4&&min(data{cell_ind})==0,['Cell types must span 0-4 and they do not in file ' filenames(thisfile).name])
    % 0 is rubbish
    % 1 is tumour
    % 2 is lymphocyte
    % 3 is stroma
    % 4 is normal
   
 
% First compute the number of each cell type
num_total = size(data_trimmed{cell_ind},1);
num_tum_cells = sum(data_trimmed{cell_ind}==1);
num_ly_cells = sum(data_trimmed{cell_ind}==2);
num_str_cells = sum(data_trimmed{cell_ind}==3);

% Then compute the proportion of each cell type
prop_tum_cells = num_tum_cells/num_total;
prop_ly_cells = num_ly_cells/num_total;
prop_str_cells = num_str_cells/num_total;
    
 

%% 
% Now compute the euclidian distances between each cell type and its
    % nearest neighbour of another cell type
    av_bootstrap_distance = zeros(length(cluster_size),size(all_combinations,2));
    iqr_bootstrap_distance = zeros(length(cluster_size),size(all_combinations,2));
    av_real_distance = zeros(length(cluster_size),size(all_combinations,2));
    iqr_real_distance = zeros(length(cluster_size),size(all_combinations,2));
  
    
     this_comb_inside = 1:size(all_cells,2)
        sprintf(['Working on file ' ' combination ' num2str(this_comb_inside) ])
        base_cells = data_trimmed{end}==all_combinations(1,this_comb_inside);
        neighbour_cells = data_trimmed{end}==all_combinations(2,this_comb_inside);
%         if sum(base_cells)==0||sum(neighbour_cells)==0
%             av_bootstrap_distance(:,this_comb_inside) = NaN;
%             av_real_distance(:,this_comb_inside) = NaN;
%             iqr_bootstrap_distance(:,this_comb_inside) = NaN;
%             iqr_real_distance(:,this_comb_inside) = NaN;
%         else
            
            % First compute the real distances between the cell of interest and its
            % second nearest neighbour (to allow same cell type distance
            % computation)
            %av_real_distance(this_comb_inside) = mean(max(pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',2)));
            %clear all_real_distances all_multi_real_distances
            [all_multi_real_distances, all_indexes] = pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',max(cluster_size)+1);
            i = 0;
            j = 0;
            if size(all_multi_real_distances,1)<max(cluster_size)+1
                all_multi_real_distances(size(all_multi_real_distances,1)+1:max(cluster_size)+1,:)=NaN;
            end
            for this_clustsize = cluster_size
                i = i+1;
                if any(all_multi_real_distances(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                    all_mean_real_distances = mean(all_multi_real_distances(2:(this_clustsize+1),:),1);
                    all_real_distances = all_multi_real_distances(2:(this_clustsize+1),:);
                else
                    all_mean_real_distances = mean(all_multi_real_distances(1:this_clustsize,:),1);
                    all_real_distances = all_multi_real_distances(1:(this_clustsize),:);
                end
                av_real_distance(i,this_comb_inside) = median(all_mean_real_distances);
                iqr_real_distance(i,this_comb_inside) = iqr(all_mean_real_distances);
                
                
                %Now calculate the clusters using DBScan if same cell type
                if all_combinations(1,this_comb_inside) == all_combinations(2,this_comb_inside)
                    if this_clustsize > 1
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
                                Neighbors=setdiff(all_indexes(all_multi_real_distances(:,this_cell)<=epsilon,this_cell),this_cell)';
                                if numel(Neighbors)<this_clustsize
                                    % X(i,:) is NOISE
                                    isnoise(this_cell)=true;
                                else
                                    C=C+1;
                                    cluster_membership{i}(this_cell)=C;
                                    k = 1;
                                    while true
                                        this_neighbour_cell = Neighbors(k);
                                        
                                        if ~visited(this_neighbour_cell)
                                            visited(this_neighbour_cell)=true;
                                            Neighbors2=setdiff(all_indexes(all_multi_real_distances(:,this_neighbour_cell)<=epsilon,this_neighbour_cell),this_neighbour_cell)';
                                            if numel(Neighbors2)>=this_clustsize
                                                Neighbors=[Neighbors Neighbors2];   %#ok
                                            end
                                        end
                                        if cluster_membership{i}(this_neighbour_cell)==0
                                            cluster_membership{i}(this_neighbour_cell)=C;
                                        end
                                        
                                        k = k + 1;
                                        if k > numel(Neighbors)
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
                            this_cluster_boundary_numbers = boundary(this_x_subset,this_y_subset, 1);
                            this_cluster_boundary{i}{this_cluster_number} = [this_x_subset(this_cluster_boundary_numbers),this_y_subset(this_cluster_boundary_numbers)];
                        end
                        figure
                        ax = gca();
                        hold on
                        for this_cluster_number = 1:length(this_cluster_boundary{i})
                            if ~isempty(this_cluster_boundary{i}{this_cluster_number})
                                plot(this_cluster_boundary{i}{this_cluster_number}(:,1),this_cluster_boundary{i}{this_cluster_number}(:,2),'k-')
                            end
                        end
                        hold(ax, 'on');
                        imh = imshow(large_thumbnail_io);
                        hold(ax, 'off');
                        uistack(imh, 'bottom')
                        title(['Minimum cluster size ' num2str(this_clustsize)]);
                        drawnow
                        thisfname = ['./clustfigs/' filenames(thisfile).name '_' num2str(all_combinations(1,this_comb_inside)) 'to' num2str(all_combinations(2,this_comb_inside)) ' minimum cluster size ' num2str(this_clustsize) '.png'];
%                         eval(['export_fig ''' thisfname ''' -transparent'])

                    end
                end
                %                     figure
                %                     ax = gca();
                %                     scatter(data{X_ind}(data{cell_ind}==1),data{Y_ind}(data{cell_ind}==1),1,cluster_membership{2}) %Ignore cell type 0
                %                     hold(ax, 'on');
                %                     imh = imshow(large_thumbnail_io);
                %                     hold(ax, 'off');
                %                     uistack(imh, 'bottom')
                %
                %                     test_cluster = cluster_membership{i};
                %                     test_cluster(test_cluster~=0) = 1;
                %
                %                     figure
                %                     ax = gca();
                %                     scatter(data{X_ind}(data{cell_ind}==1),data{Y_ind}(data{cell_ind}==1),1,test_cluster) %Ignore cell type 0
                %                     hold(ax, 'on');
                %                     imh = imshow(large_thumbnail_io);
                %                     hold(ax, 'off');
                %                     uistack(imh, 'bottom')
                %
                
                
                
                
            end
                      
           
            
        end
    end


