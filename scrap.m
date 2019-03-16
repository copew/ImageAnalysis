% trying out drawing boundaries around the cores

data = fitsread(['./IT_PT_zone/625951.fits'],'binarytable');
info = fitsinfo(['./IT_PT_zone/625951.fits']);
image_path = ['./IT_PT_zone/625951.svs'];

X_ind = 3;
Y_ind = 4;
cell_ind = size(data,2)

data_trimmed = data;
    for i = 1:size(data, 2)
        data_trimmed{i} = data_trimmed{i}(data{cell_ind}~=0);
    end
 data_trimmed{end+1} = 1+zeros(size(data_trimmed{1}, 1),1);
    
% and combine tumour and normal
data_trimmed{cell_ind}(data_trimmed{cell_ind} == 4) = 1;

this_clustsize = 50

%%
% getting the thumbnail for later
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
%this is the bit to identify clustering of all cell-type
all_cells = combvec (1, 1)

 % Now compute the euclidian distances between each cell type and its
    % nearest neighbour of another cell type
    av_bootstrap_distance = zeros(length(cluster_size),size(all_cells,2));
    iqr_bootstrap_distance = zeros(length(cluster_size),size(all_cells,2));
    av_real_distance = zeros(length(cluster_size),size(all_cells,2));
    iqr_real_distance = zeros(length(cluster_size),size(all_cells,2));
    
   this_comb_inside = 1
        sprintf(['Working on file ' filenames(thisfile).name ' combination ' num2str(this_comb_inside) ])
        base_cells = data_trimmed{63}==all_cells(1,this_comb_inside);
        neighbour_cells = data_trimmed{63}==all_cells(2,this_comb_inside);
%         if sum(base_cells)==0||sum(neighbour_cells)==0
%             av_bootstrap_distance(:,this_comb_inside) = NaN;
%             av_real_distance(:,this_comb_inside) = NaN;
%             iqr_bootstrap_distance(:,this_comb_inside) = NaN;
%             iqr_real_distance(:,this_comb_inside) = NaN;
%         else
      
        [all_multi_real_distances, all_indexes] = pdist2([data_trimmed{X_ind}(neighbour_cells) data_trimmed{Y_ind}(neighbour_cells)],[data_trimmed{X_ind}(base_cells) data_trimmed{Y_ind}(base_cells)],'euclidean','Smallest',50);
        %calculating 50 nearest distances
        i = 0;
        j = 0;
            if size(all_multi_real_distances,1)< 50
                all_multi_real_distances(size(all_multi_real_distances,1)+1:50,:)=NaN;
            end
            
            i = i+1;
            if any(all_multi_real_distances(1,:)==0) %Now exclude case where zero distance indicates it is the same point
                all_mean_real_distances = mean(all_multi_real_distances(2:50,:),1);
                all_real_distances = all_multi_real_distances(2:50,:);
            else
                all_mean_real_distances = mean(all_multi_real_distances(1:this_clustsize,:),1);
                all_real_distances = all_multi_real_distances(1:(this_clustsize),:);
            end
            av_real_distance(i,this_comb_inside) = median(all_mean_real_distances);
            iqr_real_distance(i,this_comb_inside) = iqr(all_mean_real_distances);
     
                
                %Now calculate the clusters using DBScan if same cell type
                if all_cells(1,this_comb_inside) == all_cells(2,this_comb_inside)
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
                        thisfname = ['./clustfigs/' filenames(thisfile).name '_' num2str(all_cells(1,this_comb_inside)) 'to' num2str(all_cells(2,this_comb_inside)) ' minimum cluster size ' num2str(this_clustsize) '.png'];
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
                      
          