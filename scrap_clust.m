
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
        rand_cell_order = randperm(n);
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
                            %if numel(Neighbours2)>=this_clustsize
                                Neighbours=[Neighbours setdiff(Neighbours2, Neighbours)];   %#ok
                            %end
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
