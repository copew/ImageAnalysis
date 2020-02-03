% this is to calculate pair wise distances between different cell types
%
%this is a total list
%image_list =[603283,603279,603271,603262,603269,603263,603257,603253,602994,603006,603245,603007,603055,602983,602976,602971,602962,602952,602958,602966,602945,602951,602942,602935,602927,602921,602915,593960,593708,593978,593971,593987,594006,594000,594017,594035,594044,594051,594060,594080,594088,594095,594105,594116,594110,595806,594137,597777,597768,597794,597800,597812,597820,597786,597830,597838,599803,599797,599831,599850,599813,604094,603952,604089,604098,604083,607165,603944,604078,604057,607166,604071,604077,603911,603919,604105,603925,603922,603928,603941,603904,604061,604066,604055,604047,603987,604045,603976,603982,603963,603970,603966,603956,619678,619884,619809,619942,619868,619556,619896,619937,619926,619922,619932,619916,619911,619533,619852,619863,619859,619842,619837,619832,619847,619583,619815,619803,619799,619793,619579,619571,619538,619550,619877,619539,619460,619954,619465,619470,619476,619483,619489,619503,619508,619496,619515,619527,625891,625322,625333,625338,626172,625865,626171,625876,625887,625895,625908,625911,619953,625916,626166,625923,625930,625936,625958,625946,625951,626018,626047,626102,626103,626160, 603288, 603298,603292, 619872, 619857, 619905];

%%%%% this is not done!!!! %%%%%% 3/2/2020

image_list = [603283,603279,603271];

%% for debugging this can be run for loading
for image = 1:size(image_list,2)
    image_filenumber = image_list(image);
    image_path_stem = '/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/MATLAB/ImageAnalysis';
    
    data = load([image_path_stem '/mat_file_new/' num2str(image_filenumber) '.mat']);
    % % % %     image_path = [image_path_stem '/IT_PT_zone/' num2str(image_filenumber) '.svs'];
    
    load(['/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/tum_lymph_overlap/workspace1/' num2str(image_filenumber) '_workspace1.mat']);
    
    
    %% load images and fits files
    
    % % % % % % % % function tum_lymph_cluster_hpc_diff_parameters(image_filenumber_fullpath)
    % % % % % % % % disp(image_filenumber_fullpath)
    % % % % % % % % disp(class(image_filenumber_fullpath))
    % % % % % % % % if ischar(image_filenumber_fullpath)
    % % % % % % % %     [image_path_stem,image_filenumber,~] = fileparts(image_filenumber_fullpath);
    % % % % % % % %     image_filenumber = str2num(image_filenumber);
    % % % % % % % %     if isempty(image_path_stem)
    % % % % % % % %         image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt';
    % % % % % % % %         warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
    % % % % % % % %     end
    % % % % % % % % elseif isnumeric(image_filenumber_fullpath)
    % % % % % % % %     image_filenumber = image_filenumber_fullpath;
    % % % % % % % %     image_path_stem = '/rds-d4/user/ww234/hpc-work/itpt';
    % % % % % % % %     warning('The input did not give a full path so assuming the data are in /rds-d4/user/ww234/hpc-work/itpt')
    % % % % % % % % end
    % % % % % % % %
    % % % % % % % % %%%loading files
    % % % % % % % % % % % data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
    % % % % % % % % % % % info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
    % % % % % % % % % % % image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];
    % % % % % % % % % %
    % % % % % % % % data= load([image_path_stem '/' num2str(image_filenumber) '.mat']);
    % % % % % % % % %image_path = [image_path_stem '/' num2str(image_filenumber) '.svs'];
    
    %%%%%%%%%%%%%%%% comment from function to here for debugging  %%%%%%%%%%%%%%%
    
    
    
    % make a folder to save files
    mkdir(num2str(image_filenumber))
    
    
    %% loading up data
    
    %trim data
    data_trimmed = data.data_mini;
    data_trimmed = num2cell(data_trimmed, 1); %so that the cell index works later on
    
    %removed non cells, overlaps, low signal to noise ratio (set threshold at 1.3 - discussed with A Dariush) and combine tumour and normal
    
    
    
    
    % set parameters
    
    %now the cluster size is the number of measurement that we are taking
    cluster_size = [1000];
    
    
    % and index
    X_ind = 1;
    Y_ind = 2;
    cell_ind = 3;
    
    
    
    %% Now compute the tumour and lymphocyte clusters with subsetting
    
    core_list = [];
    core_list=[1:1:size(core_polygon,2)];
    
    
    for this_core = 1:size(core_list,2)
        
        %tumour - tumour
        
        base_cells = in_core{this_core}{cell_ind}==1;
        neighbour_cells = in_core{this_core}{cell_ind}==1;
        
        [all_multi_real_distances_tt{this_core}, all_indexes_tt{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
        
        
        
        % lymphocyte - lymphocyte
        
        base_cells = in_core{this_core}{cell_ind}==2;
        neighbour_cells = in_core{this_core}{cell_ind}==2;
        
        [all_multi_real_distances_ll{this_core}, all_indexes_ll{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
        
        
        
        % stroma - stroma
        
        base_cells = in_core{this_core}{cell_ind}==3;
        neighbour_cells = in_core{this_core}{cell_ind}==3;
        
        [all_multi_real_distances_ll{this_core}, all_indexes_ll{this_core}] = pdist2([in_core{this_core}{X_ind}(neighbour_cells) in_core{this_core}{Y_ind}(neighbour_cells)],[in_core{this_core}{X_ind}(base_cells) in_core{this_core}{Y_ind}(base_cells)],'euclidean','Smallest', cluster_size+1);
        
        
        
        
        
        
        
        
        
    end
end



