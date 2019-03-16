
filenames = dir('./IT_PT_zone/*.fits');

%Create output file
outfile = ['./clustering_data.csv'];
fileID = fopen(outfile,'w');

%Write output file header
all_combinations = combvec(0:4,0:4);
key{1} = 'rubbish';
key{2} = 'tumour';
key{3} = 'lymphocyte';
key{4} = 'stroma';
key{5} = 'normal';
header_string = [];

fprintf(fileID,['Slide_ID,Cluster_Size,Num_Total,Num_Rubbish,Num_Tumour,Num_Lymphs,Num_Stroma,Num_Normal,Prop_Rubbish,Prop_Tumour,Prop_Lymphs,Prop_Stroma,Prop_Normal' header_string '\n']);
fclose(fileID);


for thisfile = 1:size(filenames,1) %could parallelise here by subject
   
    try
        sprintf(['Working on file ' filenames(thisfile).name])
        
        data = fitsread(['./IT_PT_zone/' filenames(thisfile).name],'binarytable');
        info = fitsinfo(['./IT_PT_zone/' filenames(thisfile).name]);
        
        % Assume that X_global is always the third field, Y_global the fourth
        % field, and cell type the final field.
        
        X_ind = 3;
        Y_ind = 4;
        cell_ind = size(data,2);
        
        %assert(max(data{cell_ind})==4&&min(data{cell_ind})==0,['Cell types must span 0-4 and they do not in file ' filenames(thisfile).name])
        % 0 is rubbish
        % 1 is tumour
        % 2 is lymphocyte
        % 3 is stroma
        % 4 is normal
        
        % %Optionally visualise the slide
        % figure
        % scatter(data{X_ind}(data{cell_ind}~=0),data{Y_ind}(data{cell_ind}~=0),1,data{cell_ind}(data{cell_ind}~=0)) %Ignore cell type 0
        
        % First compute the number of each cell type
        num_total = size(data{cell_ind},1);
        num_rubbish = sum(data{cell_ind}==0);
        num_tum_cells = sum(data{cell_ind}==1);
        num_ly_cells = sum(data{cell_ind}==2);
        num_str_cells = sum(data{cell_ind}==3);
        num_norm_cells = sum(data{cell_ind}==4);
        
        % Then compute the proportion of each cell type
        prop_rubbish = num_rubbish/num_total;
        prop_tum_cells = num_tum_cells/num_total;
        prop_ly_cells = num_ly_cells/num_total;
        prop_str_cells = num_str_cells/num_total;
        prop_norm_cells = num_norm_cells/num_total;
        
        % Now compute the euclidian distances between each cell type and its
        % nearest neighbour of another cell type
        
        data_string = [];
        
        
        outfile = ['./clustering_data.csv'];
        fileID = fopen(outfile,'a');
        fprintf(fileID,[filenames(thisfile).name(1:end-5) ',' num2str(cluster_size) ',' num2str(num_total) ',' num2str(num_rubbish) ',' num2str(num_tum_cells) ',' num2str(num_ly_cells) ',' num2str(num_str_cells) ',' num2str(num_norm_cells) ',' num2str(prop_rubbish) ',' num2str(prop_tum_cells) ',' num2str(prop_ly_cells) ',' num2str(prop_str_cells) ',' num2str(prop_norm_cells) data_string '\n']);
        fclose(fileID);
        sprintf(['Finished file ' filenames(thisfile).name])
    catch
        outfile = ['./clustering_data.csv'];
        fileID = fopen(outfile,'a');
        fprintf(fileID,[filenames(thisfile).name(1:end-5) ',failed at ' num2str(this_comb) '\n']);
        fclose(fileID);
        sprintf(['Failed file ' filenames(thisfile).name ' moving on'])
    end
end

%end