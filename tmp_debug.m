
for i=1:size(tumour_buffer_in, 2) %looking at the cores with tumours
    tbuffer_l_intersection{i}=[];
    tbuffer_l_intersection_area{i}=[];
    tbuffer_l_intersection_lymph_count{i} = [];
    
    for j=1:size(tumour_buffer_in{i}, 2)
        if isempty(tumour_buffer_in{i}{j})
            continue
        end
        for k = 1:size(lymphocyte_polygon{i}, 2)
            if isempty(lymphocyte_polygon{i}{k})
                continue
            end
            tmp = intersect(tumour_buffer_in{i}{j}, lymphocyte_polygon{i}{k});
            if tmp.NumRegions == 0
                continue
            end
            tmp_lymph = inpolygon(data_mini{X_ind}(data_mini{cell_ind}==2), data_mini{Y_ind}(data_mini{cell_ind}==2), tmp.Vertices(:,1), tmp.Vertices(:,2));
            tmp_lymph_count = sum(tmp_lymph);
            tbuffer_l_intersection{i} = [tbuffer_l_intersection{i} tmp];
            tbuffer_l_intersection_area{i} = [tbuffer_l_intersection_area{i} area(tmp)];
            tbuffer_l_intersection_lymph_count{i} = [tbuffer_l_intersection_lymph_count{i} tmp_lymph_count];
            
            if length(tbuffer_l_intersection_area) ~= length(tbuffer_l_intersection_lymph_count)
                pause
            end
              
            
        end
    end
end