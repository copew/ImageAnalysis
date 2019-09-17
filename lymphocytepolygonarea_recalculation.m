%%found a bug in lymphocyte polygon area calculation, now redoing this one small bit


for i = 1:size(lymphocyte_core_list, 2)
    for l = 1:size(lymphocyte_polygon{lymphocyte_core_list(i)}, 2)
        lymphocyte_polygon_area{lymphocyte_core_list(i)}{l} = area(lymphocyte_polygon{lymphocyte_core_list(i)}{l});
    end
end

%create an output file

lymphocyte_polygon_area_combined(cellfun(@isempty, lymphocyte_polygon_area_combined))=[];

%turn into csv files

csvwrite(['./' num2str(image_filenumber) '/' num2str(image_filenumber) '_lymphocyte_area.csv'], lymphocyte_polygon_area_combined);
