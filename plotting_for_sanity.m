% this is to check the polygons are correct and plot it

%list of images
image_list=[603288, 593987, 619857, 619872, 619905, 625951];

% image = 2 and 6, these are the weird ones
for image = 1:size(image_list, 2)
    image_filenumber=image_list(image);
    
% load image and image info
%data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
image_path = ['./IT_PT_zone/' num2str(image_filenumber) '.svs'];

% load the polygons
load(['./' num2str(image_filenumber) '/tumour_polygon_in.mat'])
load(['./' num2str(image_filenumber) '/lymphocyte_polygon.mat'])

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


%% plotting

figure
ax=gca();
hold on;
imshow(large_thumbnail_io);
hold(ax, 'on');
% draw the tumour cluster first
for i = 1:size(tumour_polygon_in, 2)
    if isempty(tumour_polygon_in{i}) == 1
        continue
    end
    for j = 1:size(tumour_polygon_in{i}, 2)
        % hold on;
        if isempty(tumour_polygon_in{i}{j}) == 1
            continue
        end
        plot(tumour_polygon_in{i}{j}, 'Facecolor', 'blue'); %,'k-');
    end
    hold on;
end
    
% now draw the lymphocyte cluster
for i = 1:size(lymphocyte_polygon, 2)
    if isempty(lymphocyte_polygon{i}) == 1
        continue
    end
    for j = 1:size(lymphocyte_polygon{i}, 2)
        hold on;
        if isempty(lymphocyte_polygon{i}{j}) == 1
            continue
        end
        plot(lymphocyte_polygon{i}{j}, 'Facecolor', 'red'); %,'r-');
    end
    hold on;
end
hold(ax, 'off')

%
end

