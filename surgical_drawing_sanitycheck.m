% this script is to check that the annotations from qupath is recognised/loaded correctly
% onto matlab before further downstream analysis
%599871

load('/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Neoadjuvant/R_scripts/testlist7.mat')
testlist7 = struct2cell(testlist7);
testlist7 = reshape(testlist7, [1,2]); 





%load image
image_filenumber = 599871;

% loading files
%data = fitsread(['./IT_PT_zone/' num2str(image_filenumber) '.fits'],'binarytable');
%info = fitsinfo(['./IT_PT_zone/' num2str(image_filenumber) '.fits']);
image_path = ['./surgical_test/' num2str(image_filenumber) '.svs'];


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


%plotting
figure % Plot the results for sanity
imshow(large_thumbnail_io)
hold on;
% 
for i = 1:size(testlist7,2)
    plot(testlist7{i}(:,1),testlist7{i}(:,2),'g','LineWidth',1);
end