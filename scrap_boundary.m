
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
figure % Plot the results for sanity
imshow(large_thumbnail_io)
hold on;
for i = 1:size(grain,2)
    plot(this_boundary{i}(:,2),this_boundary{i}(:,1),'g','LineWidth',3);
end
