load('list_of_matrices.mat');
hold off

matObj = matfile('list_of_matrices.mat');
details = whos(matObj);

for i = 1:length(details)
    all_matrices{i} = [eval([details(i).name '.X1']); eval([details(i).name, '.X2'])]';
end
%plotting the polygon
for i=1:length(all_matrices)
    pg{i}=polyshape(all_matrices{i});
    plot(pg{i});
    hold on;
end

for i = 1:length(all_matrices)
    polyout{i}= polybuffer(pg{i},0.015);
    plot(polyout{i})
    hold on;
end

