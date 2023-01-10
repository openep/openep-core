function show_pointCloud(xyz, pointData)

pcshow(xyz, pointData, 'MarkerSize', 15)

% n_data = length(pointData);
% pointData_range = linspace(min(pointData), max(pointData), 1000);
% cmap = parula(1000);
% pointData_cmap = zeros(n_data, 3);
% for i = 1:n_data
%     for j = 1:1000
%         if pointData_range(j) > pointData(i)
%             pointData_cmap(i, :) = cmap(j, :);
%             break
%         end
%     end
% end
% pcshow(xyz, pointData_cmap);

% n_colours = range(pointData)/0.01;
% cmap = parula(n_data);
% cmap_to_data = zeros(n_colours,1);
% for i=1:n_colours
%     cmap_to_data(i) = min(pointData)+(i-1)*0.01;
% end
% cmapData = zeros(n_data, 3);
% for i = 1:n_data
%     continue
% end
% figure;
% pcshow(xyz);

end