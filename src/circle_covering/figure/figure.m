figure('Color', 'w');
hold on;
% load container
if exist("container.csv",'file')
    disp('Container found:');
    container = csvread('container.csv');
    x_container = container(:,1); 
    y_container = container(:,2);
    % plot(x_container,y_container,'black');
    fill(x_container, y_container, [255,196,112]/255, 'EdgeColor', 'k','LineWidth',1.5,'FaceAlpha',1);
end
poly_container = polyshape(x_container,y_container);

% load voronoi
if exist("voronoi",'file')
    folderPath = 'voronoi';
    fileList = dir(fullfile(folderPath, '*.csv'));
end

for i  = 1:length(fileList)
    hold on;
    filePath = fullfile(folderPath, fileList(i).name);
    data = csvread(filePath);
    x_point_data = data(1,1);
    y_point_data = data(1,2);
    % scatter(x_point_data, y_point_data,'.', 'MarkerEdgeColor', [0, 0, 0]);
    scatter(x_point_data, y_point_data,'.');
    disp(i);
    disp(x_point_data);
    x_data = data(2:end,1);
    y_data = data(2:end,2);
    poly = polyshape(x_data,y_data);
    poly_intersect = intersect(poly,poly_container);
    % plot(poly_intersect,'FaceColor',rand(1,3),'FaceAlpha',0.1);
    plot(poly_intersect,'FaceColor','none',LineStyle=':');
end

%load holes
hole_files = dir(fullfile("holes", '*.csv'));
if ~isempty(hole_files)
    disp('Hole found.');
    for i = 1:length(hole_files)
        hole_file=fullfile("holes", hole_files(i).name);
        hole = csvread(hole_file)
        x_hole = hole(:,1);
        y_hole = hole(:,2); 
        hole_plot = fill(x_hole,y_hole,'w','FaceAlpha',1,'LineWidth',1.5,"EdgeColor",'k');
    end
else
    disp('No Hole found.');
end

points = csvread('initial_points.csv');
scatter(points(:, 1), points(:, 2),100,'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

circle_data = csvread('Results.csv');
circle_centers = circle_data(:, 1:2);
circle_radii = circle_data(:,3);
% scatter(circle_centers(:, 1), circle_centers(:, 2), '.', 'MarkerFaceColor', 'green');

for i = 1:size(circle_centers, 1)
    theta = linspace(0, 2*pi, 100);
    x = circle_centers(i, 1) + circle_radii(i) * cos(theta);
    y = circle_centers(i, 2) + circle_radii(i) * sin(theta);
    patch(x, y, [0.5, 0.5, 0.5], 'EdgeColor', [0, 0.5, 1], 'LineWidth', 1, 'FaceAlpha', 0.1);
end
axis equal;
hold off;