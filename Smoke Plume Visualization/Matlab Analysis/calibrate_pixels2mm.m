%% Defaults
clear, close all;
side_length = 24; % mm

%% 11th July
% Import calibration
calibration = imread(fullfile(pwd,'Calibration',...
    'calibration_image_11th_July.tif'));
[imagePoints,~] = detectCheckerboardPoints(calibration);
% Plot detected points
imshow(calibration);
axis ij
hold on, plot(imagePoints(:,1),imagePoints(:,2),'or');

% Sort the points
xpts = reshape(imagePoints(:,1),4,size(imagePoints,1)/4);
ypts = reshape(imagePoints(:,2),4,size(imagePoints,1)/4);

dist_pixel = sqrt(diff(xpts).^2 + diff(ypts).^2);
pixel_length = mean(mean(dist_pixel));

mm_per_pixel_11_july = side_length/pixel_length;
fprintf('mm per pixel for the 11th july is: %f \n', mm_per_pixel_11_july);

saveas(gcf, fullfile(pwd,'Calibration',...
    'calibration_image_11th_July_detected_points.tif'));
close(gcf);

%% 12th July
% Import calibration
calibration = imread(fullfile(pwd,'Calibration',...
    'calibration_image_12th_July.tif'));
[imagePoints,~] = detectCheckerboardPoints(calibration);
% Plot detected points
imshow(calibration);
axis ij
hold on, plot(imagePoints(:,1),imagePoints(:,2),'or');

% Sort the points
xpts = reshape(imagePoints(:,1),4,size(imagePoints,1)/4);
ypts = reshape(imagePoints(:,2),4,size(imagePoints,1)/4);

dist_pixel = sqrt(diff(xpts).^2 + diff(ypts).^2);
pixel_length = mean(mean(dist_pixel));

mm_per_pixel_12_july = side_length/pixel_length;
fprintf('mm per pixel for the 12th july is: %f \n', mm_per_pixel_12_july);

saveas(gcf, fullfile(pwd,'Calibration',...
    'calibration_image_12th_July_detected_points.tif'));
close(gcf);

%% 19th July
% Import calibration
calibration = imread(fullfile(pwd,'Calibration',...
    'calibration_image_19th_July.tif'));
[imagePoints,~] = detectCheckerboardPoints(calibration);
% Plot detected points
imshow(calibration);
axis ij
hold on, plot(imagePoints(:,1),imagePoints(:,2),'or');

% Sort the points
xpts = reshape(imagePoints(:,1),4,size(imagePoints,1)/4);
ypts = reshape(imagePoints(:,2),4,size(imagePoints,1)/4);

dist_pixel = sqrt(diff(xpts).^2 + diff(ypts).^2);
pixel_length = mean(mean(dist_pixel));

mm_per_pixel_19_july = side_length/pixel_length;
fprintf('mm per pixel for the 19th july is: %f \n', mm_per_pixel_19_july);

saveas(gcf, fullfile(pwd,'Calibration',...
    'calibration_image_19th_July_detected_points.tif'));
close(gcf);

