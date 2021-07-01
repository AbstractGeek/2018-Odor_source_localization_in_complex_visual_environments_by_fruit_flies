%% Defaults
clear, close all;
% calibration
mm_per_pixel_11_july = 0.1241;
mm_per_pixel_12_july = 0.1188;
mm_per_pixel_19_july = 0.1166;
% image details
img_size = [700,900];
background_image_ind = [(1:2), (4:14),(16:25)];
% source and destination
object_center_csv_name = fullfile(pwd,'Center_Digitization',...
    'object_location_xypts.csv');
object_centers = csvread(object_center_csv_name,1,0);
cropped_folder = fullfile(pwd,'Analyzed','cropped_images_extra');
if ~isdir(cropped_folder)
    mkdir(cropped_folder);
end

%% Load Steady state average
base_image = fullfile('Data(Raw)','Steady_state_average_combined.tif');
base_image_info = imfinfo(base_image);
num_base_images = numel(base_image_info);
new_object_centers = nan(size(object_centers,1),4);

for i=1:num_base_images
    % Load image
    base = imread(base_image,i,'Info',base_image_info);
    
    % Obtain center
    curr_center = round([object_centers(i,1),size(base,1) - object_centers(i,2)]);
    
    % Obtain 20 mm and 80 mm distance from the center
    if i==1
        pixels_20mm = ceil(20 * 1/mm_per_pixel_11_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_11_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_11_july);
    elseif (i<16)
        pixels_20mm = ceil(20 * 1/mm_per_pixel_12_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_12_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_12_july);
    else 
        pixels_20mm = ceil(20 * 1/mm_per_pixel_19_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_19_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_19_july);
    end
    
    % Obtain new corners from the data
    left_corner = [curr_center(1)-pixels_20mm, curr_center(2)-pixels_40mm];
    right_corner = [curr_center(1)+pixels_80mm, curr_center(2)+pixels_40mm];
    
    % Obtain new image
    base_new = base(left_corner(2):right_corner(2),left_corner(1):right_corner(1));
    
    % Save it in the template
    new_start = (img_size - size(base_new)) + [1,1];
    base_fixed = uint8(zeros(img_size));
    base_fixed(new_start(1):end,new_start(2):end) = base_new;
    
    % Save new center
    new_center = [new_start(2) + pixels_20mm, new_start(1) + pixels_40mm];
    new_object_centers(i,1:2) = new_center;
    new_object_centers(i,3:4) = new_start;
    
    % Save images as tiff
    if i==1        
        % Write in
        imwrite(base_fixed,...
            fullfile('Data(Cropped)','Cropped_steady_state_average.tif'),...
            'tiff', 'Compression','none');
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Steady_state_average_with_center_marked.tif'),...
            'tiff', 'Compression','none');
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_steady_state_average_with_center_marked.tif'),...
            'tiff', 'Compression','none');        
    else
        imwrite(base_fixed,...
            fullfile('Data(Cropped)','Cropped_steady_state_average.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Steady_state_average_with_center_marked.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_steady_state_average_with_center_marked.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
    end
    
end


%% Load background subtracted averages
base_image = fullfile('Data(Raw)','Background_subtracted_average_combined.tif');
base_image_info = imfinfo(base_image);
num_base_images = numel(base_image_info);
new_object_centers = nan(size(object_centers,1),4);

for i=1:num_base_images
    % Load image
    base = imread(base_image,i,'Info',base_image_info);
    
    % Obtain center
    curr_center = round([object_centers(i,1),size(base,1) - object_centers(i,2)]);
    
    % Obtain 20 mm and 80 mm distance from the center
    if i==1
        pixels_20mm = ceil(20 * 1/mm_per_pixel_11_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_11_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_11_july);
    elseif (i<16)
        pixels_20mm = ceil(20 * 1/mm_per_pixel_12_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_12_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_12_july);
    else 
        pixels_20mm = ceil(20 * 1/mm_per_pixel_19_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_19_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_19_july);
    end
    
    % Obtain new corners from the data
    left_corner = [curr_center(1)-pixels_20mm, curr_center(2)-pixels_40mm];
    right_corner = [curr_center(1)+pixels_80mm, curr_center(2)+pixels_40mm];
    
    % Obtain new image
    base_new = base(left_corner(2):right_corner(2),left_corner(1):right_corner(1));
    
    % Save it in the template
    new_start = (img_size - size(base_new)) + [1,1];
    base_fixed = uint8(zeros(img_size));
    base_fixed(new_start(1):end,new_start(2):end) = base_new;
    
    % Save new center
    new_center = [new_start(2) + pixels_20mm, new_start(1) + pixels_40mm];
    new_object_centers(i,1:2) = new_center;
    new_object_centers(i,3:4) = new_start;
    
    % Save images as tiff
    if i==1        
        % Write in
        imwrite(base_fixed,...
            fullfile('Data(Cropped)','Cropped_background_subtracted_average.tif'),...
            'tiff', 'Compression','none');
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Background_subtracted_average_with_center_marked.tif'),...
            'tiff', 'Compression','none');
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_background_subtracted_average_with_center_marked.tif'),...
            'tiff', 'Compression','none');        
    else
        imwrite(base_fixed,...
            fullfile('Data(Cropped)','Cropped_background_subtracted_average.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Background_subtracted_average_with_center_marked.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_background_subtracted_average_with_center_marked.tif'),...
            'tiff','WriteMode','append', 'Compression','none');
    end
    
end

csvwrite(fullfile(pwd,'Center_Digitization','New_object_centers.csv'),...
    new_object_centers);



%% Load binary plume
binary_plume = fullfile('Data(Raw)','Binarized_plume_combined.tif');
binary_plume_info = imfinfo(binary_plume);
num_images = numel(binary_plume_info);

for i=1:num_images
    % Load image
    base = imread(binary_plume,i,'Info',binary_plume_info);
    
    % Obtain center
    curr_center = round([object_centers(background_image_ind(i),1),...
        size(base,1) - object_centers(background_image_ind(i),2)]);
    
    % Obtain 20 mm and 80 mm distance from the center
    if i==1
        pixels_20mm = ceil(20 * 1/mm_per_pixel_11_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_11_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_11_july);
    elseif (i<14)
        pixels_20mm = ceil(20 * 1/mm_per_pixel_12_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_12_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_12_july);
    else 
        pixels_20mm = ceil(20 * 1/mm_per_pixel_19_july);
        pixels_40mm = ceil(40 * 1/mm_per_pixel_19_july);
        pixels_80mm = ceil(80 * 1/mm_per_pixel_19_july);
    end
    
    % Obtain new corners from the data
    left_corner = [curr_center(1)-pixels_20mm, curr_center(2)-pixels_40mm];
    right_corner = [curr_center(1)+pixels_80mm, curr_center(2)+pixels_40mm];
    
    % Obtain new image
    base_new = base(left_corner(2):right_corner(2),left_corner(1):right_corner(1));
     
    % Save it in the template
    new_start = new_object_centers(background_image_ind(i),3:4);
    base_fixed = uint8(zeros(img_size));
    base_fixed(new_start(1):end,new_start(2):end) = base_new;
    
    % Obtain new center
    new_center = new_object_centers(background_image_ind(i),1:2);
    
    % Save images as tiff
    if i==1        
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Binarized_plume_with_center_marked.tif'));
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_binarized_plume_with_center_marked.tif'));
        imwrite(base_fixed,...
            fullfile('Data(Cropped)','Cropped_binarized_plume.tif'),...
            'tiff', 'Compression','none');
    else
        imwrite(insertMarker(base,curr_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Binarized_plume_with_center_marked.tif'),...
            'WriteMode','append');
        imwrite(insertMarker(base_fixed,new_center,'o','size',10),...
            fullfile(cropped_folder, ...
            'Cropped_binarized_plume_with_center_marked.tif'),...
            'WriteMode','append');
        imwrite(base_fixed,...            
            fullfile('Data(Cropped)','Cropped_binarized_plume.tif'),...
            'tiff', 'WriteMode','append', 'Compression','none');
    end
    
end

%csvwrite('New_object_centers.csv',new_object_centers);