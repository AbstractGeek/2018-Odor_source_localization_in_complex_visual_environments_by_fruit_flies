%% Defaults
clear, close all;
%calibrations
mm_per_pixel_11_july = 0.1241;
mm_per_pixel_12_july = 0.1188;
mm_per_pixel_19_july = 0.1166;
% image details
img_size = [700,900];
background_image_ind = [(1:2), (4:14),(16:25)];
object_centers = csvread(fullfile(pwd,'Center_Digitization','New_object_centers.csv'));
% color
color = [1,0,0];
color_alpha = 0.4;
% load data and set out folder
load('plume_data');
mask_folder = fullfile(pwd,'Analyzed','masked_images');
if ~isdir(mask_folder)
    mkdir(mask_folder);
end
base_folder = fullfile(pwd,'Analyzed','base_images');
if ~isdir(base_folder)
    mkdir(base_folder);
end

%% Load both the raw and the overlaid image
base_image = fullfile(pwd,'Data(Cropped)','Cropped_steady_state_average.tif');
binary_plume = fullfile(pwd,'Data(Cropped)','Cropped_binarized_plume.tif');
base_image_info = imfinfo(base_image);
num_base_images = numel(base_image_info);
binary_plume_info = imfinfo(binary_plume);
num_binary = numel(binary_plume_info);

if num_binary == length(background_image_ind)
    
    for i=1:num_binary
        base = imread(base_image,background_image_ind(i),...
            'Info',base_image_info);
        binary = imread(binary_plume,i,...
            'Info',binary_plume_info);
        obj_center = object_centers(background_image_ind(i),1:2);
        img_start = object_centers(background_image_ind(i),3:4);
        plume_1 = plume_data.(sprintf('x%d',i)).plume_1;
        plume_2 = plume_data.(sprintf('x%d',i)).plume_2;
        
        % Make base rgb
        base_rgb = cat(3, base, base, base);
        % Obtain colored binary
        binary_mask = double(binary>128);
        binary_mask = cat(3, binary_mask.*color(1),...
            binary_mask.*color(2), binary_mask.*color(3));
        
        % Add base and rgb
        masked_image = im2uint8(im2double(base_rgb) + color_alpha.*binary_mask);
        
        % Add discretized plume boundaries
        temp_img = insertMarker(masked_image,obj_center,'o','size',10);
        temp_img = insertShape(temp_img,'Line',...
            [plume_1(1:end-1,:),plume_1(2:end,:)],...
            'LineWidth',3.0);
        temp_img = insertShape(temp_img,'Line',...
            [plume_2(1:end-1,:),plume_2(2:end,:)],...
            'LineWidth',3.0);
        
        % shorten image based on start
        final_img = temp_img(img_start(1):end,img_start(2):end,:);
        final_base_img = base_rgb(img_start(1):end,img_start(2):end,:);
        
        % save image
        imwrite(final_img, fullfile(mask_folder,sprintf('%d.png',i)));
        
        % save the raw image too (with the center)
        imwrite(final_base_img, fullfile(base_folder,sprintf('%d.png',i)));
        
    end
end

%% Create discretized plume video
