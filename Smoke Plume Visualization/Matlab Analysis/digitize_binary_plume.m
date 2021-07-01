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
% digitization details
plume_window = 75;
half_window = round((plume_window-1)/2);
start_offset = 6;   % mm
% save folder
digitized_folder = fullfile(pwd,'Analyzed','digitized_images');
if ~isdir(digitized_folder)
    mkdir(digitized_folder);
end

%% Discretize the plume
binary_plume = fullfile(pwd,'Data(Cropped)','Cropped_binarized_plume.tif');
binary_plume_info = imfinfo(binary_plume);
num_images = numel(binary_plume_info);
plume_data = struct;

for i=1:num_images
    binary = imread(binary_plume,i,'Info',binary_plume_info);
    obj_center = object_centers(background_image_ind(i),1:2);
    
    % starting x point
    if i==1
        x_curr = obj_center(1) + round(start_offset * 1/mm_per_pixel_11_july);
    elseif i<14
        x_curr = obj_center(1) + round(start_offset * 1/mm_per_pixel_12_july);
    else 
        x_curr = obj_center(1) + round(start_offset * 1/mm_per_pixel_19_july);
    end
    
    % starting y point
    y_curr_1 = obj_center(2) - half_window - 1;
    y_curr_2 = obj_center(2) + half_window + 1;
    y_init_1 = y_curr_1;
    y_init_2 = y_curr_2;
    
    % Create plume 
    x_values = (x_curr:1:size(binary,2))';
    plume_1 = nan(length(x_values),2);
    plume_2 = nan(length(x_values),2);
    
    for j=1:length(x_values)
        % find both plume boundaries
        plume_1_curr = round(median(find(binary(y_curr_1-half_window:y_curr_1+half_window,x_values(j,1))>128)));
        plume_2_curr = round(median(find(binary(y_curr_2-half_window:y_curr_2+half_window,x_values(j,1))>128)));
        
        % append to the plume trajectory
        if ~isnan(plume_1_curr)
            plume_1(j,1) = x_values(j,1);
            plume_1(j,2) = plume_1_curr + y_curr_1-half_window - 1;
            if (plume_1(j,2) < y_init_1)
                y_curr_1 = plume_1(j,2);
            end
        end
        % append to the plume trajectory
        if ~isnan(plume_2_curr)
            plume_2(j,1) = x_values(j,1);
            plume_2(j,2) = plume_2_curr + y_curr_2-half_window - 1;
            if (plume_2(j,2) > y_init_2)
                y_curr_2= plume_2(j,2);
            end
        end
                
    end
    
    % save the digitized plume into a structure
    plume_data.(sprintf('x%d',i)).raw_plume_1 = plume_1;
    plume_data.(sprintf('x%d',i)).raw_plume_2 = plume_2;
    
    % Obtain the non nan data
    plume_1_without_nan = plume_1(~any(isnan(plume_1),2),:);
    plume_2_without_nan = plume_2(~any(isnan(plume_2),2),:);
    
    % Save as a temporary tiff
    temp_img = insertMarker(binary,obj_center,'o','size',10);
    temp_img = insertShape(temp_img,'Line',...
        [plume_1_without_nan(1:end-1,:),plume_1_without_nan(2:end,:)],...
        'LineWidth',3.0);
    temp_img = insertShape(temp_img,'Line',...
        [plume_2_without_nan(1:end-1,:),plume_2_without_nan(2:end,:)],...
        'LineWidth',3.0);    
    
    % Find out the start and the end indexes
    plume_1_vertices = [find(diff(isnan(plume_1(:,2)))==-1,1,'first')+1,...
        size(plume_1(:,2),1)];
    plume_2_vertices = [find(diff(isnan(plume_2(:,2)))==-1,1,'first')+1,...
        size(plume_2(:,2),1)];
    
    if length(plume_1_vertices)==1
        plume_1_vertices = [1,plume_1_vertices];
    end
    
    if length(plume_2_vertices)==1
        plume_2_vertices = [1,plume_2_vertices];
    end
    
    if i==1
        imwrite(temp_img,...
            fullfile(digitized_folder, ...
            'Raw_binarized_plume_with_center_marked.tif'));                
    else
        imwrite(temp_img,...
            fullfile(digitized_folder, ...
            'Raw_binarized_plume_with_center_marked.tif'),...
            'WriteMode','append');        
    end
    
    % spline approximate the data
    plume_1 = [x_values(plume_1_vertices(1):plume_1_vertices(2)),...
        pchip(plume_1_without_nan(:,1),... %piecewise cubic interpolating polynomial
        plume_1_without_nan(:,2),...
        x_values(plume_1_vertices(1):plume_1_vertices(2)))];
    plume_2 = [x_values(plume_2_vertices(1):plume_2_vertices(2)),...
        pchip(plume_2_without_nan(:,1),...
        plume_2_without_nan(:,2),...
        x_values(plume_2_vertices(1):plume_2_vertices(2)))];
    
    % save the digitized plume into a structure
    plume_data.(sprintf('x%d',i)).plume_1 = plume_1;
    plume_data.(sprintf('x%d',i)).plume_2 = plume_2;
    
    % convert to mm and add image
    if i==1
        odor_plume = [x_values(1:plume_1_vertices(2)),...
            [nan(plume_1_vertices(1)-1,1);...
                plume_1(:,2)],...
            [nan(plume_2_vertices(1)-1,1);...
                plume_2(:,2)]];    
    else        
        odor_plume = [x_values,...
            [nan(plume_1_vertices(1)-1,1);...
                plume_1(:,2)],...
            [nan(plume_2_vertices(1)-1,1);...
                plume_2(:,2)]];    
    end
    
    plume_data.(sprintf('x%d',i)).odor_plume = odor_plume;
    
    % Save digitized image
    temp_img = insertMarker(binary,obj_center,'o','size',10);
    temp_img = insertShape(temp_img,'Line',...
        [plume_1(1:end-1,:),plume_1(2:end,:)],...
        'LineWidth',3.0,'Color','red');
    temp_img = insertShape(temp_img,'Line',...
        [plume_2(1:end-1,:),plume_2(2:end,:)],...
        'LineWidth',3.0,'Color','red');    
    imwrite(temp_img,fullfile(digitized_folder, ...
            sprintf('%d.png',i)));          
    
end

save('plume_data','plume_data');