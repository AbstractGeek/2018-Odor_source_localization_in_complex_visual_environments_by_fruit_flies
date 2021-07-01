%% Defaults
clear, close all;

% Load plume data
load('plume_data');

% Categorize treatments
treatment.incense_stick = [1];
treatment.bead = [2,20,21,22];
treatment.capillary_bead_1cm = [3,4,14,15];
treatment.capillary_bead_2cm = [16,17,18,19];
treatment.capillary = [5,6,7,23];
treatment.capillary_beads = [8,9];
treatment.beads = [10,11,12,13];

% Calibrate
mm_per_pixel_11_july = 0.1241;
mm_per_pixel_12_july = 0.1188;
mm_per_pixel_19_july = 0.1166;
mm_per_pixel = [mm_per_pixel_11_july;...
    mm_per_pixel_12_july.*ones(12,1);...
    mm_per_pixel_19_july.*ones(10,1);];

% Plotting colors (mean, stderr plots)
line_color = [0.8, 0.8, 0.8];
mean_color = [0, 0.45, 0.74];
std_error_color = [0, 0.45, 0.74];

% Plotting colors (combined plot)
colors = [rgb('Amethyst');...
    rgb('Blue');...
    rgb('Navy');...
    rgb('Sunflower');...
    rgb('Red')];

% image details
object_centers = csvread(fullfile(pwd,'Center_Digitization','New_object_centers.csv'));
background_image_ind = [(1:2), (4:14),(16:25)];

% save folder
plots_folder = fullfile(pwd,'Analyzed','plots');
if ~isdir(plots_folder)
    mkdir(plots_folder);
end

% create a table for mean data and fill it with plume distance for now
mean_odor_plume = table((1:100)','VariableNames',{'SourceDistance'});
stderr_odor_plume = table((1:100)','VariableNames',{'SourceDistance'});

%% Plot mean and standard deviation plots

trts = fieldnames(treatment);

for i=1:length(trts)
    inds = treatment.(trts{i});
    plume_distance_width = nan(100,length(inds)+1);
    plume_distance_width(:,1) = (1:100)';
    
    for j=1:length(inds)
        odor_plume = plume_data.(sprintf('x%d',inds(j))).odor_plume;        
        
        obj_center = object_centers(background_image_ind(inds(j)),1:2);
    
        curr_plume = round([odor_plume(:,1)-obj_center(1),...
            odor_plume(:,3) - odor_plume(:,2)]);
        
        % Convert pixel to mm
        curr_plume = curr_plume .* mm_per_pixel(inds(j));
        
        % Run through the unique dataset
        curr_plume(:,1) = round(curr_plume(:,1));
        unique_distances = unique(curr_plume(:,1));
        
        for k=1:length(unique_distances)
            % save the mean of the data
            plume_distance_width(unique_distances(k),j+1) = ...
                nanmean(curr_plume(curr_plume(:,1)==unique_distances(k),2));
        end
        
    end
    
    % Obtain mean and standard deviation
    mean_plume_width = nanmean(plume_distance_width(:,2:end),2);
    std_error_temp = (nanstd(plume_distance_width(:,2:end),0,2)./sqrt(size(plume_distance_width,2)-1));
    
    % Obtain indexes of interest
    imp_ind = ~isnan(mean_plume_width);
    plume_distance = plume_distance_width(imp_ind,1);
    
    % Obtain standard error flanks
    std_error_top = mean_plume_width(imp_ind) + std_error_temp(imp_ind);
    std_error_down = mean_plume_width(imp_ind) - std_error_temp(imp_ind);
    
    % Plot mean and standard deviation
    figure,hold on;
    for k=2:size(plume_distance_width,2)
        curr_ind = ~isnan(plume_distance_width(:,k));
        plot(plume_distance_width(curr_ind,1),...
            plume_distance_width(curr_ind,k),...
            'Color', line_color);
    end

    fill([plume_distance;plume_distance(end:-1:1)],...
        [std_error_top; std_error_down(end:-1:1)], std_error_color,...
        'Edgecolor', 'none', 'FaceAlpha', 0.25);
    plot(plume_distance, mean_plume_width(imp_ind), 'Color', mean_color);
    
    % Add labels
    xlabel('Source Distance (mm)');
    ylabel('Plume Width (mm)');
    
    % Save figure
    saveas(gcf, fullfile(plots_folder,...
        sprintf('%s.fig',trts{i})));
    saveas(gcf, fullfile(plots_folder,...
        sprintf('%s.png',trts{i})));
    close(gcf);
    
    % Save mean plume width
    mean_odor_plume = [mean_odor_plume,...
        table(mean_plume_width,'VariableNames',trts(i))];   % Save
    stderr_odor_plume = [stderr_odor_plume,...
        table(std_error_temp,'VariableNames',trts(i))];   % Save
    
end

save('odor_plume_data','plume_data','mean_odor_plume','stderr_odor_plume');

%% Plot all categories
treatments = {'capillary','capillary_bead_1cm','capillary_bead_2cm',...
    'bead','beads'};
figure, hold on;
for i=1:length(treatments)    
    % Obtain source distance and plume width
    source_distance = mean_odor_plume{:,1};
    plume_width = mean_odor_plume.(treatments{i});
    plume_error = stderr_odor_plume.(treatments{i});
    % Get not nans
    not_nans = ~isnan(plume_width);
    source_distance = source_distance(not_nans);
    plume_width = plume_width(not_nans);
    plume_error = plume_error(not_nans);
    std_error_top = plume_width + plume_error;
    std_error_down = plume_width - plume_error;
    % Add standard error
    fill([source_distance;source_distance(end:-1:1)],...
        [std_error_top; std_error_down(end:-1:1)], colors(i,:),...
        'Edgecolor', 'none', 'FaceAlpha', 0.25,'HandleVisibility','off');
    % plot mean
    plot(source_distance, plume_width,'Color',colors(i,:),...        
        'DisplayName',strjoin(strsplit(treatments{i},'_'),'-'));        
end

legend show;
% Save figure
saveas(gcf, fullfile(plots_folder,...
    sprintf('%s.fig','all_plumes')));
saveas(gcf, fullfile(plots_folder,...
    sprintf('%s.png','all_plumes')));

