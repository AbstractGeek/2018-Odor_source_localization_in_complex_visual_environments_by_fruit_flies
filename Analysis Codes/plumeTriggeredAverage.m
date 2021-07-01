% A simple script to obtain and plot plume triggered averages of FlySpeed,
% CastSpeed and SurgeSpeed

%% Initialize
% Clear and load the analyzed mat file
clear,clc;
close all;

% Find speed change dataset
matlist = dir('PlumeEncounterDataset*.mat');
if size(matlist,1)>1
    filename = uigetfile('PlumeEncounterDataset*.mat');
    load(filename);
else
    load(matlist(1).name);
end
treatments = fieldnames(assay);

%% Defaults
parameterNames = {'FlySpeed','CastSpeed','SurgeSpeed', 'AngVel',...
    'TurnAng', 'OrientationAngle'};
bias_list = [-102, -51, 0, 51, 102];

% Add orientation data
assay = append_orientation_angle(assay);


for n=1:length(bias_list)    
    bias = bias_list(n);
    
    % Initialize
    timeWindow = 0.25;   % Seconds (set after literature scan)
    odor_plume = assay.(treatments{1}).odor_plume;
    sam_freq = assay.(treatments{1}).sam_freq;
    
    halfWindow = timeWindow * sam_freq;
    windowSize = halfWindow * 2 + 1;
    trt_size = 30 * 100;  % Max treatment size * maximum entries
    all_size = trt_size * 20; % Multiplied by number of treatments
    % Output structures
    first_entry = struct();     % plume triggered data
    first_entry_names = struct();
    entries = struct();
    exits = struct();

    
    % For each parameter, run through all the treatments
    for i=1:length(parameterNames)
        % Initialize all vectors
        first_entry.(parameterNames{i}).all = NaN(all_size,windowSize);
        entries.(parameterNames{i}).all = NaN(all_size,windowSize);
        exits.(parameterNames{i}).all = NaN(all_size,windowSize);
        % Initialize indices
        all_first_entry = 1;
        all_ind_entry = 1;
        all_ind_exit = 1;
        
        for j=1:length(treatments)
            % Initialize treatment vectors
            first_entry.(parameterNames{i}).(treatments{j}) = NaN(trt_size,windowSize);
            first_entry_names.(parameterNames{i}).(treatments{j}) = {};
            entries.(parameterNames{i}).(treatments{j}) = NaN(trt_size,windowSize);
            exits.(parameterNames{i}).(treatments{j}) = NaN(trt_size,windowSize);
            % Initialize indices
            curr_ind_entry = 1;
            curr_ind_exit = 1;
            
            % Append plume triggered data into the vectors
            for k=1:size(assay.(treatments{j}).trajData,1)
                dataLength = size(assay.(treatments{j}).trajData{k,1}.expdata.pts,1);
                contacts = round(...
                    assay.(treatments{j}).trajData{k,1}.globalpara.OdorAxisDist,1)...
                    <= odor_plume;
                entries_exits = diff(contacts);
                entry_ind = find(entries_exits==1) + 1;
                exit_ind = find(entries_exits==-1) + 1;
                
                % Add biases
                entry_ind = entry_ind + bias;
                exit_ind = exit_ind + bias;
                              
                if ~isempty(entry_ind)
                    
                    % Process first_entry
                    first_entry_ind = entry_ind(1);
                    % Check (and ignore) the first entry if it is too near to
                    % the start or the end of the trajectory.
                    if (first_entry_ind - halfWindow > 1) && ...
                            (first_entry_ind + halfWindow < dataLength-1)
                        % Entry point is well within the trajectory.
                        % If they are velocities, assume that the subsequent
                        % point is the trigger. That way the trigger indices
                        % remain the same.
                        first_entry.(parameterNames{i}).(treatments{j})(k,:) = ...
                            assay.(treatments{j}).trajData{k,1}.globalpara. ...
                            (parameterNames{i})(first_entry_ind-halfWindow:first_entry_ind+halfWindow);
                        
                        first_entry_names.(parameterNames{i}).(treatments{j}) = ...
                            [first_entry_names.(parameterNames{i}).(treatments{j});...
                            {assay.(treatments{j}).trajData{k,1}.name}];
                    end
                    
                    % Process all entries
                    for l=1:length(entry_ind)
                        if (entry_ind(l) - halfWindow > 1) && ...
                                (entry_ind(l) + halfWindow < dataLength-1)
                            
                            entries.(parameterNames{i}).(treatments{j})(curr_ind_entry,:) = ...
                                assay.(treatments{j}).trajData{k,1}.globalpara. ...
                                (parameterNames{i})(entry_ind(l)-halfWindow:entry_ind(l)+halfWindow);
                            
                            % Increment curr_ind_entry
                            curr_ind_entry = curr_ind_entry + 1;
                        end
                    end
                    
                end
                
                if ~isempty(exit_ind)
                    
                    % Process all exits
                    for l=1:length(exit_ind)
                        if (exit_ind(l) - halfWindow > 1) && ...
                                (exit_ind(l) + halfWindow < dataLength-1)
                            
                            exits.(parameterNames{i}).(treatments{j})(curr_ind_exit,:) = ...
                                assay.(treatments{j}).trajData{k,1}.globalpara. ...
                                (parameterNames{i})(exit_ind(l)-halfWindow:exit_ind(l)+halfWindow);
                            
                            % Increment curr_ind_exit
                            curr_ind_exit = curr_ind_exit + 1;
                        end
                    end
                    
                end
                
                
            end
            
            % Treatment complete
            % Clean up (remove NaNs)
            first_entry.(parameterNames{i}).(treatments{j})(...
                isnan(first_entry.(parameterNames{i}).(treatments{j})(:,1)),:) = [];
            entries.(parameterNames{i}).(treatments{j})(...
                isnan(entries.(parameterNames{i}).(treatments{j})(:,1)),:) = [];
            exits.(parameterNames{i}).(treatments{j})(...
                isnan(exits.(parameterNames{i}).(treatments{j})(:,1)),:) = [];
            
            % Save the entries to all
            % first_entry
            curr_first_entry = ...
                size(first_entry.(parameterNames{i}).(treatments{j}),1)+1;
            first_entry.(parameterNames{i}). ...
                all(all_first_entry:all_first_entry+curr_first_entry-2,:) = ...
                first_entry.(parameterNames{i}).(treatments{j});
            all_first_entry = all_first_entry + curr_first_entry-1;
            % entries
            entries.(parameterNames{i}). ...
                all(all_ind_entry:all_ind_entry+curr_ind_entry-2,:) = ...
                entries.(parameterNames{i}).(treatments{j});
            all_ind_entry = all_ind_entry + curr_ind_entry-1;
            % exits
            exits.(parameterNames{i}). ...
                all(all_ind_exit:all_ind_exit+curr_ind_exit-2,:) = ...
                exits.(parameterNames{i}).(treatments{j});
            all_ind_exit = all_ind_exit + curr_ind_exit-1;
            
            
        end
        
        % Clean up alls (remove NaNs)
        first_entry.(parameterNames{i}).all(...
            isnan(first_entry.(parameterNames{i}).all(:,1)),:) = [];
        entries.(parameterNames{i}).all(...
            isnan(entries.(parameterNames{i}).all(:,1)),:) = [];
        exits.(parameterNames{i}).all(...
            isnan(exits.(parameterNames{i}).all(:,1)),:) = [];
        
    end
    
    %  Plume triggered data collection done. Save it into a nice structure
    ptd.first_entry = first_entry;
    ptd.first_entry_names = first_entry_names;
    ptd.entries = entries;
    ptd.exits = exits;
    ptd.treatments = treatments;
    ptd.odor_plume = odor_plume;
    ptd.timeWindow = timeWindow;
    ptd.sam_freq = sam_freq;
    ptd.halfWindow = halfWindow;
    ptd.windowSize = windowSize;
    ptd.encounterInd = (windowSize+1)/2;
    ptd.bias = bias;
    
    
    save(sprintf('PlumeTriggeredData (ptd)[window=%g](pw=%g)[bias=%d].mat',...
        ptd.timeWindow,odor_plume,bias),'ptd','odor_plume');
    clearvars -except ptd parameterNames treatments assay odor_plume bias bias_list;
    
    
    % Plot PlumeTriggeredData
    line_color = [0.8, 0.8, 0.8];
    mean_color = [0, 0.45, 0.74];
    std_error_color = [0, 0.45, 0.74];
    
    % Initialize-figure 1
    % h1=figure('Units','Pixels','Position',[60 60 1260 960]);
    h1=figure();
    time = (-ptd.timeWindow:1/ptd.sam_freq:ptd.timeWindow)';
    outfolder = fullfile(pwd, sprintf('pta (pw=%g)[bias=%d]',odor_plume, bias));
    if ~isdir(outfolder)
        mkdir(outfolder)
    end
    
    for i=1:length(parameterNames)
        
        if isempty(ptd.first_entry.(parameterNames{i}).all)
            continue;
        end
        
        % first entry
        % Plot individual trajectories
        plot(time,ptd.first_entry.(parameterNames{i}).all',...
            'Color', line_color);
        hold on;
        
        % Obtain mean and std
        if strcmp(parameterNames{i}, 'OrientationAngle')
            % Calculate circular mean and ste for orientation angle
            mean_temp = circ_mean(ptd.first_entry.(parameterNames{i}).all * pi/180, [], 1) * 180/pi;
            std_error_temp = circ_std(ptd.first_entry.(parameterNames{i}).all * pi/180, [], 1) * 180/pi ./sqrt(size(ptd.first_entry.(parameterNames{i}).all,1));            
        else
            % Calculate mean and standard error
            mean_temp = mean(ptd.first_entry.(parameterNames{i}).all,1);
            std_error_temp = std(ptd.first_entry.(parameterNames{i}).all,1)./sqrt(size(ptd.first_entry.(parameterNames{i}).all,1));
        end
        
        % Obtain borders
        std_error_top = mean_temp + std_error_temp;
        std_error_down = mean_temp - std_error_temp;
        
        % Plot standard error and the mean
        fill([time;time(end:-1:1)],...
            [std_error_top, std_error_down(end:-1:1)]', std_error_color,...
            'Edgecolor', 'none', 'FaceAlpha', 0.25);
        plot(time,mean_temp', 'Color', mean_color, 'LineWidth', 2.0);
        % Clean up axis
        axis tight;
        v = round(axis);
        plot([0,0],[v(3),v(4)],'r');
        xlabel('time (s)');
        ylabel(sprintf('%s (cm/s)', parameterNames{i}));
        title(sprintf('All %s first entry', parameterNames{i}));
        %TufteStyle(gca);
        % Save the figure
        saveas(h1,fullfile(outfolder,...
            sprintf('%s_first_entry (window=%g).fig',parameterNames{i},ptd.timeWindow)));
        set(h1,'PaperPositionMode','auto');
        print(fullfile(outfolder,...
            sprintf('%s_first_entry (window=%g).png',parameterNames{i},ptd.timeWindow)),'-dpng');
        clf(h1);
        
        % Entries
        % Plot individual trajectories
        plot(time,ptd.entries.(parameterNames{i}).all',...
            'Color', line_color);
        hold on;
        % Calculate mean and standard error
        mean_temp = mean(ptd.entries.(parameterNames{i}).all,1);
        std_error_temp = std(ptd.entries.(parameterNames{i}).all,1)./sqrt(size(ptd.entries.(parameterNames{i}).all,1));
        std_error_top = mean_temp + std_error_temp;
        std_error_down = mean_temp - std_error_temp;
        % Plot standard error and the mean
        fill([time;time(end:-1:1)],...
            [std_error_top, std_error_down(end:-1:1)]', std_error_color,...
            'Edgecolor', 'none', 'FaceAlpha', 0.25);
        plot(time,mean_temp', 'Color', mean_color, 'LineWidth', 2.0);
        % Clean up axis
        axis tight;
        v = round(axis);
        plot([0,0],[v(3),v(4)],'r');
        xlabel('time (s)');
        ylabel(sprintf('%s (cm/s)', parameterNames{i}));
        title(sprintf('All %s entries', parameterNames{i}));
        %TufteStyle(gca);
        saveas(h1,fullfile(outfolder,...
            sprintf('%s_entries (window=%g).fig',parameterNames{i},ptd.timeWindow)));
        set(h1,'PaperPositionMode','auto');
        print(fullfile(outfolder,...
            sprintf('%s_entries (window=%g).png',parameterNames{i},ptd.timeWindow)),'-dpng');
        clf(h1);
        
        % Exits
        % Plot individual trajectories
        plot(time,ptd.exits.(parameterNames{i}).all',...
            'Color', line_color);
        hold on;
        % Calculate mean and standard error
        mean_temp = mean(ptd.exits.(parameterNames{i}).all,1);
        std_error_temp = std(ptd.exits.(parameterNames{i}).all,1)./sqrt(size(ptd.exits.(parameterNames{i}).all,1));
        std_error_top = mean_temp + std_error_temp;
        std_error_down = mean_temp - std_error_temp;
        % Plot standard error and the mean
        fill([time;time(end:-1:1)],...
            [std_error_top, std_error_down(end:-1:1)]', std_error_color,...
            'Edgecolor', 'none', 'FaceAlpha', 0.25);
        plot(time,mean_temp', 'Color', mean_color, 'LineWidth', 2.0);
        % Clean up axis
        axis tight;
        v = round(axis);
        plot([0,0],[v(3),v(4)],'r');
        xlabel('time (s)');
        ylabel(sprintf('%s (cm/s)', parameterNames{i}));
        title(sprintf('All %s exits', parameterNames{i}));
        % TufteStyle(gca);
        saveas(h1,fullfile(outfolder,...
            sprintf('%s_exits (window=%g).fig',parameterNames{i},ptd.timeWindow)));
        set(h1,'PaperPositionMode','auto');
        print(fullfile(outfolder,...
            sprintf('%s_exits (window=%g).png',parameterNames{i},ptd.timeWindow)),'-dpng');
        clf(h1);
    end
    close(h1);
    
end
