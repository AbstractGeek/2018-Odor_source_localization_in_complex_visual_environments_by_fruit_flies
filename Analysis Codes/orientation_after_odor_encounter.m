% A script file to quantify orientation change after odor encounter
% Dinesh Natesan
% 5th Oct, 2017

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
bias_list = [-102, 0, 102];

for n=1:length(bias_list)
    
    % Current bias
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
    entries = struct();
    exits = struct();
    
    for j=1:length(treatments)
        % Initialize treatment vectors
        first_entry.(treatments{j}) = NaN(trt_size,windowSize);
        entries.(treatments{j}) = NaN(trt_size,windowSize);
        exits.(treatments{j}) = NaN(trt_size,windowSize);
        
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
                    first_entry.(treatments{j})(k,:) = ...
                        assay.(treatments{j}).trajData{k,1}.globalpara. ...
                        (first_entry_ind-halfWindow:first_entry_ind+halfWindow);
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
