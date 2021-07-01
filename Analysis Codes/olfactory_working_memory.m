%% A simple script to quantify olfactory working memory.
% Dinesh Natesan
% 5th Oct, 2017

% Variables
treatment_name = cell(10000,1);
trial_name = cell(10000,1);
landing_type = cell(10000,1);
memory_time = nan(10000,1);
curr_ind = 1;

plume_width = 0.8;  %cm (radius)
sam = 100;

% Iterate through treatments
treatments = fieldnames(assay);

for i=2:12
    
    % Iterate through trajectories
    for j=1:length(assay.(treatments{i}).trajData)
        
        % Process individual trajectories
        currTraj = assay.(treatments{i}).trajData{j,1};
        [~,landingPoint] = min(sqrt((...
            currTraj.expdata.cs(:,1)-currTraj.expdata.pts(end,1)).^2+...
            (currTraj.expdata.cs(:,2)-currTraj.expdata.pts(end,2)).^2+...
            (currTraj.expdata.cs(:,3)-currTraj.expdata.pts(end,3)).^2));
        
        % Find the last odor encounter
        if (currTraj.globalpara.OdorAxisDist(end)>plume_width)
            
            % Check if correct or incorrect landing
            if landingPoint==1
                landing_type{curr_ind} = "correct";
            else
                landing_type{curr_ind} = "incorrect";              
                
                % Find odor encounter
                odor_encounter = diff(currTraj.globalpara.OdorAxisDist<=0.8);
                ind = find(odor_encounter==1,1,'last');
                if isempty(ind)
                    ind = 1;
                end
                
                fprintf("memory time: %f\n",(length(currTraj.globalpara.OdorAxisDist) - (ind + 1)) / sam);
                memory_time(curr_ind) = (length(currTraj.globalpara.OdorAxisDist) - (ind + 1)) / sam;
                
                % Add other stuff
                treatment_name{curr_ind} = treatments{i};
                trial_name{curr_ind} = currTraj.name;
                
                % Increment current index
                curr_ind = curr_ind + 1;
            end
        end
        
    end
    
end

% Table for last encounter
last_encounter = table(treatment_name(1:curr_ind-1),...
    trial_name(1:curr_ind-1),landing_type(1:curr_ind-1),...
    memory_time(1:curr_ind-1),'VariableNames',...
    {'treatment','trial','landing','time'});


% Plot histogram
bin_size = 2 * iqr(last_encounter.time) * length(last_encounter.time)^(-1/3);
figure, histogram(last_encounter.time,0:bin_size:6);
xlabel('Time between last odor encounter and landing (s)');
ylabel('Occurence (number)');

% writetable(last_encounter, 'last_encounter.csv');
