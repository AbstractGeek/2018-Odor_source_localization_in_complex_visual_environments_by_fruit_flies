% A simple script to plot velocity at odor encounters, and calculate the
% mean change before and after the velocity encounter.

clear, close all;
timeWindow = 0.1;
odor_plume = 0.5;
load(sprintf('PlumeTriggeredData (ptd)[window=%g](pw=%g).mat',...
    timeWindow,odor_plume));
% Pull out necessary variables
outfolder_overall = 'z_all';
treatments = ptd.treatments;
types = fieldnames(ptd);
types = types(cell2mat(cellfun(@(x) isstruct(ptd.(x)),...
    types,'UniformOutput',false)));
parameters = fieldnames(ptd.(types{1}));

time = (-ptd.timeWindow:1/ptd.sam_freq:ptd.timeWindow)';
h1=figure('Units','Pixels','Position',[10 10 1260 960]);

% get outfolder
outfolder = fullfile(pwd,outfolder_overall,...
    sprintf('plume_triggered_speed_change(window=%g)[pw=%g]',...
    ptd.timeWindow,odor_plume));
if ~isdir(outfolder)
    mkdir(outfolder);
end

for i=1:length(treatments)
    
    for j=1:length(types)
        % Get pdf folder
        pdffolder = fullfile(pwd, treatments{i}(2:end),...
            sprintf('%s_plume_triggered_speed_change(window=%g)[pw=%g]',...
            treatments{i}(2:end),ptd.timeWindow,odor_plume),...
            types{j});
        
        for k=1:length(parameters)
            % Savefolder
            savefolder = fullfile(pdffolder, parameters{k});
            % Create folders if necessary
            if ~isdir(savefolder)
                mkdir(savefolder);
            end
            % Get number of individual trajectoies
            entry_num = size(ptd.(types{j}).(parameters{k}).(treatments{i}),1);
            
            % Process each entry individually
            for l = 1:entry_num
                current_entry = (ptd.(types{j}).(parameters{k}).(treatments{i})(l,:))';
                % calculate mean and std
                mean_before_encounter = mean(current_entry(1:ptd.encounterInd));
                std_before_encounter = std(current_entry(1:ptd.encounterInd));
                mean_after_encounter = mean(current_entry(ptd.encounterInd:end));
                std_after_encounter = std(current_entry(ptd.encounterInd:end));
                
                % plot the data                                
                subplot(1,2,1), hold on;
                plot(time,current_entry,'k');
                v = round(axis);
                plot([0,0],[v(3),v(4)],'Color',[0.0118 0.2627 0.8745]);
                xlabel('time (s)');
                ylabel(sprintf('%s (cm/s)', parameters{k}));
                
                subplot(1,2,2);
                errorbar([-1,1], [mean_before_encounter, mean_after_encounter],...
                    [std_before_encounter,std_after_encounter]);
                axis([-2 2 v(3) v(4)]);
                xlabel('Position w.r.t encounter');
                ylabel(sprintf('%s (cm/s)', parameters{k}));
                
                % save as fig and as a pdf
                figtitle(sprintf('%s %s %s - %d', treatments{i},...
                    parameters{k}, strjoin(strsplit(types{j},'_'),'-'), l),...
                    'fontweight','bold');
                saveas(h1,fullfile(savefolder,sprintf('%d.fig',l)));
                orient landscape;
                print(gcf, '-dpdf', '-bestfit', fullfile(savefolder,...
                    sprintf('%d.pdf',l)));
                clf(h1);
                
            end
            % Done - continue onto the next parameter
        end
        % Generate pdfs and copy them into the outfolder
        commandstring = sprintf('~/Scripts/combinePDFs.py "%s"',pdffolder);
        system(commandstring);
        % copy the pdfs into the overall folder
        for k=1:length(parameters)
            copyfolder = fullfile(outfolder, sprintf('%s_%s',...
                parameters{k}, types{j}));
            currfile = fullfile(pdffolder,sprintf('%s.pdf',parameters{k}));
            
            if ~isdir(copyfolder)
                mkdir(copyfolder);
            end
            
            if exist(currfile,'file')==2
                copyfile(currfile,fullfile(copyfolder,...
                    sprintf('%s.pdf',treatments{i})));
            end
            
        end
        % Done. Continue onto the next type
    end
    % Done. Continue onto the next treatment
end
close(h1);