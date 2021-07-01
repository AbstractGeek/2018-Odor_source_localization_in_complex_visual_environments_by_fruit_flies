% A simple script to plot speed distributions before and after odor contact.

% Initialize
clear, close all;
stats_tablename = {'Difference', 'Comparision', 'SE', 'qobs', 'qexp','H','Conclusion'};

% Find speed change dataset
matlist = dir('PlumeTriggeredData*.mat');

h1=figure();

for n = 1:size(matlist,1)    
    % Load data
    load(matlist(n).name);
    % Pull out necessary variables
    timeWindow = ptd.timeWindow;
    odor_plume = ptd.odor_plume;
    bias = ptd.bias;
    treatments = ptd.treatments;
    types = fieldnames(ptd);
    types = types(cell2mat(cellfun(@(x) isstruct(ptd.(x)),...
        types,'UniformOutput',false)));
    parameters = fieldnames(ptd.(types{1}));
    
    time = (-ptd.timeWindow:1/ptd.sam_freq:ptd.timeWindow)';    
    distributionStats = struct;
    
    % Create an rmd file for the plots
    rmdfile = sprintf('plume_triggered_distributions(window=%g)[pw=%g][bias=%d].rmd',...
        ptd.timeWindow,odor_plume, bias);
    sfid = fopen(fullfile(pwd,rmdfile),'w');   % Disregard contents if any
    fprintf(sfid,'---\n');
    fprintf(sfid,'title: "%s"\n', ...
        sprintf('Plume Triggered Distribution Plots (window=%g) [pw=%g] [bias=%d]',...
        ptd.timeWindow,odor_plume, bias));
    fprintf(sfid,'author: "%s"\n','Dinesh');
    fprintf(sfid,'date: "%s"\n',datestr(datetime('now')));
    fprintf(sfid,'---\n\n');
    
    %% Histogram compare
    % Probability distribution compare
    % get outfolder
    outfolder = fullfile(pwd,...
        sprintf('plume_triggered_distributions(window=%g)[pw=%g][bias=%d]',...
        ptd.timeWindow,odor_plume, bias));
    if ~isdir(outfolder)
        mkdir(outfolder);
    end
    
    for i=1:length(treatments)
        
        for j=1:length(types)
            
            fprintf(sfid,'# Type - %s\n\n', types{j});
            
            for k=1:length(parameters)
                fprintf(sfid,'## %s\n\n', parameters{k});
                % Do seperate things based on the parameter
                switch parameters{k}
                    case 'FlySpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Fly Speed';
                    case 'CastSpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Cast Speed';
                    case 'SurgeSpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Surge Speed';
                    case 'TurnAng'
                        hist_edges = 0:1:150;
                        xlabel_str = 'Turn Angle';
                    case 'AngVel'
                        hist_edges = 0:100:17000;
                        xlabel_str = 'Angular Velocity';
                    case 'Time'
                        hist_edges = 0:0.5:18;
                        xlabel_str = 'Trajectory duration';
                    case 'FirstContactTime'
                        hist_edges = 0:0.5:12;
                        xlabel_str = 'Time of first contact w.r.t landing';
                    case 'FirstContactXPos'
                        hist_edges = 0:0.5:12;
                        xlabel_str = 'X position of first contact w.r.t landing';
                    case 'FirstContactSPos'
                        hist_edges = 0:0.5:15;
                        xlabel_str = 'Distance of first contact w.r.t landing';
                    case 'PlumeTime'
                        hist_edges = 0:0.1:2;
                        xlabel_str = 'Time spent in Plume';
                    otherwise
                        warning('Unexpected parameter type. Skipping.')
                        continue;
                end
                
                % Get before contact parameters
                beforeContact = ptd.(types{j}).(parameters{k}).(treatments{i})(:,1:ptd.encounterInd);
                meanBeforeContact = mean(beforeContact,2);
                beforeContact = beforeContact(:);
                
                % Get after contact parameters
                afterContact = ptd.(types{j}).(parameters{k}).(treatments{i})(:,ptd.encounterInd:end);
                meanAfterContact = mean(afterContact,2);
                afterContact = afterContact(:);
                                
                % Plot occurence histograms
                histogram(beforeContact,...
                    hist_edges,'DisplayName','Before Odor Contact');
                hold on;
                histogram(afterContact,...
                    hist_edges,'DisplayName','After Odor Contact');
                
                legend show;
                title(sprintf('%s-%s-%s',treatments{i},parameters{k},...
                    strjoin(strsplit(types{j},'_'),'')));
                xlabel(xlabel_str);
                ylabel('Occurences');
                %TufteStyle(gca);
                saveas(h1,fullfile(outfolder,...
                    sprintf('%s_%s_occurences.fig',parameters{k},types{j})));
                set(h1,'PaperPositionMode','auto');
                print(fullfile(outfolder,...
                    sprintf('%s_%s_occurences.png',parameters{k},types{j})),'-dpng');
                clf(h1);                
                fprintf(sfid,'### Occurence Histograms (%s %s)\n',...
                    parameters{k}, types{j});
                fprintf(sfid,'![](%s)\n\n',...
                    fullfile(outfolder,...
                    sprintf('%s_%s_occurences.png',parameters{k},types{j})));                
                % Done with occurence histograms
                
                % Plot probability histograms
                histogram(beforeContact,...
                    hist_edges,'Normalization','probability',...
                    'DisplayName','Before Odor Contact');
                hold on;
                histogram(afterContact,...
                    hist_edges,'Normalization','probability',...
                    'DisplayName','After Odor Contact');
                
                legend show;
                title(sprintf('%s-%s-%s',treatments{i},parameters{k},...
                    strjoin(strsplit(types{j},'_'),'')));
                xlabel(xlabel_str);
                ylabel('Probability');
                %TufteStyle(gca);
                saveas(h1,fullfile(outfolder,...
                    sprintf('%s_%s_probability.fig',parameters{k},types{j})));
                set(h1,'PaperPositionMode','auto');
                print(fullfile(outfolder,...
                    sprintf('%s_%s_probability.png',parameters{k},types{j})),'-dpng');
                clf(h1);                
                % Done with probability histograms
                
                % box plots
                boxplot([meanBeforeContact, meanAfterContact], ...
                    'Labels', {'Before Contact', 'After Contact'});
                title(sprintf('%s-%s-%s',treatments{i},parameters{k},...
                    strjoin(strsplit(types{j},'_'),'')));                
                ylabel(xlabel_str);
                %TufteStyle(gca);
                saveas(h1,fullfile(outfolder,...
                    sprintf('%s_%s_boxplot.fig',parameters{k},types{j})));
                set(h1,'PaperPositionMode','auto');
                print(fullfile(outfolder,...
                    sprintf('%s_%s_boxplot.png',parameters{k},types{j})),'-dpng');
                clf(h1);
                
                % categorical scatter plots
                X = [meanBeforeContact; meanAfterContact];
                Groups = [ones(size(meanBeforeContact));...
                    2.*ones(size(meanAfterContact))];
                CategoricalScatterplot(X, Groups,...
                    'Labels', {'Before Contact', 'After Contact'});
                title(sprintf('%s-%s-%s',treatments{i},parameters{k},...
                    strjoin(strsplit(types{j},'_'),'')));                
                ylabel(xlabel_str);
                %TufteStyle(gca);
                saveas(h1,fullfile(outfolder,...
                    sprintf('%s_%s_CategoricalScatterplot.fig',parameters{k},types{j})));
                set(h1,'PaperPositionMode','auto');
                print(fullfile(outfolder,...
                    sprintf('%s_%s_CategoricalScatterplot.png',parameters{k},types{j})),'-dpng');
                clf(h1);
                fprintf(sfid,'### Categorical Scatter Plots [mean] (%s %s)\n',...
                    parameters{k}, types{j});
                fprintf(sfid,'![](%s)\n\n',...
                    fullfile(outfolder,...
                    sprintf('%s_%s_CategoricalScatterplot.png',parameters{k},types{j})));        
                % Done with categorical scatter plots
                
                % Perform stats test on them
                [h, P, stats] = PerformStats(X, Groups);
                
                % Save into a structure                
                distributionStats.(types{j}).(parameters{k}).h = h;
                distributionStats.(types{j}).(parameters{k}).P = P;
                distributionStats.(types{j}).(parameters{k}).stats = stats;
                % distributionStats.(types{j}).(parameters{k}).stats_all = stats_all;                    
                distributionStats.(types{j}).(parameters{k}).X = X;
                distributionStats.(types{j}).(parameters{k}).Group = Groups;
                
                % Save it into the rmd file
                fprintf(sfid,'#### Kruskal Wallis test : H=%d, p=%g\n',...
                    h, P);
                stats_table = cell2table( ...
                    stats.table(2:end,:),...
                    'VariableNames', stats_tablename);
                fprintf(sfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                
                % continue onto the next parameter
                
            end
            
        end
        % Done. Continue onto the next type
    end
    % Done. Continue onto the next treatment
    
    distributionStats.bias = ptd.bias;    
    % save the structure
    save(sprintf('PlumeTriggeredDistributionStats[window=%g](pw=%g)[bias=%d].mat',...
        ptd.timeWindow,odor_plume,bias),'distributionStats');
    % Close and convert rmd file
    fclose(sfid);
    system(strcat('Rscript -e "rmarkdown::render(''',fullfile(pwd,rmdfile),''')"'));
    
    
    %% Probability distribution compare
    % get outfolder
    outfolder = fullfile(pwd,...
        sprintf('plume_triggered_pdf(window=%g)[pw=%g][bias=%d]',...
        ptd.timeWindow,odor_plume,bias));
    if ~isdir(outfolder)
        mkdir(outfolder);
    end
    
    for i=1:length(treatments)
        
        for j=1:length(types)            
            
            for k=1:length(parameters)
                
                % Do seperate things based on the parameter
                switch parameters{k}
                    case 'FlySpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Fly Speed';
                    case 'CastSpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Cast Speed';
                    case 'SurgeSpeed'
                        hist_edges = 0:1:70;
                        xlabel_str = 'Surge Speed';
                    case 'TurnAng'
                        hist_edges = 0:1:150;
                        xlabel_str = 'Turn Angle';
                    case 'AngVel'
                        hist_edges = 0:100:17000;
                        xlabel_str = 'Angular Velocity';
                    case 'Time'
                        hist_edges = 0:0.5:18;
                        xlabel_str = 'Trajectory duration';
                    case 'FirstContactTime'
                        hist_edges = 0:0.5:12;
                        xlabel_str = 'Time of first contact w.r.t landing';
                    case 'FirstContactXPos'
                        hist_edges = 0:0.5:12;
                        xlabel_str = 'X position of first contact w.r.t landing';
                    case 'FirstContactSPos'
                        hist_edges = 0:0.5:15;
                        xlabel_str = 'Distance of first contact w.r.t landing';
                    case 'PlumeTime'
                        hist_edges = 0:0.1:2;
                        xlabel_str = 'Time spent in Plume';
                    otherwise
                        warning('Unexpected parameter type. Skipping.')
                        continue;
                end
                
                % Get before contact parameters
                beforeContact = ptd.(types{j}).(parameters{k}).(treatments{i})(:,1:ptd.encounterInd);
                beforeContact = beforeContact(:);
                beforeContact(beforeContact == 0) = [];
                if isempty(beforeContact)
                    continue;
                end
                % Get after contact parameters
                afterContact = ptd.(types{j}).(parameters{k}).(treatments{i})(:,ptd.encounterInd:end);
                afterContact = afterContact(:);
                afterContact(afterContact == 0) = [];
                if isempty(afterContact)
                    continue;
                end
                
                % Fit a log normal distribution - beforeContact
                lnpd_bC = fitdist(beforeContact,'LogNormal');
                ln_fit_bC = lognpdf(hist_edges,lnpd_bC.mu,lnpd_bC.sigma);
                % Fit a log normal distribution - afterContact
                lnpd_aC = fitdist(afterContact,'LogNormal');
                ln_fit_aC = lognpdf(hist_edges,lnpd_aC.mu,lnpd_aC.sigma);
                
                % Plot histograms
                histogram(beforeContact,...
                    hist_edges,'Normalization','probability',...
                    'DisplayName','Before Odor Contact');
                hold on;
                histogram(afterContact,...
                    hist_edges,'Normalization','probability',...
                    'DisplayName','After Odor Contact');
                
                % Plot fits
                plot(hist_edges,ln_fit_bC,'Color',rgb('blue'),...
                    'DisplayName', sprintf('%s (mu-%0.2f, sigma-%0.2f)',...
                    'BeforeContact', lnpd_bC.mu, lnpd_bC.sigma),...
                    'LineWidth',2.0);
                plot(hist_edges,ln_fit_aC,'Color',rgb('red'),...
                    'DisplayName', sprintf('%s (mu-%0.2f, sigma-%0.2f)',...
                    'After Contact', lnpd_aC.mu, lnpd_aC.sigma),...
                    'LineWidth',2.0);
                
                legend show;                
                title(sprintf('%s-%s-%s',treatments{i},parameters{k},...
                    strjoin(strsplit(types{j},'_'),'')));
                xlabel(xlabel_str);
                ylabel('Probability');
                saveas(h1,fullfile(outfolder,...
                    sprintf('%s_%s.fig',parameters{k},types{j})));
                set(h1,'PaperPositionMode','auto');
                print(fullfile(outfolder,...
                    sprintf('%s_%s.png',parameters{k},types{j})),'-dpng');                
                clf(h1);
                
                % Done - continue onto the next parameter
            end
            
            % Done. Continue onto the next type
        end
        % Done. Continue onto the next treatment
    end
        
end

% Done - close figure
close(h1);