function [result] = analyzeSortedMatFiles(varargin)
% function [] = analyzeSortedMatFiles()
%
%
%

%% Initialization
% Defaults
stats_tablename = {'Difference', 'Comparision', 'SE', 'qobs', 'qexp','H','Conclusion'};
params_nodisp = {'CurveArea', 'HorzArea', 'VertArea', 'PlumeCast'};
bin_nodisp = {'SourceDistBin', 'OdorAxisBin'};
% Get the mat file that has all the Odor-Tracking data
if isempty(varargin)
    % Get input Directory
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the whole odor-tracking data mat file');
    % if Cancel button is pressed
    if (FilterIndex == 0)
        disp('Action cancelled, quiting');
        return;
    end
    
else
    inputFile = GetFullPath(varargin{1});
    selection = varargin{2};
    
    [PathName,FileName,ext] = fileparts(inputFile);
    FileName = strcat(FileName, ext);
end

% Enter and load mat files
cd (PathName);
load (FileName);

% Get the treatments to compare and analyze
treatmentNames = fieldnames(assay);
if isempty(varargin)
    [selection,ok] = listdlg('ListString',treatmentNames,'SelectionMode','multiple','Name','Select Treatment','PromptString','Select Treatments to Compare and Analyze');
    % selection = [] or if Cancel button is pressed
    if (ok == 0)
        disp('No Selection / Action cancelled, quiting');
        return;
    end
end

%Selected Treatments
selectedTreatments = treatmentNames(selection);
%Save basic details for the binning from the selection
treatmentStats.binValues = assay.(selectedTreatments{1}).binValues;

%% Start processing selected data
% Get sizes of treatments & obtain a list of selected treatments
treatmentSize = ones(length(selection),1);
% Start for loop to make a group and find size
treatmentSize(1) = size(assay.(selectedTreatments{1}).OverallSum.Distance,1);
treatmentStats.Group = ones(treatmentSize(1),1);
names = regexp(selectedTreatments{1},'t','split');
savefolder = names{2};
for i=2:length(selection)
    treatmentSize(i) = size(assay.(selectedTreatments{i}).OverallSum.Distance,1);
    treatmentStats.Group = [treatmentStats.Group; ones(treatmentSize(i),1).*i];
    names = regexp(selectedTreatments{i},'t','split');
    savefolder = strcat(savefolder,'-',names{2});
end

% Create the dest-folder based on obtained fields
destpath = fullfile(PathName,strcat(savefolder,'-Stats Results'));
insigpath = fullfile(PathName,strcat(savefolder,'-Stats Results'),'The Insignificants');
% destpath = fullfile(PathName);
% insigpath = fullfile(PathName,'The Insignificants');

if ~isdir(destpath)
    mkdir(destpath);
end

% Setup the stats structure
treatmentStats.selectedTreatments = selectedTreatments;
treatmentStats.treatmentSize = treatmentSize;

h1 = figure;

%% Write into a markdown file
% Open plot file and enter significant plots
rmdfile = 'significant_plots.rmd';
sfid = fopen(fullfile(destpath,rmdfile),'w');   % Disregard contents if any
fprintf(sfid,'---\n');
fprintf(sfid,'title: "%s"\n', ...
    strjoin(cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '),...
    selectedTreatments, 'UniformOutput', false),' - '));
fprintf(sfid,'author: "%s"\n','Dinesh');
fprintf(sfid,'date: "%s"\n',datestr(datetime('now')));
fprintf(sfid,'---\n\n');

% Open plot file and enter insignificant plots
rmdfile = 'insignificant_plots.rmd';
insfid = fopen(fullfile(destpath,rmdfile),'w');   % Disregard contents if any
fprintf(insfid,'---\n');
fprintf(insfid,'title: "%s"\n', ...
    strjoin(cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '),...
    selectedTreatments, 'UniformOutput', false),' - '));
fprintf(insfid,'author: "%s"\n','Dinesh');
fprintf(insfid,'date: "%s"\n',datestr(datetime('now')));
fprintf(insfid,'---\n\n');

% Open plot file and enter significant plots
rmdfile = 'other_parameter_plots.rmd';
opfid = fopen(fullfile(destpath,rmdfile),'w');   % Disregard contents if any
fprintf(opfid,'---\n');
fprintf(opfid,'title: "%s"\n', ...
    strjoin(cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '),...
    selectedTreatments, 'UniformOutput', false),' - '));
fprintf(opfid,'author: "%s"\n','Dinesh');
fprintf(opfid,'date: "%s"\n',datestr(datetime('now')));
fprintf(opfid,'---\n\n');


%% Process Overalls
% OverallSums
sumNames = fieldnames(assay.(selectedTreatments{1}).OverallSum);
savepath = fullfile(destpath,'OverallSum');
if ~isdir(savepath)
    mkdir(savepath);
end
insigSavePath = fullfile(insigpath,'OverallSum');
if ~isdir(insigSavePath)
    mkdir(insigSavePath);
end

fprintf(sfid,'# Overall Sums\n');
fprintf(insfid,'# Overall Sums\n');
fprintf(opfid,'# Overall Sums\n');

for i=1:size(sumNames)
    % Get X
    treatmentStats.OverallSum.(char(sumNames(i))).X = assay.(selectedTreatments{1}).OverallSum.(char(sumNames(i)));
    for j=2:length(selection)
        treatmentStats.OverallSum.(char(sumNames(i))).X = [treatmentStats.OverallSum.(char(sumNames(i))).X; assay.(selectedTreatments{j}).OverallSum.(char(sumNames(i)))];
    end
    % Perform Stats
    [treatmentStats.OverallSum.(char(sumNames(i))).h,treatmentStats.OverallSum.(char(sumNames(i))).P,...
        treatmentStats.OverallSum.(char(sumNames(i))).stats,treatmentStats.OverallSum.(char(sumNames(i))).stats_all] = ...
        PerformStats(treatmentStats.OverallSum.(char(sumNames(i))).X,treatmentStats.Group);
    
    % Save if H0 is rejected
    if strcmp(treatmentStats.OverallSum.(char(sumNames(i))).stats.table{2,7},'Reject H0')
        fprintf('OverallSum - %s \n',char(sumNames(i)));
        boxplot(treatmentStats.OverallSum.(char(sumNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(sumNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(savepath,strcat(char(sumNames(i)),'_BoxPlot')),'fig');
        clf(h1);
        % Add categorical scatter plots
        categoricalscatterplot(treatmentStats.OverallSum.(char(sumNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(sumNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(savepath,strcat(char(sumNames(i)),'_CategoricalScatterPlot')),'fig');
        set(h1,'PaperPositionMode','auto');
        print(fullfile(savepath,strcat(char(sumNames(i)),'_CategoricalScatterPlot')),'-dpng');
        % add into the markdown file
        if any(ismember(params_nodisp, char(sumNames(i))))
            fprintf(opfid, '## %s (%s-Significant)\n',char(sumNames(i)), 'Overall Sum');
            fprintf(opfid,'![](%s)\n\n',...
                fullfile(savepath,strcat(char(sumNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallSum.(char(sumNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        else
            fprintf(sfid, '## %s (%s)\n',char(sumNames(i)), 'Overall Sum');
            fprintf(sfid,'![](%s)\n\n',...
                fullfile(savepath,strcat(char(sumNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallSum.(char(sumNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(sfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        end
        clf(h1);
        
    else
        boxplot(treatmentStats.OverallSum.(char(sumNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(sumNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(insigSavePath,strcat(char(sumNames(i)),'_BoxPlot')),'fig');
        clf(h1);
        % Add categorical scatter plots
        categoricalscatterplot(treatmentStats.OverallSum.(char(sumNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(sumNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(insigSavePath,strcat(char(sumNames(i)),'_CategoricalScatterPlot')),'fig');
        set(h1,'PaperPositionMode','auto');
        print(fullfile(insigSavePath,strcat(char(sumNames(i)),'_CategoricalScatterPlot')),'-dpng');
        % add into the markdown file
        if any(ismember(params_nodisp, char(sumNames(i))))
            fprintf(opfid, '## %s (%s Insignificant)\n',char(sumNames(i)), 'Overall Sum');
            fprintf(opfid,'![](%s)\n\n',...
                fullfile(insigSavePath,strcat(char(sumNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallSum.(char(sumNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        else
            fprintf(insfid, '## %s (%s)\n',char(sumNames(i)), 'Overall Sum');
            fprintf(insfid,'![](%s)\n\n',...
                fullfile(insigSavePath,strcat(char(sumNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallSum.(char(sumNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(insfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        end
        
        clf(h1);
    end
end

% OverallAverages
avgNames = fieldnames(assay.(selectedTreatments{1}).OverallAverage);
savepath = fullfile(destpath,'OverallAverage');
if ~isdir(savepath)
    mkdir(savepath);
end
insigSavePath = fullfile(insigpath,'OverallAverage');
if ~isdir(insigSavePath)
    mkdir(insigSavePath);
end

fprintf(sfid,'# Overall Averages\n');
fprintf(insfid,'# Overall Averages\n');
fprintf(opfid,'# Overall Averages\n');

for i=1:size(avgNames)
    % Get X
    treatmentStats.OverallAverage.(char(avgNames(i))).X = assay.(selectedTreatments{1}).OverallAverage.(char(avgNames(i)));
    for j=2:length(selection)
        treatmentStats.OverallAverage.(char(avgNames(i))).X = [treatmentStats.OverallAverage.(char(avgNames(i))).X; assay.(selectedTreatments{j}).OverallAverage.(char(avgNames(i)))];
    end
    % Perform Stats
    [treatmentStats.OverallAverage.(char(avgNames(i))).h,treatmentStats.OverallAverage.(char(avgNames(i))).P,...
        treatmentStats.OverallAverage.(char(avgNames(i))).stats,treatmentStats.OverallAverage.(char(avgNames(i))).stats_all] = ...
        PerformStats(treatmentStats.OverallAverage.(char(avgNames(i))).X,treatmentStats.Group);
    
    % Save if H0 is rejected
    if strcmp(treatmentStats.OverallAverage.(char(avgNames(i))).stats.table{2,7},'Reject H0')
        fprintf('OverallAverage - %s \n',char(avgNames(i)));
        boxplot(treatmentStats.OverallAverage.(char(avgNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(avgNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(savepath,strcat(char(avgNames(i)),'_BoxPlot')),'fig');
        clf(h1);
        % Add categorical scatter plot
        categoricalscatterplot(treatmentStats.OverallAverage.(char(avgNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(avgNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(savepath,strcat(char(avgNames(i)),'_CategoricalScatterPlot')),'fig');
        set(h1,'PaperPositionMode','auto');
        print(fullfile(savepath,strcat(char(avgNames(i)),'_CategoricalScatterPlot')),'-dpng');
        % add into the markdown file
        if any(ismember(params_nodisp, char(avgNames(i))))
            fprintf(opfid, '## %s (%s Significant)\n',char(avgNames(i)), 'Overall Average');
            fprintf(opfid,'![](%s)\n\n',...
                fullfile(savepath,strcat(char(avgNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallAverage.(char(avgNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
            
        else
            fprintf(sfid, '## %s (%s)\n',char(avgNames(i)), 'Overall Average');
            fprintf(sfid,'![](%s)\n\n',...
                fullfile(savepath,strcat(char(avgNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallAverage.(char(avgNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(sfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        end
        
        clf(h1);
    else
        boxplot(treatmentStats.OverallAverage.(char(avgNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(avgNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(insigSavePath,strcat(char(avgNames(i)),'_BoxPlot')),'fig');
        clf(h1);
        % Add categorical scatter plot
        categoricalscatterplot(treatmentStats.OverallAverage.(char(avgNames(i))).X,treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
        ylabel(char(avgNames(i)));
        TufteStyle(gca);
        saveas(h1,fullfile(insigSavePath,strcat(char(avgNames(i)),'_CategoricalScatterPlot')),'fig');
        set(h1,'PaperPositionMode','auto');
        print(fullfile(insigSavePath,strcat(char(avgNames(i)),'_CategoricalScatterPlot')),'-dpng');
        % add into the markdown file
        if any(ismember(params_nodisp, char(avgNames(i))))
            fprintf(opfid, '## %s (%s Insignificant)\n',char(avgNames(i)), 'Overall Average');
            fprintf(opfid,'![](%s)\n\n',...
                fullfile(insigSavePath,strcat(char(avgNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallAverage.(char(avgNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        else
            
            fprintf(insfid, '## %s (%s)\n',char(avgNames(i)), 'Overall Average');
            fprintf(insfid,'![](%s)\n\n',...
                fullfile(insigSavePath,strcat(char(avgNames(i)),'_CategoricalScatterPlot.png')));
            stats_table = cell2table( ...
                treatmentStats.OverallAverage.(char(avgNames(i))).stats.table(2:end,:),...
                'VariableNames', stats_tablename);
            fprintf(insfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
        end
        
        clf(h1);
    end
end

%%%% Overalls Processed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Inter-treatment Comparision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sortedBinNames = fieldnames(assay.(selectedTreatments{1}).sortedBins);      % SourceDistBin,OdorAxisBin or RectBin
variableNames =  fieldnames(assay.(selectedTreatments{1}).sortedBins.(sortedBinNames{1}));

fprintf(sfid,'# Inter-treatment Comparision\n');
fprintf(insfid,'# Inter-treatment Comparision\n');
fprintf(opfid,'# Inter-treatment Comparision\n');

for i=1:size(sortedBinNames)
    binSize = assay.(selectedTreatments{1}).binNumbers.(strcat(char(sortedBinNames(i)),'Number'));
    
    for m=1:size(variableNames)
        
        savepath = fullfile(destpath,'Inter-Treatment Comparision',char(sortedBinNames(i)),char(variableNames(m)));
        
        if ~isdir(savepath)
            mkdir(savepath);
        end
        
        insigSavePath = fullfile(insigpath,'Inter-Treatment Comparision',char(sortedBinNames(i)),char(variableNames(m)));
        if ~isdir(insigSavePath)
            mkdir(insigSavePath);
        end
        
        for j=1:binSize
            % Get X
            treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X = assay.(selectedTreatments{1}).sortedBins.(char(sortedBinNames(i))).(char(variableNames(m)))(j,:)';
            
            for k=2:length(selection)
                treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X = ...
                    [treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X; assay.(selectedTreatments{k}).sortedBins.(char(sortedBinNames(i))).(char(variableNames(m)))(j,:)'];
            end
            
            [treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).h,...
                treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).P,...
                treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats,...
                treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats_all] = ...
                PerformStats(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X,treatmentStats.Group);
            
            % Save if H0 is rejected
            if strcmp(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats.table{2,7},'Reject H0')
                fprintf('InterTreatment Plots-%s \n',strcat(char(sortedBinNames(i)),'-',char(variableNames(m)),'-DistanceBin-',int2str(j)));
                boxplot(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X,...
                    treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
                ylabel(sprintf('%s-%s(%d)',char(variableNames(m)),char(sortedBinNames(i)),j));
                TufteStyle(gca);
                saveas(h1,fullfile(savepath,strcat('DistanceBin_',int2str(j),'_BoxPlot')),'fig');
                clf(h1);
                % Add categorical scatter plot
                categoricalscatterplot(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X,...
                    treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
                ylabel(sprintf('%s-%s(%d)',char(variableNames(m)),char(sortedBinNames(i)),j));
                TufteStyle(gca);
                saveas(h1,fullfile(savepath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot')),'fig');
                set(h1,'PaperPositionMode','auto');
                print(fullfile(savepath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot')),'-dpng');
                % add into the markdown file
                if any(ismember(params_nodisp, char(variableNames(m)))) || any(ismember(bin_nodisp, char(sortedBinNames(i))))
                    fprintf(opfid, '## %s (%s-%d-Significant)\n',char(variableNames(m)),char(sortedBinNames(i)),j);
                    fprintf(opfid,'![](%s)\n\n',...
                        fullfile(savepath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                else
                    fprintf(sfid, '## %s (%s-%d)\n',char(variableNames(m)),char(sortedBinNames(i)),j);
                    fprintf(sfid,'![](%s)\n\n',...
                        fullfile(savepath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(sfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                end
                
                clf(h1);
            else
                boxplot(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X,...
                    treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
                ylabel(sprintf('%s-%s(%d)',char(variableNames(m)),char(sortedBinNames(i)),j));
                TufteStyle(gca);
                saveas(h1,fullfile(insigSavePath,strcat('DistanceBin_',int2str(j),'_BoxPlot')),'fig');
                clf(h1);
                % Add categorical scatter plot
                categoricalscatterplot(treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).X,...
                    treatmentStats.Group,'labels',cellfun(@(x) strjoin(strsplit(x(2:end),'_'), ' '), selectedTreatments, 'UniformOutput', false));
                ylabel(sprintf('%s-%s(%d)',char(variableNames(m)),char(sortedBinNames(i)),j));
                TufteStyle(gca);
                saveas(h1,fullfile(insigSavePath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot')),'fig');
                set(h1,'PaperPositionMode','auto');
                print(fullfile(insigSavePath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot')),'-dpng');
                % add into the markdown file
                if any(ismember(params_nodisp, char(variableNames(m)))) || any(ismember(bin_nodisp, char(sortedBinNames(i))))
                    fprintf(opfid, '## %s (%s-%d-Insignificant)\n',char(variableNames(m)),char(sortedBinNames(i)),j);
                    fprintf(opfid,'![](%s)\n\n',...
                        fullfile(insigSavePath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                else
                    fprintf(insfid, '## %s (%s-%d)\n',char(variableNames(m)),char(sortedBinNames(i)),j);
                    fprintf(insfid,'![](%s)\n\n',...
                        fullfile(insigSavePath,strcat('DistanceBin_',int2str(j),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.interTreatment.sortedBins.(char(sortedBinNames(i))).(char(variableNames(m))).(strcat('DistanceBin',int2str(j))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(insfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                end
                
                clf(h1);
            end
            
        end
        
        
    end
    
end

%%%% Intra-treatment Comparision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(sfid,'# Intra-treatment Comparision\n');
fprintf(insfid,'# Intra-treatment Comparision\n');
fprintf(opfid,'# Intra-treatment Comparision\n');

for i=1:size(sortedBinNames)
    
    for j=1:size(treatmentStats.selectedTreatments,1)
        [binSize,treatmentSize] = size(assay.(selectedTreatments{j}).sortedBins.(char(sortedBinNames(i))).(char(variableNames(1))));
        Group = ones(binSize,treatmentSize);
        for n=1:binSize
            Group(n,:) = Group(n,:).*n;
        end
        Group = Group';
        Group = Group(:);
        
        savepath = fullfile(destpath,'Intra-Treatment Comparision',(selectedTreatments{j}),char(sortedBinNames(i)));
        if ~isdir(savepath)
            mkdir(savepath);
        end
        
        insigSavePath = fullfile(insigpath,'Intra-Treatment Comparision',(selectedTreatments{j}),char(sortedBinNames(i)));
        if ~isdir(insigSavePath)
            mkdir(insigSavePath);
        end
        
        for m=1:size(variableNames)
            
            temp = assay.(selectedTreatments{j}).sortedBins.(char(sortedBinNames(i))).(char(variableNames(m)))';
            treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X = temp(:);
            
            [treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).h,...
                treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).P,...
                treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats,...
                treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats_all]=...
                PerformStats(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X,Group);
            
            if strcmp(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats.table{2,7},'Reject H0')
                fprintf('IntraTreatment Plots-%s \n',strcat(char(sortedBinNames(i)),'-',(selectedTreatments{j}),'-',char(variableNames(m))));
                boxplot(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X,Group);
                ylabel(sprintf('%s-%s',char(variableNames(m)),selectedTreatments{j}(2:end)));
                xlabel(char(sortedBinNames(i)));
                TufteStyle(gca);
                saveas(h1,fullfile(savepath,strcat(char(variableNames(m)),'_BoxPlot')),'fig');
                clf(h1);
                % Add categorical scatter plot
                categoricalscatterplot(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X,Group);
                ylabel(sprintf('%s-%s',char(variableNames(m)),selectedTreatments{j}(2:end)));
                xlabel(char(sortedBinNames(i)));
                TufteStyle(gca);
                saveas(h1,fullfile(savepath,strcat(char(variableNames(m)),'_CategoricalScatterPlot')),'fig');
                set(h1,'PaperPositionMode','auto');
                print(fullfile(savepath,strcat(char(variableNames(m)),'_CategoricalScatterPlot')),'-dpng');
                % add into the markdown file
                if any(ismember(params_nodisp, char(variableNames(m)))) || any(ismember(bin_nodisp, char(sortedBinNames(i))))
                    fprintf(opfid, '## %s (%s-%s-Significant)\n',char(variableNames(m)), selectedTreatments{j},char(sortedBinNames(i)),j);
                    fprintf(opfid,'![](%s)\n\n',...
                        fullfile(savepath,strcat(char(variableNames(m)),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                else
                    fprintf(sfid, '## %s (%s-%s)\n',char(variableNames(m)), selectedTreatments{j},char(sortedBinNames(i)),j);
                    fprintf(sfid,'![](%s)\n\n',...
                        fullfile(savepath,strcat(char(variableNames(m)),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(sfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                end
                
                clf(h1);
            else
                boxplot(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X,Group);
                ylabel(sprintf('%s-%s',char(variableNames(m)),selectedTreatments{j}(2:end)));
                xlabel(char(sortedBinNames(i)));
                TufteStyle(gca);
                saveas(h1,fullfile(insigSavePath,strcat(char(variableNames(m)),'_BoxPlot')),'fig');
                clf(h1)
                % Add categorical scatter plot
                categoricalscatterplot(treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).X,Group);
                ylabel(sprintf('%s-%s',char(variableNames(m)),selectedTreatments{j}(2:end)));
                xlabel(char(sortedBinNames(i)));
                TufteStyle(gca);
                saveas(h1,fullfile(insigSavePath,strcat(char(variableNames(m)),'_CategoricalScatterPlot')),'fig');
                set(h1,'PaperPositionMode','auto');
                print(fullfile(insigSavePath,strcat(char(variableNames(m)),'_CategoricalScatterPlot')),'-dpng');
                % add into the markdown file
                if any(ismember(params_nodisp, char(variableNames(m)))) || any(ismember(bin_nodisp, char(sortedBinNames(i))))
                    fprintf(opfid, '## %s (%s-%s-Insignificant)\n',char(variableNames(m)), selectedTreatments{j},char(sortedBinNames(i)),j);
                    fprintf(opfid,'![](%s)\n\n',...
                        fullfile(insigSavePath,strcat(char(variableNames(m)),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(opfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                else
                    
                    fprintf(insfid, '## %s (%s-%s)\n',char(variableNames(m)), selectedTreatments{j},char(sortedBinNames(i)),j);
                    fprintf(insfid,'![](%s)\n\n',...
                        fullfile(insigSavePath,strcat(char(variableNames(m)),'_CategoricalScatterPlot.png')));
                    stats_table = cell2table( ...
                        treatmentStats.intraTreatment.sortedBins.(selectedTreatments{j}).(char(sortedBinNames(i))).(char(variableNames(m))).stats.table(2:end,:),...
                        'VariableNames', stats_tablename);
                    fprintf(insfid, '%s\n\n', print_table(stats_table,'colSepStr', '  '));
                end
                
                clf(h1);
            end
        end
        
    end
end


save([destpath,'/','comparisionStatistics'],'treatmentStats');
close(h1);
fclose(sfid);
fclose(insfid);
fclose(opfid);
% Generate html files
system(strcat('Rscript -e "rmarkdown::render(''',fullfile(destpath,'significant_plots.rmd'),''')"'));
system(strcat('Rscript -e "rmarkdown::render(''',fullfile(destpath,'insignificant_plots.rmd'),''')"'));
system(strcat('Rscript -e "rmarkdown::render(''',fullfile(destpath,'other_parameter_plots.rmd'),''')"'));

result = 1;

end

function [h,P,stats,stats_all] = PerformStats(X,Group)
% function [h,stats] = PerformStats(X,Group)
%
%
%

alpha = 0.05;

% Split vector into subvectors and remove NaN's if there are any.
data = cell(max(Group),1);
for i = 1:size(X,1)
    if (isnan(X(i,1))~= 1)
        data{Group(i,1),1} = [data{Group(i,1),1};X(i,1)];
    end
end

% Rank the obtained data.
[rankedData,tiedNums] = rankData(data);
% Analysis of variance
[h,P] = kWallis(rankedData,tiedNums);
fprintf('\nHypothesis = %d; p = %f \n',h,P);
% Post-hoc analysis - Tukey test
[stats] = nonParaTukey(rankedData,tiedNums,alpha);
stats.X = X;
stats.Group = Group;
[stats_all] = nonParaTukeyAllPairwise(rankedData,tiedNums,alpha);
stats_all.Group = Group;
stats_all.X = X;

end

function [h,P] = kWallis(rankedData,tiedNums)
% [h,P] = kWallis(rankedData,tiedNums)
%
%

alpha = 0.05;       % Setting alpha to 0.05. Change if necessary
% % start
% data = {varargin{1:end-1}}';
% tiedNums = varargin{end};

numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
tempsum = 0;
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
    tempsum = tempsum + (R(i)*R(i)/N(i));
end
% Total N
sumN = sum(N);
% Kruskal-Wallis Statistic, H
H = (12/(sumN*(sumN+1)))*tempsum - 3*(sumN+1);
% Correction Factor
t = sum(tiedNums.^3 - tiedNums);
C = 1-(t/(sumN.^3-sumN));
% Corrected H
Hcor = H/C;
v = numOfInputs-1;
P = 1-chi2cdf(Hcor,v);

if (P<alpha)
    h=1;            % Reject the Null hypothesis H0 (i.e. distributions are not the same)
else
    h=0;            % Failure to reject Null hypothesis H0 (i.e. distributions are same)
end

end

function [stats] = nonParaTukey(rankedData,tiedNums,alpha)
% [stats] = nonparaTukey(rankedData,tiedNums,alpha)
%
%
%
load 'Nonparametric Tukey Test Critical Values.mat';
numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
end
% Total N
sumN = sum(N);
t = sum(tiedNums.^3 - tiedNums);
meanRanks = R./N;
stats.meanRanksOrdered = meanRanks;
% Arrange mean ranks and preserve group number
stats.meanRanks = NaN(length(meanRanks),2);
[stats.meanRanks(:,1),stats.meanRanks(:,2)] = sort(meanRanks,'ascend');
meanRanks = stats.meanRanks;
% Get Q value and get the number of comparisions
Qexp=nonParaTukeyCriticalValues(numOfInputs,nonParaTukeyCriticalValues(1,:)==alpha);          %#ok<NODEF>
numOfComp = (numOfInputs*(numOfInputs-1)/2);         % Number of comparisions

% Define stats - output variable, and make it a proper table
stats.table = cell(numOfComp+1,7);
stats.table{1,1} = 'Comparision (B vs. A)';
stats.table{1,2} = 'Difference (mean(RA) - mean(RB))';
stats.table{1,3} = 'SE';
stats.table{1,4} = 'Q';
stats.table{1,5} = 'Q (obtained from the table))';
stats.table{1,6} = 'Null Hypothesis - H0';
stats.table{1,7} = 'Conclusion';

% Loop variables
currentB = numOfInputs;     % Ranks
currentA = 1;

flag = 0;
doNotCompare = {};

for i=1:numOfComp
    % Get necessary parameters for B
    RB = meanRanks(currentB,1);
    Bindex = meanRanks(currentB,2);
    nB = N(Bindex);
    % Get necessary parameters for A
    RA = meanRanks(currentA,1);
    Aindex = meanRanks(currentA,2);
    nA = N(Aindex);
    
    for j=1:length(doNotCompare)
        tmp = doNotCompare{j};
        if (sum(tmp==currentB)&&sum(tmp==currentA))
            flag=1;
        end
    end
    
    if (flag==1)
        stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
        stats.table{i+1,7} = 'Do Not Test';
        flag = 0;
    else
        % Get SE,Q and h
        SE = sqrt(((sumN*(sumN+1)/12) - (t/(12*(sumN-1))))*(1/nA+1/nB));
        Q = (RB - RA)/SE;
        h = (Q>Qexp);
        
        if h
            conclusion = 'Reject H0';
        else
            conclusion = 'Accept H0';
            tmp = currentA:currentB;
            if(length(tmp)>2 && ~isempty(doNotCompare))
                dNClength = length( doNotCompare );
                doNotCompare{ dNClength + 1 } = tmp; %#ok<AGROW>
                %                 doNotCompare = {doNotCompare{1:end};tmp};
            elseif (length(tmp)>2 && isempty(doNotCompare))
                doNotCompare = {tmp};
            end
        end
        
        stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
        stats.table{i+1,2} = RB-RA;
        stats.table{i+1,3} = SE;
        stats.table{i+1,4} = Q;
        stats.table{i+1,5} = Qexp;
        stats.table{i+1,6} = h;
        stats.table{i+1,7} = conclusion;
        
    end
    
    % Propogate the loop
    currentA = currentA+1;
    if (currentA == currentB && currentB==1)
        % loop done
        disp('Something wrong, your loop is not ending properly! Check!!');
    elseif (currentA == currentB)
        currentB = currentB-1;
        currentA = 1;
    end
    
end

end

function [stats] = nonParaTukeyAllPairwise(rankedData,tiedNums,alpha)
% [stats] = nonparaTukeyAllPairwise(rankedData,tiedNums,alpha)
%
%
%
load 'Nonparametric Tukey Test Critical Values.mat';
numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
tempsum = 0; %#ok<NASGU>
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
end
% Total N
sumN = sum(N);
t = sum(tiedNums.^3 - tiedNums);
meanRanks = R./N;
meanRanksOrdered = meanRanks; %#ok<NASGU>
% Arrange mean ranks and preserve group number
stats.meanRanks = NaN(length(meanRanks),2);
[stats.meanRanks(:,1),stats.meanRanks(:,2)] = sort(meanRanks,'ascend');
meanRanks = stats.meanRanks;
% Get Q value and get the number of comparisions
Qexp=nonParaTukeyCriticalValues(numOfInputs,nonParaTukeyCriticalValues(1,:)==alpha);          %#ok<NODEF>
numOfComp = (numOfInputs*(numOfInputs-1)/2);         % Number of comparisions

% Define stats - output variable, and make it a proper table
stats.table = cell(numOfComp+1,7);
stats.table{1,1} = 'Comparision (B vs. A)';
stats.table{1,2} = 'Difference (mean(RA) - mean(RB))';
stats.table{1,3} = 'SE';
stats.table{1,4} = 'Q';
stats.table{1,5} = 'Q (obtained from the table))';
stats.table{1,6} = 'Null Hypothesis - H0';
stats.table{1,7} = 'Conclusion';

% Loop variables
currentB = numOfInputs;     % Ranks
currentA = 1;

for i=1:numOfComp
    % Get necessary parameters for B
    RB = meanRanks(currentB,1);
    Bindex = meanRanks(currentB,2);
    nB = N(Bindex);
    % Get necessary parameters for A
    RA = meanRanks(currentA,1);
    Aindex = meanRanks(currentA,2);
    nA = N(Aindex);
    
    % Get SE,Q and h
    SE = sqrt(((sumN*(sumN+1)/12) - (t/(12*(sumN-1))))*(1/nA+1/nB));
    Q = (RB - RA)/SE;
    h = (Q>Qexp);
    
    if h
        conclusion = 'Reject H0';
    else
        conclusion = 'Accept H0';
    end
    
    stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
    stats.table{i+1,2} = RB-RA;
    stats.table{i+1,3} = SE;
    stats.table{i+1,4} = Q;
    stats.table{i+1,5} = Qexp;
    stats.table{i+1,6} = h;
    stats.table{i+1,7} = conclusion;
    
    % Propogate the loop
    currentA = currentA+1;
    if (currentA == currentB && currentB==1)
        % loop done
        disp('Something wrong, your loop is not ending properly! Check!!');
        break;
    elseif (currentA == currentB)
        currentB = currentB-1;
        currentA = 1;
    end
    
end

end

function [rankedData,tiedNums] = rankData(data)
% [rankedData,tiedNums] = rankData(data)
% This is a precursor function for the non parametric tests - Mann-Whitney and Kruskal Wallis test
% It outputs grouped data with its ranks
%
% Input - Cell array with observations
% Output - 1)Cell array with ranks
%        - 2)tiedNums (number of ties per rank)

% data = varargin';
numOfInputs = length(data);
R = tiedrank(cell2mat(data));
rankedData = cell(numOfInputs,1);
% Break the tied ranks into groups and return the cell
temp = R;
for i = 1:numOfInputs
    cellend = size(data{i,1},1);
    %     varargout{i} = temp(1:cellend);
    rankedData{i} = temp(1:cellend);
    if i~=numOfInputs
        temp = temp(cellend+1:end);
    end
end

tiedNums = hist(R,(1:max(R))')';
% varargout{numOfInputs+1} = tiedNums;

end

function [] = categoricalscatterplot(X, Group, varargin)
% function [] = categoricalscatterplot()
%
%
%
% Dinesh Natesan
% Last Modified: 16th Aug 2016

%% Parse Inputs
p = inputParser;
% Main arguments
addRequired(p, 'X', @(X) ismatrix(X));
addRequired(p, 'Group', @(X) ismatrix(X) || iscellstr(X));
addOptional(p, 'Color', 'k', @(X) all(size(X) == [length(unique(Group)), 3]) ||...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addOptional(p, 'Labels', false, @(X) iscellstr(X));

% Scatter parameters
addParameter(p, 'binWidth', false, @(x) isnumeric(x));
addParameter(p, 'binWidthRatio', 0.05, @(x) isnumeric(x));
addParameter(p, 'spreadWidth', 0.6, @(x) isnumeric(x));
addParameter(p, 'boxWidth', 0.6, @(x) isnumeric(x));
% Plotting styles
addParameter(p, 'Marker', 'o', @(X) (ischar(X) && length(X)==1));
addParameter(p, 'MarkerSize', 25, @(x) isnumeric(x));
addParameter(p, 'FillMarker', true, @(x) islogical(x));
addParameter(p, 'BoxColor', [0.31, 0.31, 0.31], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxEdgeColor', 'none', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && (length(X)==1 || length(X)==4)));
addParameter(p, 'MedianColor', 'r', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'WhiskerColor', [0 0 0], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxAlpha', 0.50, @(x) isnumeric(x));
addParameter(p, 'BoxLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'MedianLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'WhiskerLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'BoxLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'MedianLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'WhiskerLineWidth', 1.0, @(x) isnumeric(x));
% Parse inputs and unpack structure
parse(p, X, Group, varargin{:});
% v2struct(p.Results);
parsed = p.Results;

% %% Remove NaNs if they exist
% nan_ind = isnan(X);
% X = X(~nan_ind);
% Group = Group(~nan_ind);
%  
%% Convert the groups into a cell array
groups = unique(Group);
data = cell(length(groups), 2);
new_data = cell(length(groups), 2);
Xlen = 0;
Xmax = -999;
Xmin = 999;

for i = 1:length(groups)
    data{i,1} = X(Group == groups(i));
    data{i,2} = Group(Group == groups(i));
    if (Xlen < length(data{i,1}))
        Xlen = length(data{i,1});
    end
    
    if (Xmin > floor(min(data{i,1})))
        Xmin = floor(min(data{i,1}));
    end
    
    if (Xmax < ceil(max(data{i,1})))
        Xmax = ceil(max(data{i,1}));
    end
    
end

% Get binWidth ratio is only ratio is provided
if (islogical(parsed.binWidth) && ~parsed.binWidth)
    binWidth = parsed.binWidthRatio * round(Xmax - Xmin);
else
    binWidth = parsed.binWidth;
end

%% Discretize points in a group
for i = 1:length(groups)
    Xtemp = data{i,1};
    Ytemp = data{i,2};
    if (binWidth > 0)
        [counts,~,bins] = histcounts(Xtemp, 'BinWidth', binWidth);
    else
        [counts,~,bins] = histcounts(Xtemp, 1);
    end
    inds = find(counts~=0);
    counts = counts(inds);
    
    for j=1:length(inds)
        width = parsed.spreadWidth * (1-exp(-0.1 * (counts(j)-1)));
        xpoints = linspace(-width/2, width/2, counts(j)) + i;
        Ytemp(bins==inds(j)) = xpoints;
    end
    
    new_data{i,1} = Xtemp;
    new_data{i,2} = Ytemp;
    
end

%% Plot the data beautifully
boxWidth = parsed.boxWidth;
hold on;

for i = 1:length(groups)
    imp_quantiles = quantile(new_data{i,1}, [0.25, 0.5, 0.75]);
    IQR = imp_quantiles(3) - imp_quantiles(1);
    whisker = [imp_quantiles(1) - 1.5 * IQR, ...
        imp_quantiles(3) + 1.5 * IQR];
    
    % Draw box
    patch([i-boxWidth/2, i-boxWidth/2, i+boxWidth/2, i+boxWidth/2]',...
        [imp_quantiles(3), imp_quantiles(1), imp_quantiles(1), imp_quantiles(3)]',...
        parsed.BoxColor, 'FaceAlpha', parsed.BoxAlpha,...
        'EdgeColor', parsed.BoxEdgeColor, ...
        'LineStyle', parsed.BoxLineStyle, 'LineWidth', parsed.BoxLineWidth);
    
    % Draw points
    if parsed.FillMarker
        if (ischar(parsed.Color)||size(parsed.Color,1))
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color);
        else
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color(i,:));
        end
    else
        if (ischar(parsed.Color)||size(parsed.Color,1))
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color);
        else
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color(i,:));
        end
    end
    
    % Draw median
    plot([i-boxWidth/2, i+boxWidth/2], ...
        [imp_quantiles(2), imp_quantiles(2)], ...
        'LineStyle', parsed.MedianLineStyle, ...
        'Color', parsed.MedianColor,...
        'LineWidth', parsed.MedianLineWidth);
    
    % Draw Q + 1.5 IQR
    temp = sortrows([whisker(2) - data{i,1}, (1:length(data{i,1}))'], 1);
    closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
    
    if isempty(closest_point)
        continue;
    end
    plot([i-boxWidth/5, i+boxWidth/5], ...
        [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);    
    
    % Draw top whiskers
    plot([i, i], [new_data{i,1}(closest_point), imp_quantiles(3)],...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);    
    
    % Draw Q - 1.5 IQR
    temp = sortrows([data{i,1} - whisker(1), (1:length(data{i,1}))'], 1);
    closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
    plot([i-boxWidth/5, i+boxWidth/5], ...
        [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);
    
    % Draw bottom whiskers
    plot([i, i], [new_data{i,1}(closest_point), imp_quantiles(1)],...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);
    
    
end

ax = gca;
ax.XTick = 1:length(groups);
if (islogical(parsed.Labels) && ~parsed.Labels)
    ax.XTickLabel = groups;
else
    ax.XTickLabel = parsed.Labels;
end

end

%% Print table
function varargout = print_table(dataCell, varargin)
%% PRINT_TABLE Print data in a table format (text or latex)
%
% Syntax:
%     PRINT_TABLE(dataTable)
%     PRINT_TABLE(dataCell)
%
%     PRINT_TABLE(__, dataDescCellstr)
%     PRINT_TABLE(__, dataDescCellstr, headerColumnCellstr, headerRowCellstr)
%
%     PRINT_TABLE(__, Name, Value)
%
%     tableStr = PRINT_TABLE(__)
%
% Input:
%     dataTable            - table data to print (see table)
%     dataCell             - cell with data to print (can be numeric matrix)
%
%     dataDescCellstr      - cell with sprintf syntax for elements in data
%                            Note, dataDesc is expanded to the complete table
%                            size if only a single element, row or column
%                            description is supplied.
%     headerColumnCellstr  - cell array with column header names
%     headerRowCellstr     - cell array with row header name
%
%    Note, if both headerRow/Column are supplied, one can be one element longer
%    than the dimension of the dataCell, the extra element (which should be the
%    first element in the array) is then positioned at the top left corner, i.e.
%
%       |-------|-------|-------| ... |-------|
%       | EXTRA | hCol1 | hCol2 | ... | hColN |  -> Header Columns
%       |-------|-------|-------| ... |-------|
%       | hRow1 | d(1,1)| d(1,2)| ... | d(1,N)|  \
%       |   :        :       :      :      :  |   > Data part of table
%       | hRowM | d(M,1)| d(M,2)| ... | d(M,N)|  /
%       |-------|-------|-------| ... |-------|
%          \/
%        Header
%         Rows
%
% Options, supplied as (..., Name, Value) pairs, overrides default values:
%     printHeaderCol = 1   - print header columns (if supplied)
%     printHeaderRow = 1   - print header rows (if supplied)
%     transposeTable = 0   - transpose table compared to input format
%     printMode = 'text'   - print mode, 'text' or 'latex'
%        colSepStr = '|'   - separation string between columns (if 'text')
%        rowSepStr = ''    - separation line character between rows (if 'text')
%        rowHSepStr = '-'  - separation line character between header and data
%        colHSepStr = ''   - extra separation string between col.header and data
%     textAlignment  = 'c' - text alignment in each column (alt. 'l' or 'r')
%     	Note, is possible to supply for each column as string, e.g. 'lcl...cr'.
%     numSpaceColPad = 1   - extra space padding in each column
%     spaceColPadAlign = 1 - use the extra space padding with the text alignment
%        Note, cosmetic change if we do not want the extra space padding to be
%        included in the aligned text, e.g. 'lText   ' -> ' lText  ', if false
%        and numSpaceColPad = 1 and textAlignment = 'l'.
%     printLatexFull = 1   - add tabular enviroment to latex table format
%     printBorder    = 0   - print simple border around the table (in text mode)
%       borderRowStr = '-' - border type string, should be single character
%
% Output:
%     Table printed in command window, or
%     tableStr    - string with output table, preferably printed using fprintf
%
% Comment:
%     Utility function for writing a table in either text mode or latex
%     mode. Generates a table with aligned columns by inserting spaces.
%     Can also print the transposed version of the supplied table data.
%     Provides a variety of options such as adding/removing extra space padding
%     or changing the separation characters.
%
%     Includes three additional utility functions:
%        repmat_as_needed  - repmat data into a specified size
%        rmexpzeroes       - removes unecessary zeroes from exponent string
%        cellstr2str       - concatenates cellstr with some char between parts
%
% Example usage:
%  print_table(1e2.*rand(5,3),{'%.3g'},{'a','b','c'},{'No.','1','2','3','4','5'})
%  No. |   a  |   b  |   c
% -----|------|------|------
%   1  |  45  | 28.5 | 27.5
%   2  | 20.6 | 67.3 | 71.7
%   3  |  90  | 66.4 | 28.3
%   4  | 76.3 | 12.3 | 89.6
%   5  | 88.2 | 40.7 | 82.7
%  print_table(1e2.*rand(5,3),{'%.3g'},{'a','b','c'},{'No.','1','2','3','4','5'},'printBorder',1)
% |-----|------|------|------|
% | No. |   a  |   b  |   c  |
% |-----|------|------|------|
% |  1  |  43  | 10.9 | 22.9 |
% |  2  | 69.4 |  39  | 83.4 |
% |  3  | 94.5 | 59.1 | 1.56 |
% |  4  | 78.4 | 45.9 | 86.4 |
% |  5  | 70.6 | 5.03 | 7.81 |
% |-----|------|------|------|
% print_table(1e2.*rand(5,3),{'%.3g'},{'a','b','c'},{'No.','1','2','3','4','5'},'printMode','latex')
% \begin{tabular}{|c|c|c|c|}\hline
%  No. &   a  &   b  &   c  \\ \hline
%   1  & 66.9 & 67.1 & 1.96 \\ \hline
%   2  &  50  &  60  & 43.5 \\ \hline
%   3  & 21.8 &  5.6 & 83.2 \\ \hline
%   4  & 57.2 & 5.63 & 61.7 \\ \hline
%   5  & 12.2 & 15.3 &  52  \\ \hline
% \end{tabular}
%
% See also fprintf, inputParser, table, table2cell
% repmat_as_needed, rmexpzeros, cellstr2str

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2014/10/06 14:00:00$
%   $Revision: 1.2$  $Date: 2014/10/07 10:00:00$
%     Changed name, added option input format, new options, some error control
%   $Revision: 1.3$  $Date: 2014/10/23 14:00:00$
%     Changed option input format to use inputParser object
%   $Revision: 1.4$  $Date: 2014/10/30 15:00:00$
%     Added textAlignment option for each column, and spaceColPadAlign option
%   $Revision: 1.5$  $Date: 2014/10/31 11:00:00$
%     Added support for table data type input

%% Set default input and parse input:
inpPar = inputParser;
addRequired(inpPar,'dataCell',...
    @(x) validateattributes(x,{'cell','numeric','table'},{'nonempty'}));
addOptional(inpPar,'dataDescCellstr',{'%g'},...
    @(x) validateattributes(x,{'cell'},{'nonempty'}));
% Due to bug/feature? in parse, we can not allow an optional parameter to be a
% string. Otherwise, we would change {'cell'} -> {'char','cell'}.
addOptional(inpPar,'headerColumnCellstr',[]);
addOptional(inpPar,'headerRowCellstr',[]);
addParameter(inpPar,'printMode','text',...
    @(x) any(validatestring(x,{'text','latex'})));
addParameter(inpPar,'printHeaderRow',true,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','scalar','binary'}))
addParameter(inpPar,'printHeaderCol',true,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','scalar','binary'}))
addParameter(inpPar,'transposeTable',false,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','scalar','binary'}))
addParameter(inpPar,'colSepStr','|',@ischar);
addParameter(inpPar,'rowSepStr','',@ischar);
addParameter(inpPar,'colHSepStr','',@ischar);
addParameter(inpPar,'rowHSepStr','-',@ischar);
addParameter(inpPar,'textAlignment','c',...
    @(x) all(x=='c' | x=='r' | x=='l'));
addParameter(inpPar,'spaceColPadAlign',true,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','scalar','binary'}))
addParameter(inpPar,'numSpaceColPad',1,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','integer','positive'}))
addParameter(inpPar,'printLatexFull',true,@islogical);
addParameter(inpPar,'printBorder',0,...
    @(x) validateattributes(x,{'numeric','logical'},...
    {'nonempty','scalar','binary'}))
addParameter(inpPar,'borderRowStr','-',@ischar);

% Parse input:
parse(inpPar, dataCell, varargin{:});
%% TODO-> bugged when not supplying optional input?!

% Collect some data from inpPar results:
dataDescCellstr = inpPar.Results.dataDescCellstr;
headerColumnCellstr = inpPar.Results.headerColumnCellstr;
headerRowCellstr = inpPar.Results.headerRowCellstr;

% Remake dataDescCellstr to cell if not cell
if ~iscell(dataDescCellstr)
    dataDescCellstr = {dataDescCellstr};
end

% Check if dataCell is actually in table format:
if istable(dataCell)
    dataTable   = dataCell;
    dataCell    = table2cell(dataTable);
    % If there is no input for headers, use variable names in table:
    headerRowCellstr     = dataTable.Properties.RowNames;
    headerColumnCellstr  = dataTable.Properties.VariableNames;
    % Generate dataDescCellstr if only specified as single element and not
    % matching the data types in the table:
    if ~( (size(dataDescCellstr,1) == size(dataCell,1)) || ...
            (size(dataDescCellstr,2) == size(dataCell,2)) )
        % Find out data types:
        dataIsNumLog   = cellfun(@(data) all(isnumeric(data)) | ...
            all(islogical(data)), dataCell);
        dataIsChar  = cellfun(@(data) ischar(data), dataCell);
        if ~all(dataIsNumLog(:))
            dataDescCellstr_init = dataDescCellstr;
            dataDescCellstr = cell(size(dataCell));
            dataDescCellstr(dataIsNumLog) = dataDescCellstr_init;
            dataDescCellstr(dataIsChar)   = {'%s'};
        end
    end
end

% Remake dataCell to cell if numeric:
if isnumeric(dataCell)
    dataCell = num2cell(dataCell);
end

colSepStr = inpPar.Results.colSepStr;
rowSepStr = inpPar.Results.rowSepStr;

%% Set the separation charachters to be used in the table if latex:
if strcmp(inpPar.Results.printMode,'latex')
    colSepStr = '&';
    newLineSepStr = '\\\\ \\hline \n';
    % Note, we need to repeat \ signs due to usage of sprintf.
else
    newLineSepStr = '\n';
end

% Print all data in dataCell using dataRowCellstrDesc:
dataCellstr = cellfun(@(dataPart, dataDescStr) ...
    sprintf(dataDescStr, dataPart), dataCell, ...
    repmat_as_needed( dataDescCellstr, size(dataCell) ), 'un', 0);

% Fix exponent display:
dataCellstr = cellfun(@(str) ...
    rmexpzeros(str, inpPar.Results.printMode), dataCellstr,'un',0);

% Add header if it should be printed:

% Expand dataCellStr to encompas header row/column:
if ~isempty(headerRowCellstr) && inpPar.Results.printHeaderRow
    dataCellstr = cat(2, cell(size(dataCellstr,1), 1), dataCellstr);
end
if ~isempty(headerColumnCellstr) && inpPar.Results.printHeaderCol
    dataCellstr = cat(1, cell(1, size(dataCellstr,2)), dataCellstr);
end
% Add supplied header strings:
if ~isempty(headerRowCellstr) && inpPar.Results.printHeaderRow
    dataCellstr( (size(dataCellstr,1)-length(headerRowCellstr)+1):end,1)...
        = headerRowCellstr;
end
if ~isempty(headerColumnCellstr) && inpPar.Results.printHeaderCol
    dataCellstr(1, (size(dataCellstr,2)-length(headerColumnCellstr)+1):end) = ...
        headerColumnCellstr;
end
% Note, we allow for too short headers.

% Add empty string to left corner if not specified:
if isempty(dataCellstr{1})
    dataCellstr{1} = '';
end

% Transpose table if specified
if inpPar.Results.transposeTable
    dataCellstr = dataCellstr.';
end

% Check textAlignment:
if length(inpPar.Results.textAlignment) == 1
    tmp_textAlignment = inpPar.Results.textAlignment(ones(1,size(dataCellstr,2)));
elseif length(inpPar.Results.textAlignment) == size(dataCellstr,2)
    tmp_textAlignment = inpPar.Results.textAlignment;
elseif length(inpPar.Results.textAlignment) == size(dataCellstr,2) - 1 && (...
        (inpPar.Results.transposeTable && ~isempty(headerColumnCellstr) ...
        && inpPar.Results.printHeaderCol ) || ( ...
        ~isempty(headerRowCellstr) && inpPar.Results.printHeaderRow && ...
        ~inpPar.Results.transposeTable ) )
    tmp_textAlignment = ['c', inpPar.Results.textAlignment];
else
    if (~isempty(headerColumnCellstr) && inpPar.Results.printHeaderCol && ...
            inpPar.Results.transposeTable ) || ( ...
            ~isempty(headerRowCellstr) && inpPar.Results.printHeaderRow && ...
            ~inpPar.Results.transposeTable )
        strExtra = sprintf(' (or %d)',size(dataCellstr,2)-1);
    else
        strExtra = '';
    end
    warning(['The supplied textAlignemnt string has %d columns, expected ' ...
        '%d%s number of columns. Check if the table is transposed. ' ...
        'Using the first columns value for all columns.'], ...
        length(inpPar.Results.textAlignment), size(dataCellstr,2), strExtra);
    tmp_textAlignment = inpPar.Results.textAlignment(ones(1,size(dataCellstr,2)));
end


% If latex is specified, repeat any \ sign twice, as it should PROBABLY not be
% used as an escape charachter (unless it is already repeated!):
if strcmp(inpPar.Results.printMode,'latex')
    %    dataCellstr = cellfun(@(str) strrep(str,'\','\\'), dataCellstr,'un',0);
    dataCellstr = cellfun(@(str) ...
        regexprep(str,'(?<!\\)\\(?!\\)','\\\\'), dataCellstr,'un',0);
end

% Compute lengths of each column in the dataCellstr:
dataCellstrLength = cellfun(@(str) length(sprintf(str)), dataCellstr);

% Find maximum length in each column:
if inpPar.Results.spaceColPadAlign
    columnMaxLength   = max(dataCellstrLength,[],1) + ...
        2*inpPar.Results.numSpaceColPad;
else
    columnMaxLength   = max(dataCellstrLength,[],1);
end

% Find number of spaces to pad each cell with:
numSpacePad = bsxfun(@minus, columnMaxLength, dataCellstrLength );

% Pad with spaces depending on textAlignment
for iCol = 1:length(tmp_textAlignment)
    if tmp_textAlignment(iCol) == 'r' % inpPar.Results.textAlignment == 'r'
        % Pad to the left
        dataCellstr(:,iCol) = cellfun(@(str, spaceNum) ...
            [repmat(' ',1,spaceNum), str], ...
            dataCellstr(:,iCol), num2cell(numSpacePad(:,iCol)), 'un', 0);
    elseif tmp_textAlignment(iCol) == 'l' % inpPar.Results.textAlignment == 'l'
        % Pad to the right
        dataCellstr(:,iCol) = cellfun(@(str, spaceNum) ...
            [str, repmat(' ',1,spaceNum)], ...
            dataCellstr(:,iCol), num2cell(numSpacePad(:,iCol)), 'un', 0);
    elseif tmp_textAlignment(iCol) == 'c' % inpPar.Results.textAlignment == 'c'
        % Pad equal amount to the left and right:
        dataCellstr(:,iCol) = cellfun(@(str, spaceNum) ...
            [repmat(' ',1,ceil(0.5*spaceNum)), ...
            str, repmat(' ',1,floor(0.5*spaceNum))], ...
            dataCellstr(:,iCol), num2cell(numSpacePad(:,iCol)), 'un', 0);
    end
end

% Put space padding on both sides equal to numColSpacePad:
if ~inpPar.Results.spaceColPadAlign && inpPar.Results.numSpaceColPad > 0
    dataCellstr = cellfun(@(str) ...
        [ repmat(' ', 1, inpPar.Results.numSpaceColPad), ...
        str, repmat(' ', 1, inpPar.Results.numSpaceColPad)], ...
        dataCellstr, 'un', 0);
    % Update maximum length in each column with extra padding:
    columnMaxLength   = columnMaxLength + 2*inpPar.Results.numSpaceColPad;
end

% Add header/data separation line if in text mode:
if strcmp(inpPar.Results.printMode, 'text')
    if inpPar.Results.printHeaderCol && ...
            (~inpPar.Results.transposeTable && ~isempty(headerColumnCellstr) || ...
            (inpPar.Results.transposeTable && ~isempty(headerRowCellstr) ) )
        dataCellstr = cat(1, dataCellstr(1,:), ...
            arrayfun(@(num) repmat(inpPar.Results.rowHSepStr,1,num),...
            columnMaxLength,'un',0), dataCellstr(2:end,:));
    end
    if inpPar.Results.printHeaderRow && ...
            (~inpPar.Results.transposeTable && ~isempty(headerRowCellstr) || ...
            (inpPar.Results.transposeTable && ~isempty(headerColumnCellstr) ) )
        dataCellstr(:,1) = cellfun(@(str) [str inpPar.Results.colHSepStr], ...
            dataCellstr(:,1),'un',0);
    end
end

% Add row separation lines if specified:
if ~isempty(rowSepStr) && strcmp(inpPar.Results.printMode,'text')
    dataCellstrExt = cell(size(dataCellstr)+[size(dataCellstr,1)-1,0]);
    dataCellstrExt(1:2:end,:) = dataCellstr;
    tmpRowSepLines = arrayfun(@(numChar) ...
        repmat(rowSepStr, 1, numChar),columnMaxLength,'un',0);
    if inpPar.Results.printHeaderCol && ...
            (~inpPar.Results.transposeTable && ~isempty(headerColumnCellstr) || ...
            (inpPar.Results.transposeTable && ~isempty(headerRowCellstr) ) )
        dataCellstrExt([2,4],:) = [];
        dataCellstrExt(4:2:end,:) = repmat(tmpRowSepLines,size(dataCellstr,1)-3,1);
    else
        dataCellstrExt(2:2:end,:) = repmat(tmpRowSepLines,size(dataCellstr,1)-1,1);
    end
    dataCellstr = dataCellstrExt;
end

% Add a border around the table:
if strcmp(inpPar.Results.printMode,'text') && inpPar.Results.printBorder
    % Create extended table with border around:
    dataCellstrExt = cell(size(dataCellstr)+[2,2]);
    dataCellstrExt(2:end-1, 2:end-1) = dataCellstr;
    dataCellstrExt(1,2:end-1)   = arrayfun(@(numChar) ...
        repmat(inpPar.Results.borderRowStr, 1, numChar),columnMaxLength,'un',0);
    if inpPar.Results.printHeaderRow && ...
            (~inpPar.Results.transposeTable && ~isempty(headerRowCellstr) || ...
            (inpPar.Results.transposeTable && ~isempty(headerColumnCellstr) ) )
        dataCellstrExt(1,2)   = {[dataCellstrExt{1,2} inpPar.Results.colHSepStr]};
    end
    dataCellstrExt(end,2:end-1) = dataCellstrExt(1,2:end-1);
    dataCellstr = dataCellstrExt;
end

% Generate strings for each row using cellstr2str:
dataCellstrRow = arrayfun(@(idxRow)  cellstr2str(dataCellstr(idxRow,:), ...
    colSepStr), 1:size(dataCellstr,1), 'un', 0);

% Generate strings for full table, break lines etc:
dataCellstrTable = cat(1, dataCellstrRow,...
    repmat({newLineSepStr}, size(dataCellstrRow)) );


% Create the complete table string:
tableStr = sprintf(cat(2, dataCellstrTable{:}));

% Check if latex and if we should add tabular enviroment to string;
if strcmp(inpPar.Results.printMode,'latex')
    if  inpPar.Results.printLatexFull
        % Generate \tabular enviroment to the table
        tmpColSep = cat(1,tmp_textAlignment,repmat('|',1,length(tmp_textAlignment)));
        numColStr         = [ '{|' tmpColSep(:).' '}'];
        beginTabularStr   = sprintf( '\\begin{tabular}%s\\hline\n', numColStr);
        endTabularStr     = sprintf( '\\end{tabular}\n');
        tableStr = cat(2,beginTabularStr,tableStr,endTabularStr);
    end
    % If latex is specified, repeat any \ sign again twice:
    tableStr = strrep(tableStr,'\','\\');
    % Note, this is to ensure that 'fprintf(1, tableStr)' works.
end



% Print table if no output is asked for:
if nargout >= 1
    varargout{1} = tableStr;
else
    fprintf(1,tableStr);
end
end

function repData = repmat_as_needed(dataPart, repDataSize)
%% REPMAT_AS_NEEDED Repeat data in dataPart to match repDataSize if possible
% Syntax:
%     repData = REPMAT_AS_NEEDED(dataPart, repDataSize)
%
% Comment:
%     Repeats data up repDataSize unless it is already of the correct size.

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2014/10/07 9:00:00$

% Check size difference between dataPart and repDataSize, repeat up to size:
dataSize    = size(dataPart);
dataSizeExt = ones(1,length(repDataSize));
dataSizeExt(1:length(dataSize)) = dataSize;
numRepDim   = repDataSize./dataSizeExt;
if all( (abs(round(numRepDim)-numRepDim)) <= eps('double') )
    repData     = repmat(dataPart, round(numRepDim));
else
    error('repmat_as_needed:invalidInput',...
        'The supplied dataPart and repDataSize have incompatible sizes')
end
end

function str = rmexpzeros(str, printMode)
%% RMEXPZEROS Remove extra zeroes after exponent string (i.e. 1.1e+01 remove 0)
%
% Syntax:
%     str = RMEXPZEROS(str, ( printMode='text' ) )
%
% Comment:
%     Removes zero padding in strings containing exponent expressions on the
%     form e+/- or E+/-.  Completely removes e/E if double zeros.
%     If printMode = 'latex', repaces e+/- with ${\times}10^$
%

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2012/?/? 00:00:00$
%   $Revision: 2.0$  $Date: 2014/10/06 13:00:00$
%     Cleaned function script and comments, added latex printMode

if nargin <=1
    printMode = 'text';
end

if strcmp(printMode,'text')
    
    str = strrep(str,'e+00','');
    str = strrep(str,'e+0','e+');
    str = strrep(str,'e-00','');
    str = strrep(str,'e-0','e-');
    
    str = strrep(str,'E+00','');
    str = strrep(str,'E+0','E+');
    str = strrep(str,'E-00','');
    str = strrep(str,'E-0','E-');
    
elseif strcmp(printMode,'latex')
    
    str = strrep(str,'e+00','');
    [startIndex, endIndex]  = regexp(str,'[eE][+](0*)');
    if ~isempty(startIndex)
        if strcmp(str(1:startIndex-1),'1')
            str = ['$10^{' str(endIndex+1:end) '}$'];
        else
            str = ['$' str(1:startIndex-1) ...
                '{\\times}10^{' str(endIndex+1:end) '}$'];
        end
    end
    [startIndex,endIndex]  = regexp(str,'[eE][-](0*)');
    if ~isempty(startIndex)
        if strcmp(str(1:startIndex-1),'1')
            str = ['$10^{-' str(endIndex+1:end) '}$'];
        else
            str = ['$' str(1:startIndex-1) ...
                '{\\times}10^{-' str(endIndex+1:end) '}$'];
        end
    end
end
end

function str = cellstr2str(cellstr, separationStr, numericConversionStr)
%% CELLSTR2STR Convert a cellstring to a single string with a separation string
%
% Syntax:
%     str = cellstr2str(cellstr, separationStr)
%     str = cellstr2str(numericArray, separationStr, numericConversionStr)
%
% Comment:
%     Utility function for writing a cellstr as a single string.
%     Can also convert a numeric array to a string.
%     Main purpose is add separation charachter between the parts of the
%     cellstr.

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2012/?/? 00:00:00$
%   $Revision: 2.0$  $Date: 2014/10/06 13:00:00$
%     Complete overhaul of function

% Default separationStr
if (nargin <= 1)
    separationStr = ' ';
end
% Default numericConversionStr
if nargin <= 2
    numericConversionStr = '%g';
end

% Numeric input array:
if ~iscell(cellstr)
    % Convert numeric array to cellstr:
    cellstr = cellfun(@(innum) num2str(innum, numericConversionStr),...
        num2cell(cellstr),'un',0);
end

% Concatenate into output:
cellstr        = cellstr(:).';
len_cs         = length(cellstr);
cellsty        = cat(2, repmat( { separationStr }, [1, len_cs-1]), {''});
cellstr_merge  = cat(1, cellstr, cellsty );
str            = cat(2, cellstr_merge{:} ) ;
end

%% Get Full Path
function File = GetFullPath(File, Style)
% GetFullPath - Get absolute canonical path of a file or folder
% Absolute path names are safer than relative paths, when e.g. a GUI or TIMER
% callback changes the current directory. Only canonical paths without "." and
% ".." can be recognized uniquely.
% Long path names (>259 characters) require a magic initial key "\\?\" to be
% handled by Windows API functions, e.g. for Matlab's FOPEN, DIR and EXIST.
%
% FullName = GetFullPath(Name, Style)
% INPUT:
%   Name:  String or cell string, absolute or relative name of a file or
%          folder. The path need not exist. Unicode strings, UNC paths and long
%          names are supported.
%   Style: Style of the output as string, optional, default: 'auto'.
%          'auto': Add '\\?\' or '\\?\UNC\' for long names on demand.
%          'lean': Magic string is not added.
%          'fat':  Magic string is added for short names also.
%          The Style is ignored when not running under Windows.
%
% OUTPUT:
%   FullName: Absolute canonical path name as string or cell string.
%          For empty strings the current directory is replied.
%          '\\?\' or '\\?\UNC' is added on demand.
%
% NOTE: The M- and the MEX-version create the same results, the faster MEX
%   function works under Windows only.
%   Some functions of the Windows-API still do not support long file names.
%   E.g. the Recycler and the Windows Explorer fail even with the magic '\\?\'
%   prefix. Some functions of Matlab accept 260 characters (value of MAX_PATH),
%   some at 259 already. Don't blame me.
%   The 'fat' style is useful e.g. when Matlab's DIR command is called for a
%   folder with les than 260 characters, but together with the file name this
%   limit is exceeded. Then "dir(GetFullPath([folder, '\*.*], 'fat'))" helps.
%
% EXAMPLES:
%   cd(tempdir);                    % Assumed as 'C:\Temp' here
%   GetFullPath('File.Ext')         % 'C:\Temp\File.Ext'
%   GetFullPath('..\File.Ext')      % 'C:\File.Ext'
%   GetFullPath('..\..\File.Ext')   % 'C:\File.Ext'
%   GetFullPath('.\File.Ext')       % 'C:\Temp\File.Ext'
%   GetFullPath('*.txt')            % 'C:\Temp\*.txt'
%   GetFullPath('..')               % 'C:\'
%   GetFullPath('..\..\..')         % 'C:\'
%   GetFullPath('Folder\')          % 'C:\Temp\Folder\'
%   GetFullPath('D:\A\..\B')        % 'D:\B'
%   GetFullPath('\\Server\Folder\Sub\..\File.ext')
%                                   % '\\Server\Folder\File.ext'
%   GetFullPath({'..', 'new'})      % {'C:\', 'C:\Temp\new'}
%   GetFullPath('.', 'fat')         % '\\?\C:\Temp\File.Ext'
%
% COMPILE:
%   Automatic: InstallMex GetFullPath.c uTest_GetFullPath
%   Manual:    mex -O GetFullPath.c
%   Download:  http://www.n-simon.de/mex
% Run the unit-test uTest_GetFullPath after compiling.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
% Assumed Compatibility: higher Matlab versions
% Author: Jan Simon, Heidelberg, (C) 2009-2013 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also: CD, FULLFILE, FILEPARTS.

% $JRev: R-G V:032 Sum:7Xd/JS0+yfax Date:15-Jan-2013 01:06:12 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_GetFullPath $
% $File: Tools\GLFile\GetFullPath.m $
% History:
% 001: 20-Apr-2010 22:28, Successor of Rel2AbsPath.
% 010: 27-Jul-2008 21:59, Consider leading separator in M-version also.
% 011: 24-Jan-2011 12:11, Cell strings, '~File' under linux.
%      Check of input types in the M-version.
% 015: 31-Mar-2011 10:48, BUGFIX: Accept [] as input as in the Mex version.
%      Thanks to Jiro Doke, who found this bug by running the test function for
%      the M-version.
% 020: 18-Oct-2011 00:57, BUGFIX: Linux version created bad results.
%      Thanks to Daniel.
% 024: 10-Dec-2011 14:00, Care for long names under Windows in M-version.
%      Improved the unittest function for Linux. Thanks to Paul Sexton.
% 025: 09-Aug-2012 14:00, In MEX: Paths starting with "\\" can be non-UNC.
%      The former version treated "\\?\C:\<longpath>\file" as UNC path and
%      replied "\\?\UNC\?\C:\<longpath>\file".
% 032: 12-Jan-2013 21:16, 'auto', 'lean' and 'fat' style.

% Initialize: ==================================================================
% Do the work: =================================================================

% #############################################
% ### USE THE MUCH FASTER MEX ON WINDOWS!!! ###
% #############################################

% Difference between M- and Mex-version:
% - Mex does not work under MacOS/Unix.
% - Mex calls Windows API function GetFullPath.
% - Mex is much faster.

% Magix prefix for long Windows names:
if nargin < 2
    Style = 'auto';
end

% Handle cell strings:
% NOTE: It is faster to create a function @cell\GetFullPath.m under Linux, but
% under Windows this would shadow the fast C-Mex.
if isa(File, 'cell')
    for iC = 1:numel(File)
        File{iC} = GetFullPath(File{iC}, Style);
    end
    return;
end

% Check this once only:
isWIN    = strncmpi(computer, 'PC', 2);
MAX_PATH = 260;

% Warn once per session (disable this under Linux/MacOS):
persistent hasDataRead
if isempty(hasDataRead)
    % Test this once only - there is no relation to the existence of DATAREAD!
    %if isWIN
    %   Show a warning, if the slower Matlab version is used - commented, because
    %   this is not a problem and it might be even useful when the MEX-folder is
    %   not inlcuded in the path yet.
    %   warning('JSimon:GetFullPath:NoMex', ...
    %      ['GetFullPath: Using slow Matlab-version instead of fast Mex.', ...
    %       char(10), 'Compile: InstallMex GetFullPath.c']);
    %end
    
    % DATAREAD is deprecated in 2011b, but still available. In Matlab 6.5, REGEXP
    % does not know the 'split' command, therefore DATAREAD is preferred:
    hasDataRead = ~isempty(which('dataread'));
end

if isempty(File)  % Accept empty matrix as input:
    if ischar(File) || isnumeric(File)
        File = cd;
        return;
    else
        error(['JSimon:', mfilename, ':BadTypeInput1'], ...
            ['*** ', mfilename, ': Input must be a string or cell string']);
    end
end

if ischar(File) == 0  % Non-empty inputs must be strings
    error(['JSimon:', mfilename, ':BadTypeInput1'], ...
        ['*** ', mfilename, ': Input must be a string or cell string']);
end

if isWIN  % Windows: --------------------------------------------------------
    FSep = '\';
    File = strrep(File, '/', FSep);
    
    % Remove the magic key on demand, it is appended finally again:
    if strncmp(File, '\\?\', 4)
        if strncmpi(File, '\\?\UNC\', 8)
            File = ['\', File(7:length(File))];  % Two leading backslashes!
        else
            File = File(5:length(File));
        end
    end
    
    isUNC   = strncmp(File, '\\', 2);
    FileLen = length(File);
    if isUNC == 0                        % File is not a UNC path
        % Leading file separator means relative to current drive or base folder:
        ThePath = cd;
        if File(1) == FSep
            if strncmp(ThePath, '\\', 2)   % Current directory is a UNC path
                sepInd  = strfind(ThePath, '\');
                ThePath = ThePath(1:sepInd(4));
            else
                ThePath = ThePath(1:3);     % Drive letter only
            end
        end
        
        if FileLen < 2 || File(2) ~= ':'  % Does not start with drive letter
            if ThePath(length(ThePath)) ~= FSep
                if File(1) ~= FSep
                    File = [ThePath, FSep, File];
                else                        % File starts with separator:
                    File = [ThePath, File];
                end
            else                           % Current path ends with separator:
                if File(1) ~= FSep
                    File = [ThePath, File];
                else                        % File starts with separator:
                    ThePath(length(ThePath)) = [];
                    File = [ThePath, File];
                end
            end
            
        elseif FileLen == 2 && File(2) == ':'   % "C:" current directory on C!
            % "C:" is the current directory on the C-disk, even if the current
            % directory is on another disk! This was ignored in Matlab 6.5, but
            % modern versions considers this strange behaviour.
            if strncmpi(ThePath, File, 2)
                File = ThePath;
            else
                try
                    File = cd(cd(File));
                catch    % No MException to support Matlab6.5...
                    if exist(File, 'dir')  % No idea what could cause an error then!
                        rethrow(lasterror);
                    else  % Reply "K:\" for not existing disk:
                        File = [File, FSep];
                    end
                end
            end
        end
    end
    
else         % Linux, MacOS: ---------------------------------------------------
    FSep = '/';
    File = strrep(File, '\', FSep);
    
    if strcmp(File, '~') || strncmp(File, '~/', 2)  % Home directory:
        HomeDir = getenv('HOME');
        if ~isempty(HomeDir)
            File(1) = [];
            File    = [HomeDir, File];
        end
        
    elseif strncmpi(File, FSep, 1) == 0
        % Append relative path to current folder:
        ThePath = cd;
        if ThePath(length(ThePath)) == FSep
            File = [ThePath, File];
        else
            File = [ThePath, FSep, File];
        end
    end
end

% Care for "\." and "\.." - no efficient algorithm, but the fast Mex is
% recommended at all!
if ~isempty(strfind(File, [FSep, '.']))
    if isWIN
        if strncmp(File, '\\', 2)  % UNC path
            index = strfind(File, '\');
            if length(index) < 4    % UNC path without separator after the folder:
                return;
            end
            Drive            = File(1:index(4));
            File(1:index(4)) = [];
        else
            Drive     = File(1:3);
            File(1:3) = [];
        end
    else  % Unix, MacOS:
        isUNC   = false;
        Drive   = FSep;
        File(1) = [];
    end
    
    hasTrailFSep = (File(length(File)) == FSep);
    if hasTrailFSep
        File(length(File)) = [];
    end
    
    if hasDataRead
        if isWIN  % Need "\\" as separator:
            C = dataread('string', File, '%s', 'delimiter', '\\');  %#ok<REMFF1>
        else
            C = dataread('string', File, '%s', 'delimiter', FSep);  %#ok<REMFF1>
        end
    else  % Use the slower REGEXP, when DATAREAD is not available anymore:
        C = regexp(File, FSep, 'split');
    end
    
    % Remove '\.\' directly without side effects:
    C(strcmp(C, '.')) = [];
    
    % Remove '\..' with the parent recursively:
    R = 1:length(C);
    for dd = reshape(find(strcmp(C, '..')), 1, [])
        index    = find(R == dd);
        R(index) = [];
        if index > 1
            R(index - 1) = [];
        end
    end
    
    if isempty(R)
        File = Drive;
        if isUNC && ~hasTrailFSep
            File(length(File)) = [];
        end
        
    elseif isWIN
        % If you have CStr2String, use the faster:
        %   File = CStr2String(C(R), FSep, hasTrailFSep);
        File = sprintf('%s\\', C{R});
        if hasTrailFSep
            File = [Drive, File];
        else
            File = [Drive, File(1:length(File) - 1)];
        end
        
    else  % Unix:
        File = [Drive, sprintf('%s/', C{R})];
        if ~hasTrailFSep
            File(length(File)) = [];
        end
    end
end

% "Very" long names under Windows:
if isWIN
    if ~ischar(Style)
        error(['JSimon:', mfilename, ':BadTypeInput2'], ...
            ['*** ', mfilename, ': Input must be a string or cell string']);
    end
    
    if (strncmpi(Style, 'a', 1) && length(File) >= MAX_PATH) || ...
            strncmpi(Style, 'f', 1)
        % Do not use [isUNC] here, because this concerns the input, which can
        % '.\File', while the current directory is an UNC path.
        if strncmp(File, '\\', 2)  % UNC path
            File = ['\\?\UNC', File(2:end)];
        else
            File = ['\\?\', File];
        end
    end
end

% return;
end