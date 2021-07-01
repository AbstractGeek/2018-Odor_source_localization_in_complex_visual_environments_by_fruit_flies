% A script file to extract speed change data from the whole dataset (mat
% file)

%% Initialize
% Defaults
odor_plume = 1.0;   % 1cm dia, 0.5 cm (5 mm) radius
plume_contact_cutoff = 4; % 4cm before landing

% variable initialization
qualified_files = cell(1000,1);
qualified_ind = 1;
data = struct;
data.trajData = cell(1000,1);

% Load the airflow mat file
load analyzeddata_cutoff30.mat;

% Initialize output file
logfile = 'Selected Trajectories.org';   % Save as an org file (easier to read with emacs/spacemacs)
if exist(logfile,'file')==2
   warning('logfile already present in the folder. New selections will be appended to the file with a new heading');      
end

% Open logfile and enter default stuff
fid = fopen(logfile,'a');   % Disregard contents if any
fprintf(fid,'* [%s]\n',datestr(datetime('now')));
fprintf(fid, '** Trajectory Selection Criteria\n');
fprintf(fid, 'odor plume width = %g\nplume contact cutoff=%d\n\n',...
    odor_plume, plume_contact_cutoff);
fprintf(fid, '** Selected Files\n');

%% Filter trajectories with encounter before cutoff
treatments = fieldnames(assay);
for i=1:length(treatments)
   trajData = assay.(treatments{i}).trajData;
   samplesize = size(trajData,1);
   
   for j=1:samplesize       
       % find location of first odor encounter
       encounter_ind = find(trajData{j}.globalpara.OdorAxisDist<=odor_plume, 1, 'first');
       % check if it is greater than the cutoff
       if (trajData{j}.globalpara.SourceDist(encounter_ind) > plume_contact_cutoff)           
           % save trajectory data
           data.trajData{qualified_ind, 1} = trajData{j};           
           % save the name to the list of files and increment index
           qualified_files{qualified_ind, 1} = trajData{j}.name;
           fprintf(fid, '- %s\n',trajData{j}.name);
           qualified_ind = qualified_ind+1;           
       end
   end
    
end

%% Clean up and save data
% clean up cell arrays
qualified_files(cellfun(@isempty, qualified_files)) = [];
data.trajData(cellfun(@isempty, data.trajData)) = [];

% copy basic information into the new structure
data.treatment = sprintf('contact%d',plume_contact_cutoff);
data.sam_freq = assay.(treatments{1}).sam_freq;
data.cutoff_freq = assay.(treatments{1}).cutoff_freq;
data.landingPoint = assay.(treatments{1}).landingPoint;
data.odor_plume = odor_plume;
data.plume_contact_cutoff = plume_contact_cutoff;
data.qualified_files = qualified_files;
assay = struct;
assay.(strcat('t',data.treatment)) = data;

% matname
matname = sprintf('PlumeEncounterDataset[pw-%g][cutoff=%d].mat',...
    odor_plume,plume_contact_cutoff);

% save the array
save(matname,'assay');