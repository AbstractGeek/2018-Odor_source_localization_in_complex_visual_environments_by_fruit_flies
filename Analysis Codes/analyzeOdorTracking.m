function [assay] = analyzeOdorTracking()
%'The' function to analyze odor tracking!!
% It takes in csv files as input and processes it to generate required parameters that would reveal the algorithm of odor tracking in Drosophila
% Melanogaster. It outputs a mat files for each treatment from in the same folder as the data. It also outputs a complete mat file with all data,
% which could be stored anywhere
%
% [Inputs]: A lot of inputs will be take in from the GUI. Make sure you read the instructions in the GUI.
%
% [Outputs]: The code generates a mat file, which contains a nested structure that contains parameters that can used to generate relevant
% plots. The idea is to combine all revelant data of each experimental treatment in a single mat file, which can be used later to do more analysis.
%
% StandAlone function. Check if it has the following functions:
% Datapoints.m
% DistVel.m
% FindOverallSum.m
% FindOverallAverage.m
% angular_velocity.m
% StartBinning.m
% DistBinner.m
% TimeBinner.m
% Tortuosity.m
%
%
% Version 3.5
% Last Modified sometime in Oct 2013
% by Dinesh Natesan
%
% Updated - 31 Aug 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Gathering

% Take in the number of experimental treatments to be analysed
num = listdlg('Name','Details..Gimme Details','ListSize',[300 300],'PromptString','Select Number of treatments to analyze','SelectionMode','single','ListString',int2str((1:30)'));

% Define basefields
field = {'','1','1','1','0.1','2'};

% Define Constants
% Addition: 25/12/2014
% HoverTime = movement of 5% of body length per wingbeat
HoverSpeed = 5/100*0.3*250; % 0.3 cm - body length of the fly. 250Hz - wingbeat frequency of the fly
% SlowMovementTime = movement of 10% of body length per wingbeat
SlowSpeed = 20/100*0.3*250; % 0.3 cm - body length of the fly. 250Hz - wingbeat frequency of the fly

% First Compute and save data in structures.

for m = 1:num
    
    %% Initialize required parameters
    
    % Get Path (aka folder name)
    specpath = uigetdir;
    
    % Get the lastname
    [parentdir, data.treatment, ~] = fileparts(specpath);
    
    field{1,1}=data.treatment;
    
    % Take in the Name of the Experimental treatment, Sampling rate and cut-off frequency.
    % answer = inputdlg({'Name','Source Distance Bin','Odor Axis Distance Bin','Time Bin','Bin sourcedist based on point of Landing?(Yes = 1, No = 2)'},'Experiment Details',1,{data.treatment,'1','1','0.1','2'});           
    
%         answer = inputdlg({'Name','Source Distance Bin','Odor Axis Distance Bin','Rectangular Bin','Time Bin','Bin sourcedist based on point of Landing?(Yes = 1, No = 2)'},'Experiment Details',1,{data.treatment,'1','1','1','0.1','2'});
    answer = inputdlg({'Name','Source Distance Bin','Odor Axis Distance Bin','Rectangular Bin','Time Bin','Bin sourcedist based on point of Landing?(Yes = 1, No = 2)'},'Experiment Details',1,field);
    if strcmp(answer,'Cancel')
        return;
    end
    if m==1
        field = answer;
    end    
        
    % Save it in the structure 'data'
    data.treatment  = char(answer{1,1});
    disp(data.treatment);
    data.binValues.SourceDistBin = str2double(answer{2,1});
    data.binValues.OdorAxisBin = str2double(answer{3,1});
    data.binValues.RectBin = str2double(answer{4,1});
    data.binValues.TimeBin = str2double(answer{5,1});
    data.binBasedOnLand = str2double(answer{6,1});
    
    % Set constants, change if necessary. (Only Have to change here, the code will take care of the rest)
    %%% Filtering Constants
    % data.sam_freq = 100;
    % data.cutoff_freq = 45;
    data.sam_freq = 100;
    data.cutoff_freq = 30;
    disp(strcat('Sampling Rate (fps):',num2str(data.sam_freq)));
    disp(strcat('CutOff Freq (Hz):',num2str(data.cutoff_freq)));
    %%% Bin limits
    data.binLimits.SourceDistBin = 8;
    data.binLimits.OdorAxisBin = 4;
    data.binLimits.RectBin = 8;
    
    %%% Bin sizes
    data.binNumbers.SourceDistBinNumber = ceil(data.binLimits.SourceDistBin./data.binValues.SourceDistBin);
    data.binNumbers.OdorAxisBinNumber = ceil(data.binLimits.OdorAxisBin./data.binValues.OdorAxisBin);
    data.binNumbers.RectBinNumber = ceil(data.binLimits.RectBin./data.binValues.RectBin);
    
    % Get list of csv files in it.
    specfile = dir(strcat(specpath,'/*.csv'));
    
    % Find the number of files in the directory
    sel = size(specfile,1);
    % Initalize cell variables
    data.trajData = cell(sel,1);
    % Intialize other structures
    
    
    % Start processing data and save them into the structure
    for i=1:1:sel
        
        %% Find parameters, bin and save it in the data structure
        % To make selecting convienent, get into the specpath directory
        cd(specpath);
        % Save name of the file and trajectory coordinates in the structure
        data.trajData{i,1}.name=specfile(i,1).name;
        % Change the below command later
        data.trajData{i,1}.expdata.unfiltdata = dlmread(strcat(specpath,'/',specfile(i,1).name),',',0,0);
        % Extract cs,bs and remove NaN's.
        try
            [data.trajData{i,1}.expdata.pts,data.trajData{i,1}.expdata.cs,data.trajData{i,1}.expdata.bs,data.trajData{i,1}.binnerInput.SourceDist,...
                data.landingPoint]=Datapoints(data.trajData{i,1}.expdata.unfiltdata,data.sam_freq,data.cutoff_freq,data.binBasedOnLand);
        catch
            disp(specfile(i,1).name);
            error('errormsg');
        end
        
        % Added 28/09/2014. Changed DistVel to calculate odor axis based on
        % landing
        % Start calculating first order parameters aka global parameters
        [data.trajData{i,1}.globalpara.Distance,data.trajData{i,1}.binnerInput.OdorAxisDist,data.trajData{i,1}.globalpara.CastDist,data.trajData{i,1}.globalpara.SurgeDist,...
            data.trajData{i,1}.globalpara.FlySpeed,data.trajData{i,1}.globalpara.CastSpeed,data.trajData{i,1}.globalpara.SurgeSpeed,data.trajData{i,1}.globalpara.CurveArea,...
            data.trajData{i,1}.globalpara.HorzArea,data.trajData{i,1}.globalpara.VertArea] = DistVel(data.trajData{i,1}.expdata.pts,data.trajData{i,1}.expdata.cs,data.sam_freq,data.landingPoint);
        
        % % Addition: 15/06/2015: Add a new global para - PlumeCast
        data.trajData{i,1}.globalpara.PlumeCast = abs(diff(data.trajData{i,1}.binnerInput.OdorAxisDist));
        
        % Set the rectangular bin source (just the x coordinate) - based on landing if the option is selected
        % X coordinate is negative, so making it absolute
        data.trajData{i,1}.binnerInput.RectDist = abs(data.trajData{i,1}.expdata.pts(:,1) - data.trajData{i,1}.expdata.cs(data.landingPoint,1));

        % Add SourceDist and OdorAxis Dist to the global para list
        data.trajData{i,1}.globalpara.SourceDist = data.trajData{i,1}.binnerInput.SourceDist;
        data.trajData{i,1}.globalpara.OdorAxisDist = data.trajData{i,1}.binnerInput.OdorAxisDist;
        
        % Find OverallSum
        data.trajData{i,1}.OverallSum = FindOverallSum(data.trajData{i,1}.globalpara,data.trajData{i,1}.expdata.pts,data.sam_freq);
        % Find OverallAverage
        data.trajData{i,1}.OverallAverage = FindOverallAverage(data.trajData{i,1}.globalpara,data.trajData{i,1}.expdata.pts,data.sam_freq);
        
        % Addition: 25/12/2014
        % HoverDuration and SlowSpeed Duration calculation.
        data.trajData{i,1}.OverallSum.HoverDuration = sum(data.trajData{i,1}.globalpara.FlySpeed<=HoverSpeed)./data.sam_freq;
        data.trajData{i,1}.OverallSum.SlowSpeedDuration = sum(data.trajData{i,1}.globalpara.FlySpeed<=SlowSpeed)./data.sam_freq;
        % Addition: 15/06/2015
        data.trajData{i,1}.OverallSum.PlumeCast = sum(data.trajData{i,1}.globalpara.PlumeCast);
        data.trajData{i,1}.OverallAverage.PlumeCast = mean(data.trajData{i,1}.globalpara.PlumeCast);
        
        
        % Angular velocity commented out for now
        [data.trajData{i,1}.globalpara.AngVel,data.trajData{i,1}.globalpara.TurnAng] = angular_velocity(data.trajData{i,1}.expdata.pts,data.sam_freq);
        
        % Add angular velocity to the overall average (31/08/16)
        data.trajData{i,1}.OverallAverage.AngVel = mean(data.trajData{i,1}.globalpara.AngVel);
        
        % Add more parameters here
        
        % Bin Data (Distance, OdorAxis, RectBin and Time bins)
        data.trajData{i,1}.binnedData = StartBinning (data.trajData{i,1}.binnerInput,data.trajData{i,1}.globalpara,data.binValues,[data.binLimits.SourceDistBin;data.binLimits.OdorAxisBin;data.binLimits.RectBin],data.sam_freq);
        % Find tortuosity
        [data.trajData{i,1}.binnedData.TimeBin.Tortuosity.zone,data.trajData{i,1}.binnedData.TimeBin.Tortuosity.zoneMean,data.trajData{i,1}.OverallSum.Tortuosity] = Tortuosity(data.trajData{i,1}.globalpara.Distance,data.trajData{i,1}.expdata.pts,data.binValues.TimeBin,data.sam_freq);
        
        % Addition: 25/12/2014
        % HoverDuration,SlowSpeedDuration
        % Time Spent - SourceDistBin
        data.trajData{i,1}.binnedData.SourceDistBin.TimeSpent.zoneSum = NaN(size(data.trajData{i,1}.binnedData.SourceDistBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.SourceDistBin.HoverDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.SourceDistBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.SourceDistBin.SlowSpeedDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.SourceDistBin.Distance.zone,1),1);
        for j=1:size(data.trajData{i,1}.binnedData.SourceDistBin.Distance.zone,1)
            data.trajData{i,1}.binnedData.SourceDistBin.TimeSpent.zoneSum(j,1) = size(data.trajData{i,1}.binnedData.SourceDistBin.Distance.zone{j,1},1)./data.sam_freq;
            data.trajData{i,1}.binnedData.SourceDistBin.HoverDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.SourceDistBin.FlySpeed.zone{j,1}<=HoverSpeed)./data.sam_freq;
            data.trajData{i,1}.binnedData.SourceDistBin.SlowSpeedDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.SourceDistBin.FlySpeed.zone{j,1}<SlowSpeed)./data.sam_freq;
        end
        % Time Spent - OdorAxisBin
        data.trajData{i,1}.binnedData.OdorAxisBin.TimeSpent.zoneSum = NaN(size(data.trajData{i,1}.binnedData.OdorAxisBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.OdorAxisBin.HoverDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.OdorAxisBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.OdorAxisBin.SlowSpeedDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.OdorAxisBin.Distance.zone,1),1);
        for j=1:size(data.trajData{i,1}.binnedData.OdorAxisBin.Distance.zone,1)
            data.trajData{i,1}.binnedData.OdorAxisBin.TimeSpent.zoneSum(j,1) = size(data.trajData{i,1}.binnedData.OdorAxisBin.Distance.zone{j,1},1)./data.sam_freq;
            data.trajData{i,1}.binnedData.OdorAxisBin.HoverDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.OdorAxisBin.FlySpeed.zone{j,1}<=HoverSpeed)./data.sam_freq;
            data.trajData{i,1}.binnedData.OdorAxisBin.SlowSpeedDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.OdorAxisBin.FlySpeed.zone{j,1}<SlowSpeed)./data.sam_freq;
        end        
        % Time Spent - RectBin
        data.trajData{i,1}.binnedData.RectBin.TimeSpent.zoneSum = NaN(size(data.trajData{i,1}.binnedData.RectBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.RectBin.HoverDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.RectBin.Distance.zone,1),1);
        data.trajData{i,1}.binnedData.RectBin.SlowSpeedDuration.zoneSum = NaN(size(data.trajData{i,1}.binnedData.RectBin.Distance.zone,1),1);
        for j=1:size(data.trajData{i,1}.binnedData.RectBin.Distance.zone,1)
            data.trajData{i,1}.binnedData.RectBin.TimeSpent.zoneSum(j,1) = size(data.trajData{i,1}.binnedData.RectBin.Distance.zone{j,1},1)./data.sam_freq;
            data.trajData{i,1}.binnedData.RectBin.HoverDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.RectBin.FlySpeed.zone{j,1}<=HoverSpeed)./data.sam_freq;
            data.trajData{i,1}.binnedData.RectBin.SlowSpeedDuration.zoneSum(j,1) = sum(data.trajData{i,1}.binnedData.RectBin.FlySpeed.zone{j,1}<SlowSpeed)./data.sam_freq;
        end
        % Done
        
        % Binned Data, Treatment vice
        binNames = fieldnames(data.binValues);
        
        for j=1:size(binNames(1:end-1),1)            
            data.sortedBins.(char(binNames(j))).TimeSpent(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).TimeSpent.zoneSum; % Time spent per bin
            % HoverDuration and SlowSpeedDuration Addition (25/12/2014)
            data.sortedBins.(char(binNames(j))).HoverDuration(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).HoverDuration.zoneSum; % HoverDuration per bin
            data.sortedBins.(char(binNames(j))).SlowSpeedDuration(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).SlowSpeedDuration.zoneSum; % SlowSpeedDuration per bin            
            data.sortedBins.(char(binNames(j))).Distance(:,i)=data.trajData{i,1}.binnedData.(char(binNames(j))).Distance.zoneSum; % Distance travelled per bin
            data.sortedBins.(char(binNames(j))).AvgSpeed(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).Distance.zoneSum ./data.trajData{i,1}.binnedData.(char(binNames(j))).TimeSpent.zoneSum;   % Avg Speed per bin
            data.sortedBins.(char(binNames(j))).FlySpeed(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).FlySpeed.zoneMean;
            data.sortedBins.(char(binNames(j))).CastDist(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).CastDist.zoneSum; % Cast Dist(Sum) per bin
            data.sortedBins.(char(binNames(j))).SurgeDist(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).SurgeDist.zoneSum; % Surge Dist(Sum) per bin
            data.sortedBins.(char(binNames(j))).CastSpeed(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).CastSpeed.zoneMean;    % Average Cast Speed per bin
            data.sortedBins.(char(binNames(j))).SurgeSpeed(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).SurgeSpeed.zoneMean;    % Average Surge Speed per bin
            data.sortedBins.(char(binNames(j))).CurveArea(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).CurveArea.zoneSum;    % Total Curve Area per bin
            data.sortedBins.(char(binNames(j))).HorzArea(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).HorzArea.zoneSum;  % Total Horz Area per bin
            data.sortedBins.(char(binNames(j))).VertArea(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).VertArea.zoneSum;  % Total Vert Area per bin
            data.sortedBins.(char(binNames(j))).SourceDist(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).SourceDist.zoneMean; % Average Source Dist per bin
            data.sortedBins.(char(binNames(j))).OdorAxisDist(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).OdorAxisDist.zoneMean; % Average OdorAxisDist per bin   
            % Addition: 15/06/2015 - PlumeCast
            data.sortedBins.(char(binNames(j))).PlumeCast(:,i) = data.trajData{i,1}.binnedData.(char(binNames(j))).PlumeCast.zoneSum; % Total casting per bin   
            % Addition: 31/06/2016 - AngVel, TurnAng
            data.sortedBins.(char(binNames(j))).AngVel(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).AngVel.zoneMean; % Average AngVel per bin   
            data.sortedBins.(char(binNames(j))).TurnAng(:,i) =  data.trajData{i,1}.binnedData.(char(binNames(j))).TurnAng.zoneMean; % Average TurnAng per bin   
        end
        
        % Need to add timebins here
        
        % Add incrementally add overallSums to the bin
        sumNames = fieldnames(data.trajData{i,1}.OverallSum);

        for j=1:size(sumNames)
           data.OverallSum.(char(sumNames(j)))(i,1)=data.trajData{i,1}.OverallSum.(char(sumNames(j)));          
        end
        
        % Add incrementally add overallAverages to the bin
        avgNames = fieldnames(data.trajData{i,1}.OverallAverage);
        for j=1:size(avgNames)
            data.OverallAverage.(char(avgNames(j)))(i,1)=data.trajData{i,1}.OverallAverage.(char(avgNames(j)));
        end
        
        % Add experimental parameters
%         [data.trajData{i,1}.expPara.vAngDeg,data.trajData{i,1}.expPara.vAngRad] = visualAngle(data.trajData{i,1}.expdata.pts,data.trajData{i,1}.expdata.cs);
        
    end
    % For further processing, if required
    assay.(strcat('t',data.treatment))= data;
    
    % Save individual mat files with the data
    matname = fullfile(specpath, strcat('t',data.treatment,'_distBin', int2str(data.binValues.SourceDistBin),...
        '_odorBin',int2str(data.binValues.OdorAxisBin),'_timeBin',num2str(data.binValues.TimeBin),'_cutoff',num2str(data.cutoff_freq),'.mat'));
    save(matname,'-struct','assay',strcat('t',data.treatment));
    cutoffFreq = data.cutoff_freq;        
    cd(parentdir);
    clear data;
    
    
end

% Take in directory of storage and create.
destpath = uigetdir(pwd,'Duuuuude!! Where do you want to save the whole mat file????');
% Save final Mat files
save(fullfile(destpath, strcat('analyzeddata','_cutoff',...
        num2str(cutoffFreq), '.mat')),'assay');

end


function [pts,cs,bs,SourceDist,landingPoint] = Datapoints(hspecdata,sam,cutoff,binBasedOnLand)

% This function removes NaN, filters the trajectory, trims it from the first point of entry into the 8cm zone and gives outputs them.
%
% [Inputs]:
% hspecdata -> hspecdata of the file
% sam -> sampling rate of the videos
% cutoff -> cutoff frequency of the filter
%
% [Outputs]:
% pts -> [x,y,z] points
% cs -> [x,y,z] points of centres, in sequential order
% bs -> [x,y,z] points of bases, in sequential order
% SourceDist -> Displacement from the odor source
%
% NOTE : Not doing the usual NaN removal since the data is already modified!!
%
% Version 2.02
% Last Modified on 28 Dec 2012
% by Dinesh Natesan

% Initialize required variables
col = size(hspecdata,2);
obj = (col-3)/6;

% Assign x,y,z points of the trajectory
pts = hspecdata(:,1:3);

% Initialise and Extract centre and capillary base coordinates from hspecdata
% Initialising
cs=zeros(obj,3);
bs=zeros(obj,3);

% Extracting centre and base coordinates
for n=1:obj
    i=3*n+1;
    j=3*(obj+n)+1;
    cs(n,:)=hspecdata(1,i:i+2);
    bs(n,:)=hspecdata(1,j:j+2);
end

% Filter the trajectory, if necessary
if (sam~=0 && cutoff~=0)
    pts = ButterFilt(pts,sam,cutoff);
end

% Not doing the cutting, since the data is already cut
% Take the trajectory of the fly, from the moment the 8cm zone around the visual cue. (near field cue requirement)
% Displacement from the odor source.

if (binBasedOnLand == 1) 
    [~,landingPoint]= min(sqrt((cs(:,1)-pts(end,1)).^2+(cs(:,2)-pts(end,2)).^2+(cs(:,3)-pts(end,3)).^2));
    SourceDist = sqrt((pts(:,1)-cs(landingPoint,1)).^2+(pts(:,2)-cs(landingPoint,2)).^2+(pts(:,3)-cs(landingPoint,3)).^2);
else
    landingPoint = 1;
    SourceDist = sqrt((pts(:,1)-cs(1,1)).^2+(pts(:,2)-cs(1,2)).^2+(pts(:,3)-cs(1,3)).^2);
end

end

function [Distance,OdorAxisDist,CastDist,SurgeDist,FlySpeed,CastSpeed,SurgeSpeed,CurveArea,HorzArea,VertArea] = DistVel(pts,cs,sam,landingPoint)

% This function finds the 'global parameters' from the trajectory data
%
% Inputs :
% pts -> [x,y,z] points
% cs -> [x,y,z] points of centres, in sequential order
% sam -> sampling rate (frame rate)
%
% Outputs :
% Distance -> Distance travelled by the fly between two frames
% OdorAxisDist -> Fly's distance from the odor axis (perpendicular distance).
% CastDist -> Distance travelled in the yz plane (between two frames)
% SurgeDist -> Distance travelled in the x plane (between two frames)
% FlySpeed -> Speed of the fly (Distance travelled / delta_t)
% CastSpeed -> Casting speed -> y-z plane speed (CastDist/delta_t)
% SurgeSpeed -> Surging speed -> x axis speed (SurgeDist/delta_t)
% CurveArea -> dOdorAxisDist*dx
% HorzArea -> dy*dx
% VertArea -> dz*dx
%
% Version 2.6
% Last Modified on 19th March 2013
% by Dinesh Natesan

% Distance
Distance = sqrt((diff(pts(:,1)).^2)+(diff(pts(:,2)).^2)+(diff(pts(:,3)).^2));

% Fly's distance from the odor axis (perpendicular distance).
OdorAxisDist = sqrt((pts(:,2)-cs(landingPoint,2)).^2 + (pts(:,3)-cs(landingPoint,3)).^2);

% Speed of the fly
FlySpeed= Distance.*sam;

% Area under the curve
dp=(OdorAxisDist(1:(size(OdorAxisDist)-1))+OdorAxisDist(2:size(OdorAxisDist)))/2;
dy = pts(:,2)-cs(landingPoint,2);
dy = (dy(1:(size(dy)-1))+ dy(2:size(dy)))/2;
dz = pts(:,3)-cs(landingPoint,3);
dz = (dz(1:(size(dz)-1))+ dz(2:size(dz)))/2;
dx = diff(pts(:,1));
CurveArea = abs(dp.*dx);
HorzArea = abs(dy.*dx);
VertArea = abs(dz.*dx);

% Cast Distance and Surge Distance
CastDist = sqrt((diff(pts(:,2)).^2)+(diff(pts(:,3)).^2));       % On the yz plane
SurgeDist = abs(dx);                                                 % On the x plane

% Casting and Surging speed
CastSpeed = CastDist.*sam;
SurgeSpeed = abs(dx).*sam;

end

function [Overall] = FindOverallAverage (Variable,pts,sam)

% This function finds out the overall average of all the global parameters and saves it a a seperate sub-structure.
%
% [Inputs]: A structure with global variables nested in them
%
% [Outputs]: A nested structure of overall averages of the global variables.
%
% Version 1.01
% Last Modified on 28 December 2012
% by Dinesh Natesan

Overall.Distance = mean(Variable.Distance);
Overall.CastDist = mean(Variable.CastDist);
Overall.SurgeDist = mean(Variable.SurgeDist);
Overall.CurveArea = mean(Variable.CurveArea);
Overall.HorzArea = mean(Variable.HorzArea);
Overall.VertArea = mean(Variable.VertArea);
Overall.OdorAxisDist = mean(Variable.OdorAxisDist);

% Distance versus time!
Overall.FlySpeed = mean(Variable.FlySpeed);
Overall.CastSpeed = mean(Variable.CastSpeed);
Overall.SurgeSpeed = mean(Variable.SurgeSpeed);

end

function [Overall] = FindOverallSum (Variable,pts,sam)

% This function finds out the overall average of all the global parameters and saves it a a seperate sub-structure.
%
% [Inputs]: A structure with global variables nested in them
%
% [Outputs]: A nested structure of overall sums of the global variables.
%
% Version 1.01
% Last Modified on 28 December 2012
% by Dinesh Natesan

%Precalculate contants
len = size(pts,1);
time = (len-1)/sam;         %First point is 0

Overall.Distance = sum(Variable.Distance);
Overall.CastDist = sum(Variable.CastDist);
Overall.SurgeDist = sum(Variable.SurgeDist);
Overall.CurveArea = sum(Variable.CurveArea);
Overall.HorzArea = sum(Variable.HorzArea);
Overall.VertArea = sum(Variable.VertArea);
Overall.OdorAxisDist = sum(Variable.OdorAxisDist);
% Displacement versus time!
Overall.FlySpeed = sqrt(((pts(1,1)-pts(len,1)).^2)+((pts(1,2)-pts(len,2)).^2)+((pts(1,3)-pts(len,3)).^2))/time;
Overall.CastSpeed = sqrt(((pts(1,2)-pts(len,2)).^2)+((pts(1,3)-pts(len,3)).^2))/time;
Overall.SurgeSpeed = abs(pts(1,1)-pts(len,1))/time;
Overall.TotalTime = time;

end

function [angvel,angle] = angular_velocity(pts,sam)
%ANGULAR VELOCITY Summary of this function goes here
%function [angvel] = angular_velocity(pts,sam)
%
% Inputs:
% 1) XYZ points of the trajectory (pts)
% 2) Sampling Rate
%
% Outputs: Angular velocity
%
% This function finds the angular velocity and plots it
% Also, it downsamples it
%
% Version 2.01: By Dinesh Natesan
% 18 May 2012
%
%

len = size(pts,1);
% Angle calculation between points
angle = zeros(len-1,1);

%Find angle between two points and divide it by dt
for i = 2:(len-1)
    u = pts(i,:)-pts(i-1,:);
    v = pts(i+1,:)-pts(i,:);
    %  DOT product by magnitude, by dt
%     angle(i) = acos(dot(u,v)/(norm(u)*norm(v)))*180/pi;    % Double checked. Right formula. Turn angle definition exactly same as Breugel et al.
    angle(i) = (atan2(norm(cross(u,v)),dot(u,v)))*180/pi;    % Double checked. Right formula. Turn angle definition exactly same as Breugel et al.
end
angvel = (angle)*sam;
end


function [binnedData] =  StartBinning (binnerInputs,vartobeBinned,binValues,binLimits,sam)

% This function recursively bins all the 'vartobebinned' structure based on the parameters in the 'binnerInputs' structure. It invokes two types of
% binners: 'DistBinner' which bins the variables based on source and odor axis distance and 'TimeBinner' which bins the variable based on time.
% This function outputs the binned data in a neatly arranged structure.
%
% [Inputs]
% binnerInputs -> a structure containing the inputs to the Distbinner.
% vartobeBinned -> a structure containing all the variables to be binned.
% binValues -> the structure contaning the bin size of each type of bin.
% binLimits -> the matrix bin limit of each type of bin.
% sam -> sampling rate of the videos
%
% [Outputs]
% binnedData -> a nested structure containing all the binned data
%
% Version 1.0
% Last Modified 12th May 2012
% by Dinesh Natesan

% Gets all names of the structure variables
binnerNames = fieldnames(binnerInputs);
varNames = fieldnames(vartobeBinned);
binNames = fieldnames(binValues);

% Recursively calls binners to bin all the variables in 'vartobeBinned'
for i=1:size(binNames,1)
    binValue = binValues.(char(binNames(i)));
    if (strcmp(char(binNames(i)),'TimeBin'))
        for j=1:size(varNames,1)
            tobeBinned = vartobeBinned.(char(varNames(j)));
            [binnedData.(char((binNames(i)))).(char((varNames(j)))).zoneMean,binnedData.(char((binNames(i)))).(char((varNames(j)))).zoneStdev,binnedData.(char((binNames(i)))).(char((varNames(j)))).zone] = TimeBinner(tobeBinned,binValue,sam);
        end
    else
        binnerInput = binnerInputs.(char(binnerNames(i)));
        binLimit = binLimits(i);
        for j=1:size(varNames,1)
            tobeBinned = vartobeBinned.(char(varNames(j)));
            [binnedData.(char((binNames(i)))).(char((varNames(j)))).zoneMean,binnedData.(char((binNames(i)))).(char((varNames(j)))).zoneStdev,binnedData.(char((binNames(i)))).(char((varNames(j)))).zoneSum,binnedData.(char((binNames(i)))).(char((varNames(j)))).zone] = DistBinner(tobeBinned,binnerInput,binValue,binLimit);
        end
%         [binnedData.(char((binNames(i)))).(char((binnerNames(i)))).zoneMean,binnedData.(char((binNames(i)))).(char((binnerNames(i)))).zoneStdev,binnedData.(char((binNames(i)))).(char((binnerNames(i)))).zoneSum,binnedData.(char((binNames(i)))).(char((binnerNames(i)))).zone] = DistBinner(binnerInput,[0;binnerInput],binValue,binLimit);
    end
end
end

function [zoneMean,zoneStdev,zoneSum,zone] = DistBinner(tobeBinned,binnerInput,binValue,binLimit)
%This function bins the 'tobebinned' variable based on the binner input, the bin value and the bin limit.
% [zoneMean,zoneStdev,zoneSum,zone] = DistBinner(tobeBinned,binnerInput,binValue,binLimit)
% [Inputs]: Self-explanatory
%
% [Outputs]: zoneMean, zoneStdev and the binned zone cell.
%
% Version 1.2
% Last Modified on 25th May 2012
% by Dinesh Natesan

% Find out the size of the matrix from the bin size and bin limit.
matrixSize = ceil(binLimit/binValue);
% Initialize a cell for the binned variables
zone = cell(matrixSize,2);

% Sorts out the 'tobeBinned' variable into bins
% Since all variables have a size = 1-(binnerInput size), the last point(landing) is ignored
% for j=1:(size(binnerInput)-1)
if length(binnerInput) == length(tobeBinned)
    for j=1:(size(binnerInput))
        if (binnerInput(j)>binLimit)
            % Jump if the binnerinput value is more that the limit
        else
            % Sort based on the zone
            %         index = ceil(binnerInput(j)/binValue);
            index = floor(binnerInput(j)/binValue) + 1;
            zone{index,1}=[zone{index,1};tobeBinned(j)];
            zone{index,2}=[zone{index,2};index];
        end
    end
elseif length(binnerInput) == (length(tobeBinned)+1)
    for j=2:(size(binnerInput))
        if (binnerInput(j)>binLimit)
            % Jump if the binnerinput value is more that the limit
        else
            % Sort based on the zone
            %         index = ceil(binnerInput(j)/binValue);
            index = floor(binnerInput(j)/binValue) + 1;
            if isnan(index) == 1
               flag = 1; 
            end
            zone{index,1}=[zone{index,1};tobeBinned(j-1)];
            zone{index,2}=[zone{index,2};index];
        end
    end
else
    disp('Something wrong. Have not bothered to add error here, since I am confident this should never be displayed.');
    disp('If you are reading this, run!! Either the world is coming to an end or I have made a terrible blunder.');
    disp('Mostly the former though. PS: Called from DistBinner function ');
end

% Finds the mean and the standard deviation of bins seperately
% Initialize stuff
zoneMean=NaN(matrixSize,1);
zoneStdev=NaN(matrixSize,1);
zoneSum = NaN(matrixSize,1);
% Find Mean, Standard Deviation and sum
for i = 1:matrixSize
    zoneMean(i,1) = mean(abs(zone{i,1}));
    zoneStdev(i,1) = std(abs(zone{i,1}));
    zoneSum(i,1) = sum(abs(zone{i,1}));
end

end

function [zoneMean,zoneStdev,zone] = TimeBinner (tobeBinned,binValue,sam)

%This function bins the 'tobebinned' variable into time bins of value 'binValue'. The 't=0' a.k.a begining of time is the time of landing. Hence the
%time binning happens in reverse. i.e bin 1 is time from 0 to binValue.
%
% [Inputs]: Self-explanatory
%
% [Outputs]: zoneMean, zoneStdev and the binned zone cell.
%
% Version 1.1
% Last Modified on 24th May 2012
% by Dinesh Natesan

%Get required variables
binSize = (binValue*sam);
dataSize = size(tobeBinned,1);
matrixSize = ceil(dataSize/binSize);
% Initialize a cell for the binned variables
zone = cell(matrixSize,2);

% Invert and bin from t=0
for i = 1:matrixSize
    max = dataSize-((i-1)*binSize);
    min = dataSize-((i*binSize)-1);
    if(min <=0)
        min = 1;
    end
    zone{i,1} = tobeBinned(max:-1:min,1);
    zone{i,2} = ones((max-min+1),1).*i;
end

% Finds the mean and the standard deviation of bins seperately
zoneMean=NaN(matrixSize,1);
zoneStdev=NaN(matrixSize,1);
for i = 1:matrixSize
    zoneMean(i,1) = mean(abs(zone{i,1}));
    zoneStdev(i,1) = std(abs(zone{i,1}));
    %     zoneSum(i,1) = sum(abs(zone{i,1}));
end

end

function [zone,zoneMean,avgtort] = Tortuosity (Distance,pts,Timebin,sam)

%This function finds the tortuosity of a part of the trajectory in a timebin. Also finds the overall tortuosity of the odor tracking trajectory.
%
% [Inputs]: Self-explanatory
%
% [Outputs]:
% Zone -> Tortuosity in time bins
% avgtort -> Overall tortuosity of the trajectory
%
% Version 1.0
% Last Modified on 14th May 2012
% by Dinesh Natesan

%Get required variables
binSize = (Timebin*sam);
dataSize = size(Distance,1);
matrixSize = ceil(dataSize/binSize);
len = size(pts,1);

% Initialize a cell for the binned variables
zone = cell(matrixSize,2);

% Invert and bin from t=0
for i = 1:matrixSize
    max = dataSize-((i-1)*binSize);
    min = dataSize-((i*binSize)-1);
    if(min <=0)
        min = 1;
    end
    displacement = sqrt((pts(max+1,1)-pts(min,1)).^2+(pts(max+1,2)-pts(min,2)).^2+(pts(max+1,3)-pts(min,3)).^2);
    zone{i,1} = sum(Distance(max:-1:min,1))/displacement;
    zone{i,2} = i;
    % Could have done the same thing from the time binned Distance change.. Both are same, logically and number of steps wise.
end
%Zone Mean
zoneMean = cell2mat(zone(:,1));
% Tortuosity is only defined for time. Essentially, total distance
% travelled in a timebin by displacement.

% Find the overall tortuosity
displacement = sqrt((pts(len,1)-pts(1,1)).^2+(pts(len,2)-pts(1,2)).^2+(pts(len,3)-pts(1,3)).^2);
avgtort = sum(Distance)/displacement;

% Checked!! Working!!
end

function [filt_signal] = ButterFilt( signal,sampling_freq,cutoff_freq)

% Filters in input signal at the given cutoff frequency. The cutoff freq should be atleast less than sam/2 (Nyquist).
%
% [Inputs]:
% Signal (an be a matrix)
% Sampling frequency
% Cutoff frequency
%
% [Outputs]:
% Filtered Signal
%
% Version 2.001
% Last Modified on 24th May 2012
% Last Checked 1st June 2012
% by Dinesh Natesan

if (size(signal,2)>1 || size(signal,3)>1)
    filt_signal=zeros(size(signal,1),size(signal,2));
    for j = 1:size(signal,3)
    for i=1:size(signal,2)
        filt_signal(:,i,j) = ButterFilt( signal(:,i,j),sampling_freq,cutoff_freq);
    end
    end
else
    Ny_freq = sampling_freq/2;
    [b,a]= butter(4,cutoff_freq/Ny_freq,'low');
    %The order of this above designed filter is '4'. This was found out by comparing the magnitude and phase response of each order of butter filter, and
    %selecting a order which minimizes ripple effect and maximizes the passband
    filt_signal= filtfilt(b,a,signal);
end

end