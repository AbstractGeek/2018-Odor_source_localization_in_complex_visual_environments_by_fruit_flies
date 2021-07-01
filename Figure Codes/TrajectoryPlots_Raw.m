function [] = TrajectoryPlots_Raw(light_correct, dark_correct,...
    light_incorrect, dark_incorrect)
% function [] = TrajectoryPlots_Raw(light_correct, dark_correct,...
%     light_incorrect, dark_incorrect)
% 
% 
% Dinesh Natesan
% 10th Aug 2016

% Defaults
sam = 100;
cutoff = 30;
figure_size = [60 60 1260 960];
titleName = 0;
fileNames = 0;
if nargin == 0
    light_correct = [0.5843    0.8157    0.9882];       % light blue
    dark_correct = [0.0118    0.2627    0.8745];        % blue
    light_incorrect = [1.0000    0.8196    0.8745];     % light pink
    dark_incorrect = [0.8980         0         0];      % red
end

% Find csv files
rootdirname = pwd;
csvlist = dir('*.csv');
csvlist = {csvlist(:).name}';    

% Get saveName
[~, treatmentName, ~] = fileparts(rootdirname);
saveName = fullfile(rootdirname,strcat(treatmentName,'_trajPlot'));
if ~isdir(saveName)
    mkdir(saveName)
end

% Figure out treatment
if strfind(treatmentName,'tov')
    if ~isempty(strfind(treatmentName,'noair')) || ~isempty(strfind(treatmentName,'stillair'))
        treatment = 'ov_noair';
    else
        treatment = 'ov_air';
    end
elseif strfind(treatmentName,'tvo')
    if ~isempty(strfind(treatmentName,'noair'))||~isempty(strfind(treatmentName,'stillair'))
        treatment = 'vo_noair';
    else
        treatment = 'vo_air';
    end    
elseif strfind(treatmentName,'tonlyvision')
    if ~isempty(strfind(treatmentName,'noair'))||~isempty(strfind(treatmentName,'stillair'))
        treatment = 'onlyvision_noair';
    else
        treatment = 'onlyvision_air';
    end
elseif strfind(treatmentName,'tonlyodor')
    if ~isempty(strfind(treatmentName,'noair'))||~isempty(strfind(treatmentName,'stillair'))
        treatment = 'onlyodor_noair';
    else
        treatment = 'onlyodor_air';
    end
elseif strfind(treatmentName,'tvv')
    if ~isempty(strfind(treatmentName,'noair'))||~isempty(strfind(treatmentName,'stillair'))
        treatment = 'vv_noair';
    else
        treatment = 'vv_air';
    end
else
    if ~isempty(strfind(treatmentName,'noair')) || ~isempty(strfind(treatmentName,'stillair'))
        treatment = '0_noair';
    else
        treatment = '0_air';
    end    
end


% Begin plotting

% Initialize-figure 1 (combined figure)
h1=figure('Units','Pixels','Position',figure_size,'Visible','on');
plot3(0,0,0);
hold on;
grid on;

% Initialize-figure 2 (single figure)
h2=figure('Units','Pixels','Position',figure_size,'Visible','on');

for i=1:length(csvlist)
    % Get filename
    fname = csvlist{i};
    fname = strsplit(fname,'.');
    fname = strjoin(fname(1:end-1),'_');
    % Initialize h2
    set(groot,'CurrentFigure',h2);
    plot3(0,0,0);
    hold on;
    grid on;
    
    [pts,cs,bs,landingPoint] = Datapoints(...
        dlmread(fullfile(rootdirname,csvlist{i})),sam,cutoff);
    
    % Assign color based on landing point   
    if landingPoint==1
        % Correct landing - Blue color
        dark_colors = dark_correct;
        light_colors = light_correct;
    else
        % Incorrect landing - Red color
        dark_colors = dark_incorrect;
        light_colors = light_incorrect;
    end    
    
    % Plot figure 2 first
    set(groot,'CurrentFigure',h2);    
    drawTrajPlot(pts,cs,bs,treatment,dark_colors,1,char(fname),fileNames);   
    
    % Save and close the figure
    if titleName==1
        title(fname);
    end
    
    saveas(h2,fullfile(saveName,strcat(char(fname),'.fig')));
    clf(h2);
    
    % Plot figure 1 now;
    set(groot,'CurrentFigure',h1);
    if i==length(csvlist)
        drawTrajPlot(pts,cs,bs,treatment,light_colors,1,char(fname),fileNames);
    else
        drawTrajPlot(pts,cs,bs,treatment,light_colors,0,char(fname),fileNames);
    end
    
end

saveas(h1,strcat(saveName, '.fig'));
close(h1);
close(h2);
end

function [] = drawTrajPlot(pts,cs,bs,treatment,varargin)

% Initial Treatment
% If its ov0cm - then mention, else make it 0.
% treatment = 'ov0cm';
rcolor = [1 0 1];

if (isempty(varargin))
    colors = repmat([0 0 1],length(pts)-1,1);        
    decorate = 1;
    fileNames = 0;
elseif size(varargin,2)==1
    colors = varargin{1};
    decorate = 1;
    fileNames = 0;
elseif size(varargin,2)==2
    colors = varargin{1};
    decorate = varargin{2};
    fileNames = 0;  
elseif size(varargin,2)==3
    colors = varargin{1};
    decorate = varargin{2};    
    textString =  varargin{3};  
    fileNames = 1;
else
    colors = varargin{1};
    decorate = varargin{2};    
    textString =  varargin{3};  
    fileNames = varargin{4};
end

% Initialize
hold on;
grid on;

% Invert axis, but maintain the right hand rule
pts = [-pts(:,1),pts(:,2),-pts(:,3)];
cs = [-cs(:,1),cs(:,2),-cs(:,3)];
bs = [-bs(:,1),bs(:,2),-bs(:,3)];

% Linear Transform to the sphere of interest
pts = pts-repmat(cs(1,:),size(pts,1),1);
bs=bs-repmat(cs(1,:),size(bs,1),1);
cs=cs-repmat(cs(1,:),size(cs,1),1);

% Plot trajectory
plot3(pts(:,1),pts(:,2),pts(:,3),'color',colors,...
     'LineWidth',2.0,'DisplayName', textString);

% Mark the start points with filled 'x'
plot3(pts(1,1),pts(1,2),pts(1,3),'color',rcolor,'Marker','x',...
    'MarkerFaceColor',rcolor,'LineWidth',2.0,'HandleVisibility','off');
% Mark the end points with filled 'o'
plot3(pts(length(pts),1),pts(length(pts),2),pts(length(pts),3),...
    'color',colors,'Marker','o','MarkerFaceColor',colors,'LineWidth',2.0,...
    'HandleVisibility','off');

if (fileNames==1)
    text(pts(1,1),pts(1,2),pts(1,3),strcat('   ',char(textString)));
end

% Decorate plot, i.e. add spheres and stuff
if decorate==1
    decoratePlot(cs,bs,treatment);
end

end

function [] = decoratePlot(cs,bs,treatment)
%
%
% Break down treatment
C = strsplit(treatment,'_');
treatment = C{1};
airflow = C{2};

% Initialize sphere
[x,y,z]=sphere;
[x1,y1,z1]=cylinder(0.025);
color=zeros(21,21,3);
colorr=color;
colorr(:,:,1)=ones(21,21,1);
cylcolor = zeros(2,21,3);
cylcolorr = cylcolor;
cylcolorr(:,:,1) = ones(2,21);
% cylcolor = 0.60*ones(2,21,3);
row=size(cs,1); 

if (row<3) 
    if (strcmp(treatment,'ov'))&&(row==1)
        % Plot odor,visual cue
        a=(0.3*x)+cs(1,1);
        b=(0.3*y)+cs(1,2);
        c=(0.3*z)+cs(1,3);        
        surf(a,b,c,colorr,'EdgeColor','r');
    elseif strcmp(treatment,'onlyvision')
        % Plot only visual cue
        a=(0.3*x)+cs(1,1);
        b=(0.3*y)+cs(1,2);
        c=(0.3*z)+cs(1,3);        
        surf(a,b,c,color);
    elseif strcmp(treatment,'vo')
        % Plot odor,visual cue
        a=(0.3*x)+cs(1,1);
        b=(0.3*y)+cs(1,2);
        c=(0.3*z)+cs(1,3);        
        surf(a,b,c,colorr,'EdgeColor','r');
        % Plot visual capillary
        a=x1+cs(2,1);
        b=y1+cs(2,2);
        c=z1-0.5+cs(2,3);
        surf(a,b,c,cylcolor,'FaceAlpha',0.75,'EdgeColor','None');
        
    elseif strcmp(treatment,'vv')
        % Plot odor,visual cue
        a=(0.3*x)+cs(1,1);
        b=(0.3*y)+cs(1,2);
        c=(0.3*z)+cs(1,3);        
        surf(a,b,c,colorr,'EdgeColor','r');
        % Plot odor,visual cue
        a=(0.3*x)+cs(2,1);
        b=(0.3*y)+cs(2,2);
        c=(0.3*z)+cs(2,3);        
        surf(a,b,c,color);                
    else
        % Plot odor capillary
        a=x1+cs(1,1);
        b=y1+cs(1,2);
        c=z1-0.5+cs(1,3);
        surf(a,b,c,cylcolorr,'FaceAlpha',0.75,'EdgeColor','r');
        
        if (row > 1)
            % Plot visual cue sphere
            a=(0.3*x)+cs(2,1);
            b=(0.3*y)+cs(2,2);
            c=(0.3*z)+cs(2,3);
            surf(a,b,c,color);
        end
    end
    
    for i = 1:row
        % %  Plot capillaries
        plot3([bs(i,1),cs(i,1)],[bs(i,2),cs(i,2)],[bs(i,3),cs(i,3)],'k');
    end
else
    % Plot odor,visual cue
    a=(0.3*x)+cs(1,1);
    b=(0.3*y)+cs(1,2);
    c=(0.3*z)+cs(1,3);
    surf(a,b,c,colorr,'EdgeColor','r');
    % Plot other cues
    for i = 2:row
        % Plot visual cue sphere
        a=(0.3*x)+cs(i,1);
        b=(0.3*y)+cs(i,2);
        c=(0.3*z)+cs(i,3);
        surf(a,b,c,color);
    end
    
    if row == 3
        plot3([bs(2,1),cs(1,1)],[bs(2,2),cs(1,2)],[bs(2,3),cs(1,3)],'k');
        plot3([bs(1,1),cs(2,1)],[bs(1,2),cs(2,2)],[bs(1,3),cs(2,3)],'k');
        plot3([bs(3,1),cs(3,1)],[bs(3,2),cs(3,2)],[bs(3,3),cs(3,3)],'k');
    elseif row==7
        plot3([bs(5,1),cs(1,1)],[bs(5,2),cs(1,2)],[bs(5,3),cs(1,3)],'k');
        plot3([bs(1,1),cs(2,1)],[bs(1,2),cs(2,2)],[bs(1,3),cs(2,3)],'k');
        plot3([bs(2,1),cs(3,1)],[bs(2,2),cs(3,2)],[bs(2,3),cs(3,3)],'k');
        plot3([bs(3,1),cs(4,1)],[bs(3,2),cs(4,2)],[bs(3,3),cs(4,3)],'k');
        plot3([bs(4,1),cs(5,1)],[bs(4,2),cs(5,2)],[bs(4,3),cs(5,3)],'k');
        plot3([bs(6,1),cs(6,1)],[bs(6,2),cs(6,2)],[bs(6,3),cs(6,3)],'k');
        plot3([bs(7,1),cs(7,1)],[bs(7,2),cs(7,2)],[bs(7,3),cs(7,3)],'k');
    else
        disp('Unknown treatment');
        fprintf('\n%d number of rows. Fix it - Quiting for now\n',row);
    end
    
end


% Set axis
v=axis;
pmax=round(max(v));
pmin=round(min(v));
% Make axis equal and square
v1=[v(1)-(pmax-pmin - (v(2)-v(1)))/2 v(2)+(pmax-pmin - (v(2)-v(1)))/2 ...
    v(3)-(pmax-pmin - (v(4)-v(3)))/2 v(4)+(pmax-pmin - (v(4)-v(3)))/2 ...
    v(5)-(pmax-pmin - (v(6)-v(5)))/2 v(6)+(pmax-pmin - (v(6)-v(5)))/2];
axis on;
axis manual;
axis equal;
axis(v1);
axis vis3d;

if ~strcmp(airflow,'noair') && ~strcmp(treatment,'onlyvision')
    % %  Plot odor axis
    plot3([cs(1,1),pmax],[cs(1,2),cs(1,2)],[cs(1,3),cs(1,3)],'r');
    % % Plot odor cylinder
    [y,z,x]=cylinder(0.5);  % Radius 0.5 cm
    x(2,:) = x(2,:)+pmax-1;
    surf(x,y,z,cylcolorr,'EdgeColor','none','FaceAlpha',0.20);
    theta = (0:0.1:2*pi)';
    y = 0.5.*cos(theta);
    z = 0.5.*sin(theta);
    x = pmax.*ones(size(theta));
    fill3(x,y,z,[1 0 0],'EdgeColor','none','FaceAlpha',0.20);
end

end

function [pts,cs,bs,landingPoint] = Datapoints(hspecdata,sam,cutoff)

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

[~,landingPoint]= min(sqrt((cs(:,1)-pts(end,1)).^2+(cs(:,2)-pts(end,2)).^2+(cs(:,3)-pts(end,3)).^2));

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