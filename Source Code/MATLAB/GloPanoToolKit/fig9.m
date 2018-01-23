% Script for creating the new figure 8
handles = cell(1,3);
% Load the rat data
load('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/2017_0420_Rat/Optical/07_Rat_0420.mat')
handles{2} = a;
clear a
% Load rabbit data
load('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/2017_0113_Rabbit/Optical/06_projectedData')
handles{1} = a;
clear a

%Preallocate variables for population
maxX = zeros(1,3);
minX = zeros(1,3);
maxY = zeros(1,3);
minY = zeros(1,3);
maxZ = zeros(1,3);
minZ = zeros(1,3);

for n = 1:2
    % Check min and max in X direction
    maxX(n) = max(handles{n}.pts(:,1));
    minX(n) = min(handles{n}.pts(:,1));
    % Check min and max in X direction
    maxY(n) = max(handles{n}.pts(:,2));
    minY(n) = min(handles{n}.pts(:,2));
    % Check min and max in X direction
    maxZ(n) = max(handles{n}.pts(:,3));
    minZ(n) = min(handles{n}.pts(:,3));
    
    figure
    tmp = sprintf('t%d = trisurf(handles{n}.cells,handles{n}.pts(:,1),handles{n}.pts(:,2),handles{n}.pts(:,3))',n);
    eval(tmp)
    tmp = sprintf('t%d.FaceColor = ''c'';',n);
    eval(tmp)
    set(gca,'XLim',[min(minX) max(maxX)],'YLim',[min(minY) max(maxY)],'ZLim',[min(minZ) max(maxZ)])
    view(-45,0)
end

