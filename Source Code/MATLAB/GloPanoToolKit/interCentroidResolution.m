function [distAve,distStd] = interCentroidResolution(centroids,norms,par,pos,cam)
% % % % GEOMETRY
% % % % Have user select geometry file to use
% % % cd(handles.vtk_dir)
% % % [geoName,geoPath] = uigetfile('*.vtk','Select the VTK output');
% % % geoName = geoName(1:end-8);
% % % % Load associated geometric into handles variables
% % % [handles.centroids,~,~]=readdat([geoPath geoName 'centroids.dat']);
% % % [handles.norms,~,~]=readdat([geoPath geoName 'normals.dat']);
% % % [handles.neighnum,handles.neighs,~,~]=readneighs([geoPath geoName 'neighs.dat']);   % 2/23/05 MWKay
% % % % neighs is a structure
% % % [handles.pts,handles.numparams,handles.txtparams]=readdat([geoPath geoName 'pts.dat']);
% % % [cells,~,~]=readdat([geoPath geoName 'cells.dat']);
% % % handles.cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
% % % for i=1:size(cells,1)            % VTK starts with zero, matlab starts with one
% % %     for j=1:handles.neighnum(i)
% % %         handles.neighs{i}(j)=handles.neighs{i}(j)+1;
% % %     end
% % % end
% % % % Project points onto 2D camera masks
% % % handles.X = cell(4,1);
% % % handles.newX = handles.X;
% % % handles.Xshift = handles.X;
% % % handles.Y = cell(4,1);
% % % handles.newY = handles.Y;
% % % handles.Yshift = handles.Y;
% % % handles.shift = zeros(4,2);
% % % handles.geommasks = zeros(size(handles.cmosData{1},1),...
% % %     size(handles.cmosData{1},2),4);

% FIGURE OUT HOW TO GET THIS INFORMATION FROM FIT2 TO HERE !!!!
% % % handles.mapCam = 'brainvision_ultimaL';
% % % handles.geoCam = 'iDS_UI_3220CP-M-GL_with_f1.2';

% Map 3D points to 2D mapping cameras

% Mapping camera normal
mapNorm = [0 0 1];

% Construct rotation matrix
wa=pos(4)*pi/180;  % omega 
pa=pos(5)*pi/180;  % psi   
ra=pos(6)*pi/180;  % kappa 
cw=cos(wa); sw=sin(wa);
cp=cos(pa); sp=sin(pa);
cr=cos(ra); sr=sin(ra);

Rc=zeros(3,3);
Rc(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]'; 
Rc(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
Rc(:,3)=[-sp cp*sw cp*cw]';

% Next transform from world (calibration) coordinate frame to camera frame
tc=pos(1:3);
mapNorm=(Rc'*mapNorm'+tc)';

% Identify camera facing points by performing the dot product against the
% mapping camera norm, all camera facing results will be negative
test = sum(repmat(mapNorm,[size(norms,1) 1]).*norms,2);
test = test < 0;

% Project geometry onto the mapping camera
[X,Y] = pred(centroids(test,:),par,pos,cam);

distMin = zeros(size(X,1),1);
    for n = 1:size(X,1)
        % Grab the point of interest
        poi = [X(n),Y(n)];
        poi = repmat(poi,[size(X,1)-1 1]);
        % Create a variable of all the points, excluding the current one
        tgts = [X Y];
        tgts(n,:) = [];
        % find the distance from the current point to all others
        dist = sqrt(sum((tgts-poi).^2,2));
        % find the minimum
        distMin(n) = min(dist);
    end
distAve = mean(distMin);
distStd = std(distMin);
end


% For each camera:
% 1) isolate camera facing from away facing
% 2) find the shortest distance to an adjacent centroid for each one
% 3) average and report
