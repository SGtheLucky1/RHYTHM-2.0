close all
clear all
clc

load('octvert_cuboid_try1.mat')

cd 'Calibration'
pospar=load(sprintf('cal%s_%0.3d.pospar','Geo',45));
% Compute variables need for rotation matrix
wa=pospar(4)*pi/180;  % omega, Wx rotation (world basis, x axis)
pa=pospar(5)*pi/180;  % psi, y rotation  (world basis, y axis)
ra=pospar(6)*pi/180;  % kappa, z rotation (world basis, z axis)
cw=cos(wa); sw=sin(wa);
cp=cos(pa); sp=sin(pa);
cr=cos(ra); sr=sin(ra);
R=zeros(3,3);

% Populate R with rotation matrix from the CMOS cameras
R(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]';
R(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
R(:,3)=[-sp cp*sw cp*cw]';
% Combine R with translation matrix to create tranformation matrix
Rcmap=[[R'; [0 0 0]] [pospar(1) pospar(2) pospar(3) 1]'];

% Compute xyz positions of orbit from calibration matricies
xyzmap=(Rcmap\[0 0 0 1]')';
Parmap=pospar(7:14);
Posmap=pospar(1:6);
cd ..

% Grab camera view directions
lookmap=ones(size(Rcmap,3),4);
looknmap=ones(size(Rcmap,3),3);
for i=1:size(Rcmap,3)
    lookmap(i,:)=(Rcmap(1:4,1:4,i)\[0 0 1 1]')';       % camera view direction
    lookmap(i,1:3)=lookmap(i,1:3)-xyzmap(i,1:3);            % camera view direction
    looknmap(i,:)=lookmap(i,1:3)./norm(lookmap(i,1:3));     % normalized view direction
end

% GEOMETRY
% Have user select geometry file to use

[geoName,geoPath] = uigetfile('*.vtk','Select the VTK output');
geoName = geoName(1:end-8);
% Load associated geometric into handles variables
[centroids,~,~]=readdat([geoPath geoName 'centroids.dat']);
[norms,~,~]=readdat([geoPath geoName 'normals.dat']);
[neighnum,neighs,~,~]=readneighs([geoPath geoName 'neighs.dat']);   % 2/23/05 MWKay
% neighs is a structure
[pts,numparams,txtparams]=readdat([geoPath geoName 'pts.dat']);
[cells,~,~]=readdat([geoPath geoName 'cells.dat']);
cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
for i=1:size(cells,1)            % VTK starts with zero, matlab starts with one
    for j=1:neighnum(i)
        neighs{i}(j)=neighs{i}(j)+1;
    end
end

% FIGURE OUT HOW TO GET THIS INFORMATION FROM FIT2 TO HERE !!!!
mapCam = 'brainvision_ultimaL';
geoCam = 'iDS_UI_3220CP-M-GL_with_f1.2';

% Load in silhouette images
cd 'Heart'
numImages = 72;
silhs = cell(numImages,1);
for n = 1:numImages
    silhs{n} = eval(sprintf('imread(''Cube%0.2d.tiff'')',n));
end
cd ..


% Project geometry onto silhouette images
% Map 3D points back to Geo camera
figure
[X,Y] = pred(pts,...
    Parmap,Posmap,geoCam);
subplot(1,2,1)
image(silhs{1})
hold on
scatter(X,Y,'ro')
hold off

subplot(1,2,2)
image(silhs{1})
hold on
[Vx,Vy] = pred(vert(vert(:,4)==1,1:3),Parmap,Posmap,geoCam);
scatter(Vx,Vy,'ro')
[Nx,Ny] = pred(vert(vert(:,5)==5,1:3),Parmap,Posmap,geoCam);
scatter(Nx,Ny,'go')
hold off



% Create Z-axis rotation matrix
Rz = [cosd(-20) -sind(-20) 0;
    sind(-20) cosd(-20) 0;
    0 0 1];

figure
[X,Y] = pred(pts,...
    Parmap,Posmap,geoCam);
subplot(3,6,1)
image(silhs{1})
hold on
scatter(X,Y,'ro')
hold off

% Apply rotation to points
rotPts = (Rz*pts')';

for n = 5:4:72
[X,Y] = pred(rotPts,Parmap,Posmap,geoCam);
subplot(3,6,(n+3)/4)
image(silhs{n})
hold on
scatter(X,Y,'ro')
rotPts = (Rz*rotPts')';
end


lvl5 = vert(:,5) == 6;
lvl5keep = ((vert(:,5)==6).*(vert(:,4) == 1))==1;
% Cuboid dimensions
unit = 25.4;
unitHalf = unit/2;
% Cuboid corner coordinates
corners = [-unitHalf,unitHalf,unit;
    unitHalf,unitHalf,unit;
    -unitHalf,-unitHalf,unit;
    unitHalf,-unitHalf,unit;
    -unitHalf,unitHalf,-unit;
    unitHalf,unitHalf,-unit;
    -unitHalf,-unitHalf,-unit;
    unitHalf,-unitHalf,-unit];
% Corner Order (CO)
CO= [1,1,1,4,4,4,7,7,7,6,6,6;3,2,5,3,2,8,5,8,3,5,8,2]';
figure
hold on
for n = 1:12
plot3([corners(CO(n,1),1) corners(CO(n,2),1)],...
    [corners(CO(n,1),2) corners(CO(n,2),2)],...
    [corners(CO(n,1),3) corners(CO(n,2),3)],'c')
end
scatter3(vert(lvl5,1),vert(lvl5,2),vert(lvl5,3),'go','filled')
scatter3(vert(lvl5keep,1),vert(lvl5keep,2),vert(lvl5keep,3),'ro')
hold off



