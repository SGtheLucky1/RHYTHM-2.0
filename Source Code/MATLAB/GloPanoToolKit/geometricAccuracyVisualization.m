%% Texture projection %%
 close all
 clear all
 clc
% Load geometry
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

% Visualize
f = figure;
t = trisurf(cells,pts(:,1),pts(:,2),pts(:,3));
t.FaceColor = 'c';
t.EdgeColor = 'none';
axis equal
hold on
blockXY = 12.7;
blockZ = 25.4;
blockVerts = [-blockXY,blockXY,blockZ;
    -blockXY,-blockXY,blockZ;
    blockXY,blockXY,blockZ;
    blockXY,-blockXY,blockZ;
    -blockXY+1,blockXY,-blockZ;
    -blockXY+1,-blockXY,-blockZ;
    blockXY+1,blockXY,-blockZ;
    blockXY+1,-blockXY,-blockZ];
blockEdges = [1 2; 1 3; 1 5;
    2 4; 2 6;
    3 4; 3 7;
    4 8;
    5 6; 5 7;
    6 8;
    7 8];
for n = 1:size(blockEdges,1)
   plot3([blockVerts(blockEdges(n,1),1) blockVerts(blockEdges(n,2),1)],...
       [blockVerts(blockEdges(n,1),2) blockVerts(blockEdges(n,2),2)],...
       [blockVerts(blockEdges(n,1),3) blockVerts(blockEdges(n,2),3)],...
       'r','LineWidth',4)
end
view(-225,16)
light
% light('Position',[20,20,0])

movname = 'geometryValidation';
vidObj = VideoWriter([geoPath movname],'MPEG-4');
set(vidObj,'FrameRate',20)
open(vidObj)
for n = 1:360
view([n 16])
F = getframe(f);
writeVideo(vidObj,F);
end
close(f)
close(vidObj)


% Load texture camera calibrations
cd('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/2017_0113_Rabbit/Geometry/Calibration')
camLabels = 'G';
camAng = 45;
nviews = 4;

% Grab camera parameters for each CMOS camera
pospar=load(sprintf('cal%s_%0.3d.pospar',camLabels,camAng));
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


% Grab camera view directions
lookmap=(Rcmap(1:4,1:4)\[0 0 1 1]')';       % camera view direction
lookmap=lookmap-xyzmap;            % camera view direction
looknmap=lookmap./norm(lookmap);     % normalized view direction


cd('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/2017_0113_Rabbit/Geometry/Cube')
% % % cameraMapping = 'brainvision_ultimaL';
geoCam = 'iDS_UI_3220CP-M-GL_with_f1.2';
a = imread('CubeG.tiff');
n = 5;
[Xigg,Yigg]=pred(centroids,Parmap,Posmap,geoCam);
figure
image(a)
hold on
s = scatter(Xigg,Yigg,'MarkerEdgeColor',[0 1 1]);
edgeColor = s.MarkerEdgeColor;
s.MarkerEdgeColor = [edgeColor 255*0.5];

hold off


% % % 
% % % % Load texture data
% % % scrn_size = get(0,'ScreenSize');
% % % maskMapping = zeros(200,200,4);
% % % maskGeo = zeros(200,200,4);
% % % mask = zeros(200,200,4);
% % % Xigg = cell(4,1);
% % % Yigg = cell(4,1);
% % % cameraMapping = 'brainvision_ultimaL';
% % % for n = 1:4
% % %     tmp = sprintf('image%s = imread(''20150527Rabbit_014_%s.bmp'');',...
% % %         camLabels(n),camLabels(n));
% % %     eval(tmp)
% % %     % Create figure
% % %     tmp = sprintf('f%d = figure(''Position'',[scrn_size(3)/2-400,scrn_size(4)/2-400,800,800]);',n);
% % %     eval(tmp)
% % %     % Background image
% % %     tmp = sprintf('image(image%s,''AlphaData'',0.5)',camLabels(n));
% % %     eval(tmp)
% % %     axis image
% % %     maskMapping(:,:,n) = roipoly;
% % %     tmp = sprintf('close(f%d)',n);
% % %     eval(tmp)
% % %     % Projected data points
% % %     tmp = sprintf('f%d = figure(''Position'',[scrn_size(3)/2-400,scrn_size(4)/2-400,800,800]);',n);
% % %     eval(tmp)
% % %     [Xigg{n},Yigg{n}]=pred(centroids,Parmap(:,n),Posmap(:,n),cameraMapping);
% % %     % Round values so they match pixels and create mask
% % %     Xigg{n} = round(Xigg{n});
% % %     rmX = (Xigg{n} > 0).*(Xigg{n} < 200);
% % %     rmX = unique((rmX == 0).*(1:length(Xigg{n}))');
% % %     if rmX(1) == 0
% % %         rmX(1) = [];
% % %     end
% % %     Yigg{n} = round(Yigg{n});
% % %     rmY = (Yigg{n} > 0).*(Yigg{n} < 200);
% % %     rmY = unique((rmY == 0).*(1:length(Yigg{n}))');
% % %     if rmY(1) == 0
% % %         rmY(1) = [];
% % %     end
% % %     rm = [rmX; rmY];
% % %     rm = unique(rm);
% % %     if rm ~= 0 
% % %         Xigg{n}(rm) = [];
% % %         Yigg{n}(rm) = [];
% % %     end
% % %     % Convert to index
% % %     ind = sub2ind([200 200],Yigg{n},Xigg{n});
% % %     ind = unique(ind);
% % %     % Create mask based on projected centroids
% % %     tmp = zeros(200,200);
% % %     tmp(ind) = 1;
% % %     maskGeo(:,:,n) = tmp;
% % %     maskGeo(:,:,n) = imfill(maskGeo(:,:,n));
% % %     % Combine texture and geometry masks
% % %     mask(:,:,n) = maskMapping(:,:,n).*maskGeo(:,:,n);
% % %     
% % %     % Plot
% % %     scatter(Xigg{n},Yigg{n})
% % %     hold on
% % %     % Background image
% % %     tmp = sprintf('image(image%s,''AlphaData'',0.5)',camLabels(n));
% % %     eval(tmp)
% % %     axis image    
% % %     hold off
% % % end
% % % % return to Goemertry Directory
% % % cd ..
% % % 
% % % % Project texture onto geometry
% % % cmosData = cell(4,1);
% % % cmap = colormap('gray');
% % % for n = 1:4
% % %     tmp = sprintf('cmosData{%d} = rgb2ind(image%s,cmap);',n,camLabels(n));
% % %     eval(tmp)
% % % end
% % % nr = size(cmosData{1},1);
% % % nc = size(cmosData{1},2);
% % % 
% % % [dataProj] = textureProjection(mask,norms,centroids,...
% % %     cells,pts,neighs,neighnum,Parmap,Posmap,cameraMapping,nr,nc,looknmap,...
% % %     cmosData,geoName);
% % % 
% % % % Visualize
% % % figure;
% % % trisurf(cells,pts(:,1),pts(:,2),pts(:,3),dataProj,'LineStyle','none')
% % % axis equal
% % % hold on
% % % blockXY = 12.7;
% % % blockZ = 25.4;
% % % blockVerts = [-blockXY,blockXY,blockZ;
% % %     -blockXY,-blockXY,blockZ;
% % %     blockXY,blockXY,blockZ;
% % %     blockXY,-blockXY,blockZ;
% % %     -blockXY,blockXY,-blockZ;
% % %     -blockXY,-blockXY,-blockZ;
% % %     blockXY,blockXY,-blockZ;
% % %     blockXY,-blockXY,-blockZ];
% % % blockEdges = [1 2; 1 3; 1 5;
% % %     2 4; 2 6;
% % %     3 4; 3 7;
% % %     4 8;
% % %     5 6; 5 7;
% % %     6 8;
% % %     7 8];
% % % for n = 1:size(blockEdges,1)
% % %    plot3([blockVerts(blockEdges(n,1),1) blockVerts(blockEdges(n,2),1)],...
% % %        [blockVerts(blockEdges(n,1),2) blockVerts(blockEdges(n,2),2)],...
% % %        [blockVerts(blockEdges(n,1),3) blockVerts(blockEdges(n,2),3)],'r')
% % % end
% % % light
% % % light('Position',[-20 -20 0])
% % % colormap('gray')
% % % 
% % % 
% % % 
% % % 
