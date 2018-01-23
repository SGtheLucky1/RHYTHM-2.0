% % % clear all
% % % close all
load('05_processedProjected')

xyz = handles.centroids;
tri = handles.pts;
cell = handles.cells;
% val = a.actMapGeo;
val = apdMapGeo;
pot = handles.dataProj;
norm = handles.norms;

% Create a z-axis rotation matrix for a 180 degree rotation
Rz = [cos(pi) -sin(pi) 0;
    sin(pi) cos(pi) 0;
    0 0 1];
% Apply to centroids and vertices
xyz = (Rz*xyz')';
tri = (Rz*tri')';

% Remove the NAN values
rm = find(val ~= 0);
skel = cell(val==0,:);
val = val(rm);
pot = pot(rm,:);
xyz = xyz(rm,:);
cell = cell(rm,:);
norm = norm(rm,:);


% Remove all cells not part of the largest connected component
rm = removeIslands(cell);
skel = [skel; cell(rm(:,2)~=1,:)];
rm = rm(rm(:,2) == 1,1);
val = val(rm);
pot = pot(rm,:);
xyz = xyz(rm,:);
cell = cell(rm,:);
norm = norm(rm,:);
maxCoord = max(max(xyz));
xyz = xyz/maxCoord;

maxPts = max(max(tri));
tri = tri/maxPts;

% Convert to spherical from cartesian
r = sqrt(tri(:,1).^2+tri(:,2).^2+tri(:,3).^2);
phi = atan2(tri(:,2),tri(:,1));
theta = -1*(acos(tri(:,3)./r)-pi/2);

% Convert to spherical from cartesian
xyzR = sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2);
xyzP = atan2(xyz(:,2),xyz(:,1));
xyzT = -1*(acos(xyz(:,3)./xyzR)-pi/2);

% Identify cells that have vertices that cross the map
cellX = phi(cell);
cellXDist = cellX(:,2)-cellX(:,1);
cellXDist(:,2) = cellX(:,3)-cellX(:,1);
cellXDist(:,3) = cellX(:,2)-cellX(:,3);
cellXRemove = abs(cellXDist) > 4;
cellXRemove = cellXRemove.*repmat((1:size(cellXRemove,1))',[1 3]);
cellXRemove = unique(cellXRemove);
cellXRemove = cellXRemove(2:end);

cellM = cell;
cellM(cellXRemove,:) = [];
valM = val;
valM(cellXRemove,:) = [];
% % % cellThetaDistAve = mean(reshape(cellThetaDist,[size(cellThetaDist,1)*3 1]));
% % % cellThetaDistStd = std(reshape(cellThetaDist,[size(cellThetaDist,1)*3 1]));


%% MERCATOR PROJECTION %%
figure
p = patch('Faces',cellM,'Vertices',[phi theta]);
set(gca,'CLim',[0 60])
set(p,'FaceColor','flat',...
'FaceVertexCData',round(valM),...
'EdgeColor','none',...
'CDataMapping','scaled')
axis equal
colormap('jet')
set(gca,'TickDir','out','XLim',[-pi,pi],'XTick',[-pi -pi/2 0 pi/2 pi],...
    'YLim',[-pi/2,pi/2],'YTick',[-pi/2 0 pi/2])

% Create an interpolated grid
[gridP,gridT] = meshgrid(-pi:pi/72:pi,-pi/2:pi/72:pi/2);
V = griddata(xyzP,xyzT,val,gridP,gridT);
gridR = griddata(xyzP,xyzT,xyzR,gridP,gridT);

% Convert back to cartesian
intX = gridR.*sin(gridT+pi/2).*cos(gridP+pi);
intY = gridR.*sin(gridT+pi/2).*sin(gridP+pi);
intZ = -1*gridR.*cos(gridT+pi/2);

% % % figure
% % % surf(gridP,gridT,zeros(size(gridP,1),size(gridP,2)),V,'EdgeColor','none');
% % % % % % surf(X,Y,R,V,'EdgeColor','none');
% % % % % % surf(intX,intY,intZ,V,'EdgeColor','none')
% % % view(0,90)
% % % axis equal
% % % set(gca,'TickDir','out','XLim',[-pi,pi],'XTick',[-pi -pi/2 0 pi/2 pi],...
% % %     'YLim',[-pi/2,pi/2],'YTick',[-pi/2 0 pi/2])
% % % colormap('jet')
% % % caxis([0 60])

%% HAMMER PROJECTION %%
[hammerX,hammerY]=pr_hammer(gridP,gridT,gridR);
figure
surf(hammerX,hammerY,zeros(size(hammerX,1),size(hammerX,2)),V,'EdgeColor','none')
view(0,90)
colormap('jet')
axis equal
caxis([0 60])
hold on
for n = 1:4*3:size(intX,2)
    plot(hammerX(:,n),hammerY(:,n),'k')
end
for m = 1:4*3:size(intX,1)
   plot(hammerX(m,:),hammerY(m,:),'k')
end

%% TIME COURSE OF ACTIVATION WITH MEMBRANE POTENTIAL IMAGES %%

% Rotate the heart back to its orginal position
Rz = [cos(-pi) -sin(-pi) 0;
    sin(-pi) cos(-pi) 0;
    0 0 1];
% Apply to centroids and vertices
xyz = (Rz*xyz')';
tri = (Rz*tri')';

% Time values to capture
index = [395 404 415 427 434];
time = [11 20 31 43 50];
ang = -45;
rot = [0 0;
    45 -45;
    90 -90;
    135 -135;
    180 -180];
rotV = rot+ang;
rotH = rotV-90;

camNorm = zeros(size(rotH,1),size(rotH,2),3);
camNorm(:,:,3) = zeros(size(rotH,1),size(rotH,2));
camNorm(:,:,1) = sind(rotH);
camNorm(:,:,2) = cosd(rotH);

for n = 1:length(index)
    
% % %     camFace = dot(repmat(squeeze(camNorm(n,1,:))',[size(norm,1) 1]),norm,2);
% % %     camFaceYes = find(camFace > 0);
% % %     camFaceNo = find(camFace < 0);
% % %     alpha = zeros(size(camFaceYes,1),1);
% % %     alpha(camFaceYes) = 0.2;
% % %     % Membrane potential at varios time points and angles
% % %     figure
% % %     t1 = trisurf(cell(camFaceYes,:),tri(:,1),tri(:,2),tri(:,3),pot(camFaceYes,index(n)),'LineStyle','none','FaceAlpha',0.5);
% % %     hold on
% % %     t2 = trisurf(cell(camFaceNo,:),tri(:,1),tri(:,2),tri(:,3),pot(camFaceNo,index(n)),'LineStyle','none');
% % %     % set vertex transparencies
% % %     view(rotV(n,1),0)
    figure
    trisurf(cell,tri(:,1),tri(:,2),tri(:,3),pot(:,index(n)),'LineStyle','none')
    view(rotV(n,1),0)
    axis equal
    colormap('jet')
    set(gca,'Visible','off')
    axis equal
    hold on
    trisurf(skel,tri(:,1),tri(:,2),tri(:,3),pot(:,index(n)),'FaceColor','none','EdgeColor',[0.5 0.5 0.5])
    hold off
    
    figure
    trisurf(cell,tri(:,1),tri(:,2),tri(:,3),pot(:,index(n)),'LineStyle','none')
    view(rotV(n,2),0)
    axis equal
    colormap('jet')
    set(gca,'Visible','off')
    hold on
    trisurf(skel,tri(:,1),tri(:,2),tri(:,3),pot(:,index(n)),'FaceColor','none','EdgeColor',[0.5 0.5 0.5])
    hold off
end

