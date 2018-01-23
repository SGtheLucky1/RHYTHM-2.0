% % % %% Optical resolution calculation %%
% % % 
% % % disp('Navigate to the geometry directory.')
% % % geoDir = uigetdir;
% % % cd(geoDir)

close all
clear all
clc

% CAMERA RESOLUTIONS
[camRes,camResMicro] = resolutionCalc(1);

% CALIBRATION
nviews = 1:4;
[Rcmap,xyzmap,par,pos,lookmap,looknmap] = loadCameraCalibration(nviews);

% GEOMETRY
% Have user select geometry file to use
[pts,numparams,txtparams,cells,centroids,norms,neighs,neighnum] = loadVTKmesh();

% Calculate inter centroidal resolution
mapCam = 'brainvision_ultimaL';
geoCam = 'iDS_UI_3220CP-M-GL_with_f1.2';

newDistAve = zeros(4,2);
newDistStd = zeros(4,2);
steps = [0.5 0.75 0.9 0.95 0.99 0.995];
nf = cell(2,1);
nf{1} = cells;
nv = cell(2,1);
nv{1} = pts;
newCent = cell(2,1);
newCent{1} = centroids;
newNorms = cell(2,1);
newNorms{1} = norms;
res = 1;
cnt = 1;
rec = [cnt steps(res(1))];
while res < length(steps)+1
    % Decrement step of simplification if previous step overshot
    if mean(newDistAve(:,1)) > mean(camRes(1:4))
        % Reset the values to the previous iteration to try again
        newDistAve(:,1) = newDistAve(:,2);
        newDistStd(:,1) = newDistStd(:,2);
        nf{1} = nf{2};
        nv{1} = nv{2};
        newCent{1} = newCent{2};
        newNorms{1} = newNorms{2};
        % Decrease the percentage
        res = res+1;
    else
        % Place current iteration in second slot
        newDistAve(:,2) = newDistAve(:,1);
        newDistStd(:,2) = newDistStd(:,1);
        nf{2} = nf{1};
        nv{2} = nv{1};
        newCent{2} = newCent{1};
        newNorms{2} = newNorms{1};
        
        t = trisurf(nf{1},nv{1}(:,1),nv{1}(:,2),nv{1}(:,3));
        t.FaceColor = 'c';
        
        distAve = zeros(4,1);
        distStd = zeros(4,1);
        for n = 1:4
            [distAve(n),distStd(n)] = interCentroidResolution(newCent{1},newNorms{1},par(:,n),...
                pos(:,n),mapCam);
        end
        
        [nf{1},nv{1}] = reducepatch(t,steps(res));
        
        newCent{1} = zeros(size(nf{1},1),3);
        % Calculate new centroids at these faces
        newCent{1}(:,1) = sum(reshape(nv{1}(nf{1},1),[size(nf{1},1) 3]),2)/3;
        newCent{1}(:,2) = sum(reshape(nv{1}(nf{1},2),[size(nf{1},1) 3]),2)/3;
        newCent{1}(:,3) = sum(reshape(nv{1}(nf{1},3),[size(nf{1},1) 3]),2)/3;
        % Calculate new normals at these centroids
        v1 = nv{1}(nf{1}(:,1),:)-newCent{1};
        v2 = nv{1}(nf{1}(:,2),:)-newCent{1};
        newNorms{1} = zeros(size(v1,1),3);
        for n = 1:size(v1,1)
            newNorms{1}(n,:) = cross(v2(n,:),v1(n,:));
            if sum(-newCent{1}(n,:).*newNorms{1}(n,:),2) > 0
                newNorms{1}(n,:) = cross(v1(n,:),v2(n,:));
            end
        end
        
        % % %     newDistAve = zeros(4,1);
        % % %     newDistStd = zeros(4,1);
        for n = 1:4
            [newDistAve(n,1),newDistStd(n,1)] = interCentroidResolution(newCent{1},newNorms{1},par(:,n),...
                pos(:,n),mapCam);
        end
        
        cnt = cnt + 1;
        rec = [rec; cnt steps(res)];
    end
end



% % % current_dir = pwd;
% % % cd('/Users/Chris/Documents/MATLAB/GloPanoToolkit/smoothpatch_version1b/')
% % % mex smoothpatch_curvature_double.c -v
% % % mex smoothpatch_inversedistance_double.c -v
% % % mex vertex_neighbours_double.c -v

FV.faces = cells;
FV.vertices = pts;

FV2=smoothpatch(FV,1,100);
FV3=smoothpatch(FV,0,100);

figure, 
subplot(1,3,1), patch(FV,'FaceColor',[1 0 0],'EdgeAlpha',0);  view(3); camlight
subplot(1,3,2), patch(FV2,'FaceColor',[0 1 0],'EdgeAlpha',0); view(3); camlight
subplot(1,3,3), patch(FV3,'FaceColor',[0 0 1],'EdgeAlpha',0); view(3); camlight

