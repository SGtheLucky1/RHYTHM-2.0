%% Code for generating images for figure 8

load('05_processedProjected')

handles = a;
clear a

load('tmp')
handles.apdMapGeo = ac;
clear ac

N = isnan(handles.actMapGeo);
V = N == 0;

% % % p1 = [6.3771 6.3924 2.0166];
% % % p2 = [7.0512 5.9679 -1.1187];
% % % p3 = [7.3006 4.2505 -4.8970];
p1 = [2.1076 7.6758 4.0046];
p2 = [-0.6076 7.2077 -3.1480];
p3 = [-6.6563 0.4348 -9.1168];
x = [floor(p1(1)*10)/10 ceil(p1(1)*10)/10;
    floor(p2(1)*10)/10 ceil(p2(1)*10)/10;
    floor(p3(1)*10)/10 ceil(p3(1)*10)/10];
y = [floor(p1(2)*10)/10 ceil(p1(2)*10)/10;
    floor(p2(2)*10)/10 ceil(p2(2)*10)/10;
    floor(p3(2)*10)/10 ceil(p3(2)*10)/10];
z = [floor(p1(3)*10)/10 ceil(p1(3)*10)/10;
    floor(p2(3)*10)/10 ceil(p2(3)*10)/10;
    floor(p3(3)*10)/10 ceil(p3(3)*10)/10];
% find centroid index for these points
pInd = zeros(3,1);
for n = 1:3
        test = (handles.centroids(:,1) > x(n,1)).*(handles.centroids(:,1) < x(n,2));
        test(:,2) = (handles.centroids(:,2) > y(n,1)).*(handles.centroids(:,2) < y(n,2));
        test(:,3) = (handles.centroids(:,3) > z(n,1)).*(handles.centroids(:,3) < z(n,2));
        pInd(n) = find(sum(test,2) == 3);
end

% % % % Membrane Potential
% % % figure
% % % trisurf(handles.cells(V,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),handles.dataProj(V,31),'LineStyle','none')
% % % colormap('jet')
% % % axis equal
% % % view(-135,0)
% % % caxis([0 1])
% % % hold on
% % % trisurf(handles.cells(N,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),'FaceColor','none','LineWidth',0.1);
% % % set(gca,'Visible','off')

% Signals
figure
for n = 1:3
    subplot(3,1,n)
    plot(handles.dataProj(pInd(n),1:420))
    set(gca,'XLim',[0 420])
end

% Activation
figure
trisurf(handles.cells(V,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),handles.actMapGeo(V),'LineStyle','none')
colormap('jet')
axis equal
caxis([10 45])
view(-135,0)
set(gca,'Visible','off')
hold on
trisurf(handles.cells(N,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),'FaceColor','none','LineWidth',0.1);
scatter3(handles.centroids(pInd(1),1),handles.centroids(pInd(1),2),handles.centroids(pInd(1),3),'mo','filled','SizeData',128)
scatter3(handles.centroids(pInd(2),1),handles.centroids(pInd(2),2),handles.centroids(pInd(2),3),'ko','filled','SizeData',128)
scatter3(handles.centroids(pInd(3),1),handles.centroids(pInd(3),2),handles.centroids(pInd(3),3),'co','filled','SizeData',128)

% APD80
figure
trisurf(handles.cells(V,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),handles.apdMapGeo(V),'LineStyle','none')
colormap('jet')
axis equal
caxis([60 110])
view(-135,0)
set(gca,'Visible','off')
hold on
trisurf(handles.cells(N,:),handles.pts(:,1),handles.pts(:,2),handles.pts(:,3),'FaceColor','none','LineWidth',0.1);
scatter3(handles.centroids(pInd(1),1),handles.centroids(pInd(1),2),handles.centroids(pInd(1),3),'mo','filled','SizeData',128)
scatter3(handles.centroids(pInd(2),1),handles.centroids(pInd(2),2),handles.centroids(pInd(2),3),'ko','filled','SizeData',128)
scatter3(handles.centroids(pInd(3),1),handles.centroids(pInd(3),2),handles.centroids(pInd(3),3),'co','filled','SizeData',128)


