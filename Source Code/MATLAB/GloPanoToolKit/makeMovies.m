% 3D Movie Generator %

% Locations of nan vertices
N = data.dataProj(:,1) == 0;
notInd = find(N);
% Locations of vertices with data
V = N == 0;
vertInd = find(V);
% Remove islands
rmCC = removeIslands(data.cells(V,:));
tmp = find(rmCC(:,2) ~= 1);
rmCC = rmCC(tmp,1);
notInd = [notInd;vertInd(rmCC,:)];
vertInd(rmCC,:) = [];



m = VideoWriter('test.avi');
m.FrameRate = 75;
open(m);

figure
for n = 75:420
    if n == 120
        disp(n)
    end
    
    t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        phase(vertInd,n),'LineStyle','none');
    alpha(t,0.5)
    hold on
    s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
        data.pts(:,3),'FaceColor','none','LineWidth',0.1);
    hold off
    set(gca,'Visible','off')
    caxis([-pi pi])
    colormap('jet')
    axis equal
    view(-45,0)
    frame = getframe(gcf);
    writeVideo(m,frame);
end

for n = 1:361
    view(-45+(n-1),0)
    frame = getframe(gcf);
    writeVideo(m,frame);
end

% % % for n = 421:1020
% % %     t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
% % %         phase(vertInd,n),'LineStyle','none');
% % %     alpha(t,0.5)
% % % %     hold on
% % % %     s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
% % % %         data.pts(:,3),'FaceColor','none','LineWidth',0.1);
% % % %     hold off
% % %     set(gca,'Visible','off')
% % %     caxis([-pi pi])
% % %     colormap('jet')
% % %     axis equal
% % %     view(135,0)
% % %     frame = getframe(gcf);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 1:180
% % %     view(135+n,0)
% % %     frame = getframe(gcf);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 1021:1200
% % %     t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
% % %         phase(vertInd,n),'LineStyle','none');
% % %     alpha(t,0.5)
% % % %     hold on
% % % %     s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
% % % %         data.pts(:,3),'FaceColor','none','LineWidth',0.1);
% % % %     hold off
% % %     set(gca,'Visible','off')
% % %     caxis([-pi pi])
% % %     colormap('jet')
% % %     axis equal
% % %     view(-45,0)
% % %     frame = getframe(gcf);
% % %     writeVideo(m,frame);
% % % end

close(m)