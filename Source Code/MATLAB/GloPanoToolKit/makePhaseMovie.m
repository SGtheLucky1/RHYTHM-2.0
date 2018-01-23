% 3D Movie Generator %

% Create figure for visualization
scrn_size = get(0,'ScreenSize');
f = figure('Name','PHASE ANALYSIS','Position',...
    [scrn_size(3)/2-300,scrn_size(4)/2-375,1300,500],'NumberTitle','Off');

% % % 
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

dataScrn = axes('Units','Pixels','Visible','off','XTick',[],'YTick',[],'Position',[25 50 400 400]);
phaseScrn = axes('Units','Pixels','Visible','off','XTick',[],'YTick',[],'Position',[450 50 400 400]);
ntScrn = axes('Units','Pixels','Visible','off','XTick',[],'YTick',[],'Position',[875 50 400 400]);

m = VideoWriter('test.avi');
m.FrameRate = 50;
open(m);

% % figure
for n = 1:2047
    % Data axes
    trisurf(data.cells(V,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        data.dataProj(V,n),'LineStyle','none','Parent',dataScrn)
    hold(dataScrn,'on')
    trisurf(data.cells(N,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        'FaceColor','none','Parent',dataScrn)
    hold(dataScrn,'off')
    set(dataScrn,'Visible','off')
    colormap(dataScrn,'jet')
    caxis(dataScrn,[0 1])
    axis(dataScrn,'equal')
    view(dataScrn,-180,0)
    
    
    % Phase axes
    trisurf(data.cells(V,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        phaseGeo(V,n),'LineStyle','none','Parent',phaseScrn)
    hold(phaseScrn,'on')
    trisurf(data.cells(N,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        'FaceColor','none','Parent',phaseScrn)
    hold(phaseScrn,'off')
    set(phaseScrn,'Visible','off')
    colormap(phaseScrn,'jet')
    caxis(phaseScrn,[-pi pi])
    axis(phaseScrn,'equal')
    view(phaseScrn,-180,0)


    % Singularities axes
    trisurf(data.cells(V,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        ntGeo(V,n),'LineStyle','none','Parent',ntScrn)
    hold(ntScrn,'on')
    trisurf(data.cells(N,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
        'FaceColor','none','Parent',ntScrn)
    hold(ntScrn,'off')
    set(ntScrn,'Visible','off')
    colormap(ntScrn,'jet')
    caxis(ntScrn,[-1 1])
    axis(ntScrn,'equal')
    view(ntScrn,-180,0)
    
    frame = getframe(f);
    writeVideo(m,frame);
end
close(m)


% % % % Representative OAP
% % % repScrn = axes('Units','Pixels','Visible','off','Position',[50 25 500 150]);
% % % t = 50:349;
% % % cmap = colormap('jet');
% % % c = phase(ind,t);
% % % cn = (c-min(c))/max(c-min(c));
% % % cn = ceil(cn*size(cmap,1));
% % % cn = max(cn,1);
% % % colormap('jet')
% % % hold on
% % % axes(repScrn)
% % % for n = 1:length(t)-1
% % % line(t(n:n+1),data.dataProj(ind,t(n):t(n+1)),'color',cmap(cn(n),:),'lineWidth',8)
% % % end
% % % hold off
% % % % Colorbar settings
% % % set(repScrn,'XLim',[60 345])
% % % cb = colorbar('Location','SouthOutside');
% % % set(cb,'XDir','reverse','FontSize',24,'Ticks',0:0.25:1,'TickLabels',...
% % %     {'-\pi','-\pi/2','0','\pi/2','\pi'})
% % % 
% % % phaseScrn = axes('Units','Pixels','YTick',[],'XTick',[],'Position',[50 200 500 500]);
% % % axes(phaseScrn)
% % % 
% % % m = VideoWriter('test.avi');
% % % m.FrameRate = 50;
% % % open(m);
% % % 
% % % for n = 75:420
% % %     t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
% % %         phase(vertInd,n),'LineStyle','none','Parent',phaseScrn);
% % % %     alpha(t,0.5)
% % %     hold on
% % %     s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
% % %         data.pts(:,3),'FaceColor','none','LineWidth',0.1,'Parent',phaseScrn);
% % %     hold off
% % %     set(phaseScrn,'Visible','off')
% % %     caxis([-pi pi])
% % %     colormap('jet')
% % %     axis equal
% % %     view(-45,0)
% % %     frame = getframe(f);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 1:180
% % %     view(-45+n,0)
% % %     frame = getframe(f);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 421:1020
% % %     t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
% % %         phase(vertInd,n),'LineStyle','none','Parent',phaseScrn);
% % % %     alpha(t,0.5)
% % %     hold on
% % %     s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
% % %         data.pts(:,3),'FaceColor','none','LineWidth',0.1,'Parent',phaseScrn);
% % %     hold off
% % %     set(phaseScrn,'Visible','off')
% % %     caxis([-pi pi])
% % %     colormap('jet')
% % %     axis equal
% % %     view(135,0)
% % %     frame = getframe(f);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 1:180
% % %     view(135+n,0)
% % %     frame = getframe(f);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % for n = 1021:1200
% % %     t = trisurf(data.cells(vertInd,:),data.pts(:,1),data.pts(:,2),data.pts(:,3),...
% % %         phase(vertInd,n),'LineStyle','none','Parent',phaseScrn);
% % % %     alpha(t,0.5)
% % %     hold on
% % %     s = trisurf(data.cells(notInd,:),data.pts(:,1),data.pts(:,2),...
% % %         data.pts(:,3),'FaceColor','none','LineWidth',0.1,'Parent',phaseScrn);
% % %     hold off
% % %     set(phaseScrn,'Visible','off')
% % %     caxis([-pi pi])
% % %     colormap('jet')
% % %     axis equal
% % %     view(-45,0)
% % %     frame = getframe(f);
% % %     writeVideo(m,frame);
% % % end
% % % 
% % % close(m)