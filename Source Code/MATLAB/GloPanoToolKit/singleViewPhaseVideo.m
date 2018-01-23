

scrn_size = get(0,'ScreenSize');
f = figure('Name','Phase Singularity Trajectories','Position',...
    [scrn_size(3)/2-300,scrn_size(4)/2-375,700,500],'NumberTitle','Off');

slopeMap = axes('Parent',f,'Units','Pixels','Visible','on','XTick',[],'YTick',[],'Position',[150 50 400 400]);
smPos = get(slopeMap,'Position');
colorbar('Location','eastoutside')
set(gca,'XTick',[],'YTick',[],'Position',smPos)

m = VideoWriter('test.avi');
% m.FrameRate = 150;
m.FrameRate = 20;
open(m);
for n = 10090:5:10300
    imagesc(phase{1}(:,:,n))
    colormap('jet')
    set(gca,'XTick',[],'YTick',[])
    caxis([-pi pi])   
    cb = colorbar('Ticks',-pi:pi/2:pi,'TickLabels',{['-' char(960)],['-' char(960) '/2'],'0',[char(960) '/2'],char(960)});
    set(gca,'Position',smPos,'FontSize',20)
    frame = getframe(f);
    writeVideo(m,frame);
end
close(m)