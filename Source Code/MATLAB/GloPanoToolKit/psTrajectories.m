%% maintanence visualize %%

% Create heat map figure
scrn_size = get(0,'ScreenSize');
f = figure('Name','Phase Singularity Trajectories','Position',...
    [scrn_size(3)/2-300,scrn_size(4)/2-375,700,500],'NumberTitle','Off');

slopeMap = axes('Parent',f,'Units','Pixels','Visible','on','XTick',[],'YTick',[],'Position',[150 50 400 400]);
smPos = get(slopeMap,'Position');
trajMap = axes('Parent',f,'Units','Pixels','Visible','off','XTick',[],'YTick',[],'Position',[150 50 400 400]);
tmPos = get(trajMap,'Position');

axes(slopeMap)
imagesc(maxSlopeImage(:,:,1));
caxis([0 2])
colormap(slopeMap,'jet')
colorbar('Location','eastoutside')
set(gca,'XTick',[],'YTick',[],'Position',smPos)


% Create phase singularity trajectory path
rot1 = (9160:10:9410)';
rot1(:,2) = [61,58,57,56,53,51,49,50,47,41,36,nan,nan,40,50,50,48,50,52,51,50,54,52,51,50,48]';
rot1(:,3) = [37,36,38,40,41,42,43,53,55,50,50,nan,nan,41,36,38,37,37,36,38,40,40,43,44,45,44]';

rot2 = (9410:10:9650)';
rot2(:,2) = [48,51,51,48,47,47,45,46,45,45,41,45,47,nan,nan,nan,42,41,nan,39,nan,36,36,34,48]';
rot2(:,3) = [44,46,49,55,54,55,61,63,64,64,62,63,60 ,nan,nan,nan,38,40,nan,39,nan,42,43,43,57]';

rot3 = (9660:10:9870)';
rot3(:,2) = [50,49,46,51,48,52,51,52,52,46,47,51,49,48,52,51,48,46,44,41,37,38]';
rot3(:,3) = [56,58,62,64,64,63,62,61,61,60,38,36,38,37,35,35,38,38,40,41,43,44]';

rot4 = (9880:10:10090)';
rot4(:,2) = [38,54,55,56,58,61,63,60,64,59,59,57,58,58,54,52,51,48,46,45,43,42]';
rot4(:,3) = [47,52,54,52,50,50,52,55,50,51,50,46,38,40,37,38,37,38,40,38,40,40]';

rot5 = (10100:10:10300)';
rot5(:,2) = [49,47,45,50,52,44,50,57,54,56,58,58,57,58,59,58,54,54,54,53,55]';
rot5(:,3) = [56,58,56,57,57,54,51,49,53,49,48,49,48,48,45,46,48,47,47,48,49]';


% Map path of PS onto heat map
cmap = colormap('gray');
c = rot5(:,1);
cn = (c-min(c))/max(c-min(c));
cn = ceil(cn*size(cmap,1));
cn = max(cn,1);

% Remove missing data
rm = find(sum(isnan(rot5),2));
rot5(rm,:) = [];

axes(trajMap)
set(gca,'XLim',[1 100],'YLim',[1 100],'YDir','reverse')
colormap(trajMap,'gray')
hold on
for n = 1:size(rot5,1)-1
    line(rot5(n:n+1,3),rot5(n:n+1,2),'color',cmap(cn(n),:),'lineWidth',4)
end
hold off
% Colorbar settings
cb = colorbar('Location','westoutside');
set(cb,'Ticks',0:1/(size(rot5,1)/2-1):1,'TickLabels',rot5(1:2:end,1))
set(gca,'Position',tmPos)