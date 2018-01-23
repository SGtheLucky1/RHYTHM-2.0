% Visualize Restitution Characteristics

% Find highest slope signals in RV
[topVal,topInd] = sort(maxSlope{3},'descend');
check = unique(~isnan(topVal).*(1:length(topVal))');
topInd = topInd(check(2:6));
topVal = maxSlope{3}(topInd);

figure
for n = 1:5
    subplot(5,1,n)
    plot(data{3}{6}(topInd(n),:))
    set(gca,'XLim',[0 size(data{3}{6},2)],'YLim',[-0.1 1])
    if n < 5
        set(gca,'XTickLabel',[])
    else
        xlabel('Time (s)')
        set(gca,'XTickLabel',0:0.5:4)
    end
    
    if n == 3
        ylabel('Fluorescence Intensity')
    end
end


% maxSlope vs. APD_std
figure
t = cell(1,4);
t{1} = 'LV';
t{2} = 'Anterior';
t{3} = 'RV';
t{4} = 'Posterior';
for n = 1:4
    subplot(2,2,n)
    scatter(maxSlope{n},APD_std{n}{6})
    set(gca,'XLim',[0.5 2],'YLim',[0 40])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Standard deviation of APD')
end

% Compare two curves 1) long APD and 2) short APD at long cycle lengths
[minCompVal,minCompInd] = nanmin(maxSlope{3});
[maxCompVal,maxCompInd] = nanmax(maxSlope{3});
figure
s1 = scatter(DI_ave{3}(minCompInd,:),APD_ave{3}(minCompInd,:),'bo');
hold on
s2 = scatter(DI_ave{3}(maxCompInd,:),APD_ave{3}(maxCompInd,:),'co');
set(gca,'XLim',[20 180])
fitX = 20:0.1:180;
fitYmax = f{3}{maxCompInd}(fitX);
fitYmin = f{3}{minCompInd}(fitX);
p2 = plot(fitX,fitYmax,'r');
p1 = plot(fitX,fitYmin,'m');
legend([p2,p1],[['R^{2} = ' num2str(gof{3}{maxCompInd}.rsquare)];['R^{2} = ' num2str(gof{3}{minCompInd}.rsquare)]],'Location','southeast');


% maxSlope vs. APD_ave
figure
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(APD_ave{n}(:,6));
    check(:,2) = ~isnan(maxSlope{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(maxSlope{n}(keep),APD_ave{n}(keep,6),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
    h1 = scatter(maxSlope{n},APD_ave{n}(:,6));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
    set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Average APD @ 150 ms CL')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end


