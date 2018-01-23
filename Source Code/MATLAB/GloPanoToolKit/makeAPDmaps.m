% apdMaps

% % % load('12_processedSansGeoSansAtriaAnalyzed')


% Camera view labels
t = cell(1,4);
t{1} = 'LV';
t{2} = 'Anterior';
t{3} = 'RV';
t{4} = 'Posterior';
%% Populate image variables %%
maxSlopeImage = nan(100,100,4);
goodFitImage = nan(100,100,4);
Rk = cell(1,4);
RkImage = nan(100,100,4);
for n = 1:4
    % Max slope
    ind = ~isnan(maxSlope{n});
    maxSlopeImage(dataInd{n}(ind)+(n-1)*10000) = maxSlope{n}(ind);
    % Goodness of fit
    for m = 1:length(gof{n})
        if ~isempty(gof{n}{m})
            goodFitImage(dataInd{n}(m,n)) = gof{n}{m}.rsquare;
        end
    end
    % Rate constant
    Rk{n} = APD_ave{n}(:,1)./DI_ave{n}(:,1);
    ind = ~isnan(Rk{n});
    RkImage(dataInd{n}(ind)+(n-1)*10000) = Rk{n}(ind);
end

%% Perform erode/dilate then average the internal values %%
% Create mask variables for each image
maxSlopeMask = ~isnan(maxSlopeImage);
goodFitMask = ~isnan(goodFitImage);
RkMask = ~isnan(RkImage);

% % % before = figure('Name','Before Erode/Dilate');
% % % for m = 1:3
% % %     for n = 1:4
% % %         subplot(3,4,4*(m-1)+n)
% % %         if m == 1
% % %             imagesc(maxSlopeMask(:,:,n))
% % %         elseif m == 2
% % %             imagesc(goodFitMask(:,:,n))
% % %         else
% % %             imagesc(RkMask(:,:,n))
% % %         end
% % %     end
% % % end
% Erode and dilate each one using a 
se = strel('square',3);
maxSlopeMask = imerode(maxSlopeMask,se);
maxSlopeMask = imdilate(maxSlopeMask,se);
goodFitMask = imerode(goodFitMask,se);
goodFitMask = imdilate(goodFitMask,se);
RkMask = imerode(RkMask,se);
RkMask = imdilate(RkMask,se);
% Identify largest connect component for each mask
for m = 1:3
    for n = 1:4
        if m == 1
            CC = bwconncomp(maxSlopeMask(:,:,n));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            maxSlopeMask(:,:,n) = zeros(100,100);
            maxSlopeMask(CC.PixelIdxList{idx}+(n-1)*10000) = 1;
        elseif m ==2
            CC = bwconncomp(goodFitMask(:,:,n));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            goodFitMask(:,:,n) = zeros(100,100);
            goodFitMask(CC.PixelIdxList{idx}+(n-1)*10000) = 1;
        else
            CC = bwconncomp(RkMask(:,:,n));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            RkMask(:,:,n) = zeros(100,100);
            RkMask(CC.PixelIdxList{idx}+(n-1)*10000) = 1;
        end
    end
end

% % % % After erode/dilate
% % % after = figure('Name','After Erode/Dilate');
% % % for m = 1:3
% % %     for n = 1:4
% % %         subplot(3,4,4*(m-1)+n)
% % %         if m == 1
% % %             imagesc(maxSlopeMask(:,:,n))
% % %         elseif m == 2
% % %             imagesc(goodFitMask(:,:,n))
% % %         else
% % %             imagesc(RkMask(:,:,n))
% % %         end
% % %     end
% % % end


% Max slope
maxSlopeImage(maxSlopeMask == 0) = nan;
% Goodness of fit
goodFitImage(goodFitMask == 0) = nan;
% Rate constant
RkImage(RkMask == 0) = nan;


%% Visualize %%
figMaxSlope = figure('Name','Maximum Slope of Restitution');
figGOF = figure('Name','Max Slope Goodness of Fit');
figRk = figure('Name','Rate Constant R_k');
for n = 1:4
    % Max slope
    figure(figMaxSlope)
    subplot(2,2,n)
    imagesc(maxSlopeImage(:,:,n))
    caxis([0 2])
    colormap('jet')
    colorbar
    set(gca,'XTick',[],'YTick',[])
    % Goodness of fit
    figure(figGOF)
    subplot(2,2,n)
    imagesc(goodFitImage(:,:,n))
    colormap('jet')
    caxis([0 1])
    colorbar
    set(gca,'XTick',[],'YTick',[])
    % Rate constant
    figure(figRk)
    subplot(2,2,n)
    imagesc(RkImage(:,:,n))
    colormap('jet')
    caxis([0 1.5])
    colorbar
    set(gca,'XTick',[],'YTick',[])
end

% Relationship between slope and rate constant
compSlopeRate = figure('Name','Comparison of Slope and Rate Constant');
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(maxSlope{n});
    check(:,2) = ~isnan(Rk{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(maxSlope{n}(keep),Rk{n}(keep),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(maxSlope{n}(keep))*10)/10:0.1:round(nanmax(maxSlope{n}(keep))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},Rk{n});
    h1 = scatter(maxSlope{n}(keep),Rk{n}(keep));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Rate Constant')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end

% maxSlope vs. APD_ave @ 150 ms CL
compSlopeAPD_150 = figure('Name','Comparison of Slope and APD80 @ 150 ms CL');
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
    fitX1{n} = round(nanmin(maxSlope{n}(keep))*10)/10:0.1:round(nanmax(maxSlope{n}(keep))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},APD_ave{n}(:,6));
    h1 = scatter(maxSlope{n}(keep),APD_ave{n}(keep,6));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
% % %     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Average APD @ 150 ms CL')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end

% maxSlope vs. APD_ave @ 300 ms CL
compSlopeAPD_300 = figure('Name','Comparison of Slope and APD80 @ 300 ms CL');
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(APD_ave{n}(:,1));
    check(:,2) = ~isnan(maxSlope{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(maxSlope{n}(keep),APD_ave{n}(keep,1),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(maxSlope{n}(keep))*10)/10:0.1:round(nanmax(maxSlope{n}(keep))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},APD_ave{n}(:,1));
    h1 = scatter(maxSlope{n}(keep),APD_ave{n}(keep,1));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Average APD @ 300 ms CL')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end


% maxSlope vs. DI_ave @ 150 ms CL
compSlopeDI_150 = figure('Name','Comparison of Slope and DI @ 150 ms CL');
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(DI_ave{n}(:,6));
    check(:,2) = ~isnan(maxSlope{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(maxSlope{n}(keep),DI_ave{n}(keep,6),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(maxSlope{n}(keep))*10)/10:0.1:round(nanmax(maxSlope{n}(keep))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},DI_ave{n}(:,6));
    h1 = scatter(maxSlope{n}(keep),DI_ave{n}(keep,6));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Diastolic Interval @ 150 ms CL')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end

% maxSlope vs. DI_ave @ 300 ms CL
compSlopeDI_300 = figure('Name','Comparison of Slope and DI @ 300 ms CL');
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(DI_ave{n}(:,1));
    check(:,2) = ~isnan(maxSlope{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(maxSlope{n}(keep),DI_ave{n}(keep,1),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(maxSlope{n}(keep))*10)/10:0.1:round(nanmax(maxSlope{n}(keep))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},DI_ave{n}(:,1));
    h1 = scatter(maxSlope{n}(keep),DI_ave{n}(keep,1));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Max Slope')
    ylabel('Diastolic Interval @ 300 ms CL')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end

% Rate Constant vs. APD_ave @ 300 ms CL
compRateAPD80_300 = figure('Name','Comparison of Rate Constant and APD80 @ 300 ms CL');
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(APD_ave{n}(:,1));
    check(:,2) = ~isnan(Rk{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(APD_ave{n}(keep,1),Rk{n}(keep),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(APD_ave{n}(keep,1))*10)/10:0.1:round(nanmax(APD_ave{n}(keep,1))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},Rk{n});
    h1 = scatter(APD_ave{n}(keep,1),Rk{n}(keep));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Average APD_8_0 @ 300 ms')
    ylabel('Rate Constant')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end


% Rate Constant vs. DI_ave @ 300 ms CL
compRateDI_300 = figure('Name','Comparison of Rate Constant and APD80 @ 300 ms CL');
% Fit variables
ft1 = fittype('poly1');
f1 = cell(4,1);
gof1 = cell(4,1);
output1 = cell(4,1);
fitX1 = cell(4,1);
fitY1 = cell(4,1);
for n = 1:4
    % Calculate fit
    check = ~isnan(DI_ave{n}(:,1));
    check(:,2) = ~isnan(Rk{n});
    keep = logical(check(:,1).*check(:,2));
    [f1{n},gof1{n},output1{n}] = fit(DI_ave{n}(keep,1),Rk{n}(keep),ft1);
    % Calculate the X and Y values, then calculate the maximum slope
    fitX1{n} = round(nanmin(DI_ave{n}(keep,1))*10)/10:0.1:round(nanmax(DI_ave{n}(keep,1))*10)/10;
% % %     fitX1{n} = 0.5:0.1:2;
    fitY1{n} = f1{n}(fitX1{n})';
    
    %Visualize
    subplot(2,2,n)
% % %     h1 = scatter(maxSlope{n},Rk{n});
    h1 = scatter(DI_ave{n}(keep,1),Rk{n}(keep));
    hold on
    h2 = plot(fitX1{n},fitY1{n},'r');
%     set(gca,'XLim',[0.5 2],'YLim',[80 130])
    title(t{n})
    xlabel('Average DI_8_0 @ 300 ms')
    ylabel('Rate Constant')
    legend(h2,['R^{2} = ' num2str(gof1{n}.rsquare)],'Location','southeast');
end
