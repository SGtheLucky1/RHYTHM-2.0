function [f,gof,output,fitX,fitY,maxSlope,maxSlopeInd,track] = resFit(APD_ave,DI_ave)
% Calculate a fit for each pixel
ft = cell(1,2);
f = cell(1,3);
gof = cell(1,3);
output = cell(1,3);
% % ft = fittype('exp2');
ft{1} = fittype('a*(1-exp(-b*x))+c');
f{1} = cell(size(APD_ave,1),1);
gof{1} = cell(size(APD_ave,1),1);
output{1} = cell(size(APD_ave,1),1);

ft{2} = fittype('((a-b)./(1+exp(-k*(x-c))))+b');
f{2} = cell(size(APD_ave,1),1);
gof{2} = cell(size(APD_ave,1),1);
output{2} = cell(size(APD_ave,1),1);

maxSlope = zeros(size(APD_ave,1),1);
maxSlopeInd = zeros(size(APD_ave,1),1);
fitX = cell(size(APD_ave,1),1);
fitY = cell(size(APD_ave,1),1);
h = waitbar(0,'Fitting data.');
track = zeros(size(APD_ave,1),1);
for n = 1:size(APD_ave,1)
    if ~isnan(APD_ave(n,1))
        % Remove data points whose DI is greater their predecessor
        ind = find(~isnan(DI_ave(n,:)));
        check = DI_ave(n,ind(1:end-1))-DI_ave(n,ind(2:end));
        check = check < 5;
        if sum(check) > 0
            check = ind(find(check)+1);
            APD_ave(n,check) = nan;
            DI_ave(n,check) = nan;
        end
        if sum(isnan(DI_ave(n,:))) < 3
            pts = ~isnan(APD_ave(n,:));
            [f{1}{n},gof{1}{n},output{1}{n}] = fit(DI_ave(n,pts)',APD_ave(n,pts)',ft{1},'StartPoint',[171 0.04 0]);
            [f{2}{n},gof{2}{n},output{2}{n}] = fit(DI_ave(n,pts)',APD_ave(n,pts)',ft{2},'StartPoint',[160 120 80 0.1]);
            % Calculate x-axis values
            fitX{n} = min(DI_ave(n,:)):0.1:max(DI_ave(n,:));
            % Calculate y-axis values
            if gof{1}{n}.rsquare > 0.9
                fitY{n} = f{1}{n}(fitX{n})';
                track(n) = 1;
            elseif gof{1}{n}.rsquare < gof{2}{n}.rsquare
                fitY{n} = f{2}{n}(fitX{n})';
                track(n) = 2;
            else
                fitY{n} = f{1}{n}(fitX{n})';
                track(n) = 1;
            end
            % Calculate maximum slope
            [maxSlope(n),maxSlopeInd(n)] = max((fitY{n}(2:end)-fitY{n}(1:end-1))./(fitX{n}(2:end)-fitX{n}(1:end-1)));
        else
            maxSlope(n) = nan;
            maxSlopeInd(n) = nan;
        end
    end
    % Update waitbar
    waitbar(n/size(APD_ave,1))
end
close(h)
tmp = find(track ~= 0);
f{3} = cell(length(f{2}),1);
gof{3} = cell(length(gof{2}),1);
output{3} = cell(length(output{2}),1);
for n = 1:length(tmp)
    disp(n)
    f{3}{tmp(n)} = f{track(tmp(n))}{tmp(n)};
    gof{3}{tmp(n)} = gof{track(tmp(n))}{tmp(n)};
    output{3}{tmp(n)} = output{track(tmp(n))}{tmp(n)};
end
end


