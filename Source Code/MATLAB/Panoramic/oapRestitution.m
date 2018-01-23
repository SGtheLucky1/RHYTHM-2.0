function [] = oapRestitution(CL,region)
% Description: Find the average action potential duration (apd) based on 
% the average signal shape as well as an average of all apd values with 
% standard deviation.
%
% Input:
%   CL = cycle length can be substituted with concentration (conc.) for
%   drug response analysis
%   region = left ventricle (LV) or right ventricle (RV)
%
% Output:
%   None
%
% Author: Christopher Gloschat
% Date: August 31, 2016
% 
% Modification Log:
%
%
%% Code %%
% % % CL = 250;
% % % region = 'LV';
conc = CL;
CL = 300;
name = sprintf('%sregion_%duM.mat',region,conc);
load(name)
% save data to accurately named variable
data = M;
data = reshape(data,[size(data,1)*size(data,2),size(data,3),1])';
dataInd = 1:size(data,1)*size(data,2);
dataInd = reshape(dataInd,[size(data,1) size(data,2)]);
% clear old variable
clear M
% find first derivative
dVdt = data(2:end,:)-data(1:end-1,:);
% number of anticipated peaks - USER SPECIFIED
% % % numPeaks = 8;
% % % numPeaksUsed = 7;
% beginning and end indices
bind = -round(CL*0.15);
eind = CL+bind-1;
% preallocate variable for peak data
peaks = cell(size(dVdt,1),1);
% oaps = zeros(numPeaks,CL,size(dVdt,1));
oaps = cell(size(dVdt,1),1);
oapsTimeAve = zeros(size(dVdt,2),CL);
% find peaks of activation
for n = 1:size(dVdt,2)
    % activation times
    [tmp] = peakFinderGlo(data(:,1),CL,bind,eind);
    numPeaks = length(tmp);
% % %     tmp = peakfinder(dVdt(:,n),(max(dVdt(:,n))-min(dVdt(:,n)))/1.5);
% % % 
% % %     while length(tmp) > numPeaks
% % %         interval = tmp(2:end)-tmp(1:end-1);
% % %         % Is it outside the pacing cycle length
% % %         test = (interval < CL+1).*(interval > CL-1);
% % %         % Find the index
% % %         test = find(~test);
% % %         
% % %         if n == 17
% % %             disp(n)
% % %         end
% % %         
% % %         if ~isempty(test)
% % %             if test(1) == 1
% % %                 tmp(1) = [];
% % %             else
% % %                 tmp(test(1)+1) = [];
% % %             end
% % %         elseif 0 > tmp(1)+bind
% % %             tmp(1) = [];
% % %         end
% % %         
% % %         
% % %     end
%     peaks(n,:) = dataInd(peaks(n,:),n)';
    peaks{n} = dataInd(tmp,n)';
    % index for oaps
    oaps_ind = bind:eind;
    oaps_ind = repmat(oaps_ind,[numPeaks 1])+repmat(peaks{n}',[1 size(oaps_ind,2)]);
    % single oap
    oaps{n} = data(oaps_ind);
% %     for m = 1:numPeaks-1
% %         plot(f,squeeze(oaps(m,:,1)),'Color',[rand rand rand])
% %     end
    % average oaps
    oapsTimeAve(n,:) = mean(oaps{n});
end

% zero mean the baseline
baseline = mean(oapsTimeAve(:,1:bind/2*-1),2);
baseline = repmat(baseline,[1 size(oapsTimeAve,2)]);
oapsTimeAve = oapsTimeAve - baseline;
% % % Visualize a single average over time
% % figure
% % hold on
% % for n = 1:numPeaks-1
% %     plot(squeeze(oaps(n,:,1)),'Color',[rand rand rand])
% % end
% % % overlay average
% % plot(oapsTimeAve(1,:),'Color','b','LineWidth',4)
% Average spatially
oapsSpatAve = mean(oapsTimeAve);
f = figure;
hold on
plot(oapsSpatAve,'k--','LineWidth',4)
% output apd80 value
[maxAmp,maxInd] = max(oapsSpatAve);
% values ocurring after max
apd80Ind = 1:size(oapsSpatAve,2);
apd80Ind = apd80Ind > maxInd;
tmp = oapsSpatAve < maxAmp*0.2;
apd80Ind = (apd80Ind.*tmp).*(1:length(apd80Ind));
apd80Ind = unique(apd80Ind);
apd80Ind = apd80Ind(2);
% Find activation time
[~,actInd] = max(oapsSpatAve(2:end)-oapsSpatAve(1:end-1));
% Visualize
scatter(actInd,oapsSpatAve(actInd),'go','filled')
scatter(apd80Ind,oapsSpatAve(apd80Ind),'ro','filled')
hold off
% Calculate APD80
apd80 = apd80Ind-actInd;
legend(['apd80 = ' num2str(apd80) 'ms'],'Activation','Repolarizationn')
% % % tmp = sprintf('savefig(f,''%s_%dmsPlot'')',CL,region);
tmp = sprintf('savefig(f,''%s_%duMPlot'')',region,conc);
eval(tmp)

% find statistics
[maxAmpAve,maxIndAve] = max(oapsTimeAve,[],2);
apd80IndAve = repmat(1:size(oapsTimeAve,2),[size(oapsTimeAve,1) 1]);
apd80IndAve = apd80IndAve > repmat(maxIndAve,[1 size(apd80IndAve,2)]);
tmp = oapsTimeAve < repmat(maxAmpAve,[1 size(oapsTimeAve,2)])*0.2;
apd80IndAve = (apd80IndAve.*tmp).*repmat((1:size(apd80IndAve,2)),[size(apd80IndAve,1) 1]);
apd80IndAveFinal = zeros(size(apd80IndAve,1),1);
for n = 1:size(apd80IndAve,1)
    tmp = unique(apd80IndAve(n,:));
    apd80IndAveFinal(n) = tmp(2);
end
oapsTimeAct = oapsTimeAve(:,2:end)-oapsTimeAve(:,1:end-1);
[~,oapsTimeAct] = max(oapsTimeAct,[],2);
apd80Time = apd80IndAveFinal-oapsTimeAct;
apd80TimeAve = mean(apd80Time);
fprintf('The average action potential duration was %0.2f ms.\n',apd80TimeAve)
apd80TimeStd = std(apd80Time);
fprintf('The standard devation of the action potential duration was +/- %0.2f ms.\n',apd80TimeStd)
save('Transfer','conc','apd80TimeAve','apd80TimeStd','numPeaks')

end






 