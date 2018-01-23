%% Algorithm for generating restitution curves %%
function [APD_ave,APD80,DI_ave,DI,DATA,dataInd,act_I] = restApdCalc(cmosData,ecg,nRate)

% % % load('08_processedSansGeo')

% % % cmosData = a.cmosData;
% % % nRate = a.nRate;
% % % ecg = a.ecg;

% Data dimensions
x = size(cmosData,1);
y = size(cmosData,2);
z = size(cmosData,3);

% Identify time when the paced data is equal to zero
up = find(ecg == 0);
% Each spike is actually a step with multiple values, subtract to find the
% jumps from one spike to the next
front = up(2:end)-up(1:end-1);
% Add 1 to the subtracted values to get the values of the next steps
up = up([1; find(front~=1)+1]);
% Pacing rate
% % % CL = round(mean((up(2:end)-up(1:end-1))/nRate));
% Pacing times at 1kHz timing
pacedInd = round(up/nRate);

% Grab non-NAN data indices (:,:,1) and extrapolate to 3D format
dataInd = ~isnan(cmosData(:,:,1));
% Erode the data to eliminate border noise
se = strel('square',3);
dataInd = find(imerode(dataInd,se));
dataInd = repmat(dataInd,[1 z])+repmat(0:x*y:x*y*z-1,[size(dataInd,1) 1]);

% Grab data in 2D format
data = cmosData(dataInd);
data = data(:,1:end);

% % % ecg = ecg(:,500*nRate:end);[
% % % pacedInd = pacedInd(logical((pacedInd > 499)))-499;

% Calculate noise
noise = max(data(:,1:pacedInd(1)),[],2)-min(data(:,1:pacedInd(1)),[],2);
signal = max(data(:,1:pacedInd(2)),[],2)-min(data(:,1:pacedInd(2)),[],2);
snr = signal./noise;

% Exclude all signals less than 1.5 times the standard deviation
ex = snr < mean(snr)-1.5*std(snr);
data(ex,:) = [];
dataInd(ex,:) = [];
% signal(ex) = [];
% noise(ex) = [];
% snr(ex) = [];

% Group paced beats based on cycle length
pacedCL = pacedInd(2:end)-pacedInd(1:end-1);
pacedCL = [pacedCL(1);pacedCL];
pacedTransInd = [1; find((pacedCL(2:end)-pacedCL(1:end-1)) == -30)+1];
pacedTrans = pacedInd(pacedTransInd);

% Break up the data by CL
dataByCL = cell(length(pacedTrans),1);
paceByCL = cell(length(pacedTrans),1);
for n = 1:length(pacedTrans)
    if n < length(pacedTrans)
        % Find the desired end point (ep)
% %         [~,ep] = min(data(:,pacedTrans(n+1):pacedTrans(n+1)+50),[],2);
% %         ep = round(mean(ep));
        dataByCL{n} = data(:,pacedTrans(n)-20:pacedTrans(n+1)+50);
        paceByCL{n} = pacedInd(pacedTransInd(n):pacedTransInd(n+1)-1)-pacedInd(pacedTransInd(n))+20;
    else
        dataByCL{n} = data(:,pacedTrans(n)-20:end);
        paceByCL{n} = pacedInd(pacedTransInd(n):end)-pacedInd(pacedTransInd(n))+20;
    end
end

% Find restitution values at each CL
dVdt = cell(length(pacedTrans),1);
act_I = cell(length(pacedTrans),1);
max_V = cell(length(pacedTrans),1);
max_I = cell(length(pacedTrans),1);
max_C = cell(length(pacedTrans),1);
repol = cell(length(pacedTrans),1);
APD80 = cell(length(pacedTrans),1);
DI = cell(length(pacedTrans),1);
DATA = cell(length(pacedTrans),1);
for n = 1:length(pacedTrans)
% % %     disp(n)
    % Grab data
    d = dataByCL{n};
    p = paceByCL{n};
    %  Remove the baseline
    [~,winF] = min(d(:,p(1):p(1)+49),[],2);
    winF = winF+p(1)-1;
    winFind = sub2ind([size(d,1) size(d,2)],(1:size(d,1))',winF);
    [~,winB] = min(d(:,p(end):p(end)+49),[],2);
    winB = winB+p(1)-1;
    winBind = sub2ind([size(d,1),size(d,2)],(1:size(d,1))',winB);
    slope = (d(winFind)-d(winBind))/(p(end)-p(1));
    slope = repmat(slope,[1 size(d,2)]).*repmat(0:size(d,2)-1,[size(d,1) 1]);
    d = d-slope;
    % Adjust the baseline
    % CHANGE back to 1:10 WHEN DONE
    winC = repmat(p(1:2),[1 50 size(d,1)])+repmat(0:49,[2 1 size(d,1)]);
    winR = repmat(reshape(1:size(d,1),[1 1 size(d,1)]),[size(winC,1) size(winC,2) 1]);
    winI = sub2ind([size(d,1),size(d,2)],winR,winC);
    winMin = reshape(mean(min(d(winI),[],2),1),[size(d,1) 1]);
%     base = repmat(mean(d(:,p(1:10)),2),[1 size(d,2)]);
    base = repmat(winMin,[1 size(d,2)]);
    d = d-base;
    
    % First derivative of the data
    dVdt{n} = d(:,2:end)-d(:,1:end-1);
    % Second derivative of the data
%     dV2dt = dVdt(:,2:end)-dVdt(:,1:end-1);
    % Preallocate variables
    act_I{n} = zeros(size(d,1),length(p)-1);
    max_V{n} = zeros(size(d,1),length(p)-1);
    max_I{n} = zeros(size(d,1),length(p)-1);
    max_C{n} = zeros(size(d,1),length(p)-1);
    repol{n} = zeros(size(d,1),length(p)-1);
    APD80{n} = zeros(size(d,1),length(p)-1);
    DI{n} = zeros(size(d,1),length(p)-1);

    % Find the max dVdt after each pacing spike at each pixel
    for m = 1:length(p)
        % Find the activation time (max dVdt)
        [~,act_I{n}(:,m)] = max(dVdt{n}(:,p(m):p(m)+49),[],2);
        act_I{n}(:,m) = act_I{n}(:,m)+p(m);
        if m < length(p)
% % %             % Find the maximum value
% % %             [max_V{n}(:,m),max_I] = max(d(:,p(m):p(m+1)),[],2);
% % %             max_I = max_I+p(m)-1;
            % Find the top of the upstroke by finding when the slope
            % switches its sign (+ -> -)
            ind = repmat(1:size(dVdt{n},2),[size(dVdt{n},1) 1]) >= repmat(act_I{n}(:,m),[1 size(dVdt{n},2)]);
            ind = ind.*(dVdt{n} < 0);
            ind = ind.*repmat(1:size(dVdt{n},2),[size(dVdt{n},1) 1]);
            ind(ind==0) = nan;
            [~,max_C{n}(:,m)] = nanmin(ind,[],2);
            max_I{n}(:,m) = sub2ind([size(d,1) size(d,2)],(1:size(d,1))',max_C{n}(:,m));
            % Find the APD80 repolarization index
            apdTarget = d(max_I{n}(:,m))*0.2; 
            for k = 1:size(d,1)
                % % %                 % Check to see if previously excluded
                % % %                 if m > 1
                % % %                     if isnan(repol{n}(k,m-1))
                % % %                         repol{n}(k,m) = nan;
                % % %                     end
                % % %                 end
                % % %                 % If not excluded proceed with check
                % % %                 if ~isnan(repol{n}(k,m))
                % Grab window to check
                apdWin = find(d(k,max_C{n}(k,m):p(m+1)) < apdTarget(k));
                if isempty(apdWin)
                    repol{n}(k,m) = nan;
                else
                    repol{n}(k,m) = max_C{n}(k,m)+apdWin(1)-1;
                end
                % % %                 end
            end
            % Calcuate APD80
            APD80{n}(:,m) = repol{n}(:,m)-act_I{n}(:,m);
        end
        if m > 1
            DI{n}(:,m-1) = act_I{n}(:,m)-repol{n}(:,m-1);
        end
    end
    DATA{n} = d;
end

APD_ave = zeros(size(data,1),size(APD80,1));
APD_std = zeros(size(data,1),size(APD80,1));
DI_ave = zeros(size(data,1),size(APD80,1));
% Calculate averages
for n = 1:size(APD80,1)
    APD_ave(:,n) = nanmean(APD80{n}(:,6:end),2);
    APD_std(:,n) = nanstd(APD80{n}(:,6:end),0,2);
    DI_ave(:,n) = nanmean(DI{n}(:,6:end),2);
end

% Clean up data
% All values must be within 2 std of the mean at that CL
APD_aveOverall = nanmean(APD_ave);
APD_stdOverall = nanstd(APD_ave);
APD_up = repmat(APD_aveOverall+2*APD_stdOverall,[size(APD_ave,1) 1]);
APD_dn = repmat(APD_aveOverall-2*APD_stdOverall,[size(APD_ave,1) 1]);
APD_rm = logical((APD_up < APD_ave)+(APD_dn > APD_ave));
APD_ave(APD_rm) = nan;
DI_ave(APD_rm) = nan;

DI_aveOverall = nanmean(DI_ave);
DI_stdOverall = nanstd(DI_ave);
DI_up = repmat(DI_aveOverall+2*DI_stdOverall,[size(DI_ave,1) 1]);
DI_dn = repmat(DI_aveOverall-2*DI_stdOverall,[size(DI_ave,1) 1]);
DI_rm = logical((DI_up < DI_ave)+(DI_dn > DI_ave));
DI_ave(DI_rm) = nan;
APD_ave(DI_rm) = nan;

% If 3 or more values are NAN then set all values to NAN on that row
checkNAN = repmat(sum(isnan(APD_ave),2) >= 3,[1 size(APD_ave,2)]);
APD_ave(checkNAN) = nan;
DI_ave(checkNAN) = nan;

end