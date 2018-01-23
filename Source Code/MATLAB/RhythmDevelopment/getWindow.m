function [wind_start, wind_end] = getWindow(cmosData, Fs)
% Get Activation Time for Camera A
actTime = getActTime(cmosData, 1);
% Calculate Dominant Frequency
domFreqMap = cell2mat(calDomFreqPano(cmosData,Fs));

domFreqPeriod = 1/mean2(domFreqMap(domFreqMap ~= 0));

% Establish a window with the activation time and dominant frequency
wind_start = actTime;
wind_end = actTime + domFreqPeriod;

% Iterate through the other 3 cameras
for n=2:4
    tempActTime = getActTime(cmosData,n);
    
    % Adjust start time as necessary
    while tempActTime < wind_start
        wind_start = wind_start - domFreqPeriod/100;
    end
    
    % Get APD in window
    apdMap = getAPD(cmosData,wind_start,wind_end,Fs,0.5);
    apdMapCam = apdMap{n};    
    apd = max(apdMapCam(apdMapCam ~= NaN))/1000;
    
    % Adjust end time as necessary
    endAPD = actTime + apd;
    while endAPD > wind_end
        wind_end = endAPD + domFreqPeriod/8;
    end
end
end

%% Activation Time Function %%
function actTime = getActTime(data, n)
    % Find First Derivative and time of maxium
    if size(data{n},3) == 1
        temp2 = diff(data{m},1,2);
        %[~,max_i] = max(temp2,[],2);
    else
        temp2 = diff(data{n},1,3); % first derivative
        %[~,max_i] = max(temp2,[],3); % find location of max derivative
    end
    sortedActTimeMap = sort(temp2,3,'descend');
    maxActTimeMap = sortedActTimeMap(:,:,1);
    actTime = min(maxActTimeMap(maxActTimeMap ~= NaN));
end

%% Dominant Frequency Function %%
function maxf = calDomFreqPano(data,Fs)
maxf = cell(1,length(data));
for k = 1:length(data)
    %% Window Data with Tukey Window to Minimize Edge Effects
    if size(data{k},3) == 1
        w_m = tukeywin(size(data{k},2),.05);
        win = repmat(w_m',[size(data{k},1) 1]);
    else
        w_m = tukeywin(size(data{k},3),.05);
        win = repmat(permute(w_m,[3 2 1]),[size(data{k},1),size(data{k},2)]);
    end
    data{k} = data{k}.*win;
    %% Find single-sided power spectrum of data
    if size(data{k},3) == 1
        m = size(data{k},2);           % Window length
        n = pow2(nextpow2(m));      % Transform Length
        y = fft(data{k},n,2);          % DFT of signal
    else
        m = size(data{k},3);           % Window length
        n = pow2(nextpow2(m));      % Transform Length
        y = fft(data{k},n,3);          % DFT of signal
    end
    f = Fs/2*linspace(0,1,n/2+1);   % Frequency range
    p = y.*conj(y)/n;               % Power of the DFT
    if size(data{k},3) == 1
        p_s = 2*abs(p(:,1:n/2+1));      % Single-sided power
        p_s(:,1) = [];                % Remove DC component
    else
        p_s = 2*abs(p(:,:,1:n/2+1));    % Single-sided power
        p_s(:,:,1) = [];                % Remove DC component
    end
    f(1) = [];                      % Remove DC
    
    %% Find Dominant Frequency
    if size(data{k},3) == 1
        [val,ind] = max(p_s,[],2);
        maxf{k} = f(ind).*isfinite(val)';
    else
        [val,ind] = max(p_s,[],3);
        maxf{k} = f(ind).*isfinite(val);
    end
    
    % % %     if size(data{k},3) ~= 1
    % % %         figure; imagesc(maxf{k}); colormap(cmap);colorbar
    % % %         caxis([mean2(maxf{k})*.1 mean2(maxf{k})*2])
    % % %     end
end
end

%% Action Potential Duration Function %%
function [apdMap] = getAPD(data,start,endp,Fs,percent)
%% Create initial variables
start=round(start*Fs)+1;
endp=round(endp*Fs)+1;
% Code not used in current version %
apdMap = cell(5,1);

for n = 1:4
    apd_data = data{n}(:,:,start:endp);        % window signal
    apd_data = normalize_data(apd_data); %re-normalize windowed data
    
    %%Determining activation time point
    % Find First Derivative and its index of maximum
    apd_data2 = diff(apd_data,1,3); % first derivative
    [~,max_i] = max(apd_data2,[],3); % find location of max derivative
    
    
    %%Find location of repolarization
    %%Find maximum of the signal and its index
    [~,maxValI] = max(apd_data,[],3);
    
    %locs is a temporary holding place
    locs = nan(size(apd_data,1),size(apd_data,2));
    
    %Define the baseline value you want to go down to
    requiredVal = 1.0 - percent;
    
    %%for each pixel
    for i = 1:size(apd_data,1)
        for j = 1:size(apd_data,2)
            %%starting from the peak of the signal, loop until we reach baseline
            for k = maxValI(i,j):size(apd_data,3)
                if apd_data(i,j,k) <= requiredVal
                    locs(i,j) = k; %Save the index when the baseline is reached
                    %this is the repolarizatin time point
                    break;
                end
            end
        end
    end
    
    %%account for different sampling frequencies
    unitFix = 1000.0 / Fs;
    
    % Calculate Action Potential Duration
    apd = minus(locs,max_i);
    apdMap{n} = apd * unitFix;
    apdMap{n}(apdMap{n} <= 0) = 0;
end
end