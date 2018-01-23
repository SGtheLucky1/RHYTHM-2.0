function [snr] = aveSNR(data,pace,Fs,nRate)

% Create a time variable
time = 0:1/Fs:size(data,3)/Fs;
    % Find all paced time points
    paceInd = find(pace == 0);
    % Find the beginning of each paced beat
    paceStart = paceInd(2:end)-paceInd(1:end-1);
    paceStart = [1;find(paceStart~=1)+1];
    paceStart = paceInd(paceStart);
    % Convert to the same timing as the optical cameras
    paceStart = round(paceStart/nRate);
    % Calculate the cycle length of the paced beats
    cycleLength = round(mean((time(paceStart(2:end))-time(paceStart(1:end-1)))*Fs));
    % If the first beat is within 50 ms of the beginning remove
    if time(paceStart(1)) < 0.05
        paceStart = paceStart(2:end);
    end
    % If the final beat is within 50 ms of the end remove
    if time(paceStart(end))+cycleLength/Fs > time(end)
        paceStart(end) = [];
    end
    
    % Grab time windows for calculations
    winData = repmat(paceStart,[1 cycleLength]);
    winData = winData+repmat(-50:cycleLength-51,[size(winData,1) 1]);
    
    noise = zeros(size(data,1),size(data,2),size(winData,1));
    signal = zeros(size(data,1),size(data,2),size(winData,1));
    snrN = zeros(size(data,1),size(data,2),size(winData,1));
    
    for m = 1:size(winData,1)
        disp(m)
        winDataN = reshape(winData(m,:),[1 1 size(winData,2)]);
        winDataN = repmat(winDataN,[size(data,1),size(data,2),1]);
        ind = reshape(1:size(data,1)*size(data,2),[size(data,1),size(data,2)]);
        ind = repmat(ind,[1,1,size(winDataN,3)]);
        winDataN = ind+winDataN*(size(data,1)*size(data,2));
        % Calculate the baseline noise
        winDataN = data(winDataN);
        noise(:,:,m) = max(winDataN(:,:,1:30),[],3)-min(winDataN(:,:,1:30),[],3);
        signal(:,:,m) = max(winDataN,[],3)-min(winDataN,[],3);
        snrN(:,:,m) = signal(:,:,m)./noise(:,:,m);
        
    end
    % Save out the calcualted snr values
    snr = mean(snrN,3);

end


% 1) Create 3D index  matrices for isolating the part of each signal for
% both noise and signal
% 2) Calculate signal, noise, then snr
% 3) Play with threshold removal until most noisy border pixels are removed