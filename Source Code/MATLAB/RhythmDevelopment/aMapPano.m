function [actMap1] = aMapPano(data,stat,endp,Fs)
%% aMap is the central function for creating conduction velocity maps
% [actMap1] = aMap(data,stat,endp,Fs,bg) calculates the activation map
% for a single action potential upstroke.

% INPUTS
% data = cmos data (voltage, calcium, etc.) from the micam ultima system.
% stat = start of analysis (in msec)
% endp = end of analysis (in msec)
% Fs = sampling frequency
% bg = black and white background image from the CMOS camera.  This is a
%      100X100 pixel image from the micam ultima system. bg is stored in 
%      the handles structure handles.bg.%
%
% OUTPUT
% actMap1 = activation map
%
% METHOD
% An activation map is calculated by finding the time of the maximum 
% derivative of each pixel in the specified time-windowed data.
%
% REFERENCES
%
% ADDITIONAL NOTES
%
% RELEASE VERSION 1.0.1
%
% AUTHOR: Qing Lou, Jacob Laughner (jacoblaughner@gmail.com)
%
% MAINTAINED BY: Christopher Gloschat - (cgloschat@gmail.com) - [Jan. 2015 - Present] 
%
% MODIFICATION LOG:
%
% Jan. 26, 2015 - The input cmap was added to input the colormap and code
% was added at the end of the function to set the colormap to the user
% determined values. In this case the most immediate purpose is to
% facilitate inversion of the default colormap.
%
%

%% Code
% Create initial variables
stat=round(stat*Fs)+1;
endp=round(endp*Fs)+1;
% Code not used in current version %
actMap1 = cell(length(data),1);
% % mask2 = zeros(size(data,1),size(data,2));

for n = 1:length(data)
    % identify channels that have been zero-ed out due to noise
    if size(data{n},3) == 1
        temp = data{n}(:,stat:endp);       % Windowed signal
        temp = normalize_data(temp);    % Re-normalize data in case of drift
        mask = max(temp,[],2) > 0;      % Generate mask
    else
        temp = data{n}(:,:,stat:endp);     % Windowed signal
        temp = normalize_data(temp);    % Re-normalize data in case of drift
        mask = max(temp,[],3) > 0;      % Generate mask
    end
    
    % Find First Derivative and time of maxium
    if size(data{n},3) == 1
        temp2 = diff(temp,1,2);
        [~,max_i] = max(temp2,[],2);
    else
        temp2 = diff(temp,1,3); % first derivative
        [~,max_i] = max(temp2,[],3); % find location of max derivative
    end
    
    % Activation Map Matrix
    actMap1{n} = max_i.*mask;
    actMap1{n}(actMap1{n} == 0) = 0;
    offset1 = min(min(actMap1{n}));
    actMap1{n} = actMap1{n} - offset1*ones(size(actMap1{n},1),size(actMap1{n},2));
    actMap1{n} = actMap1{n}/Fs*1000; %% time in ms
end

% identify channels that have been zero-ed out due to noise

% % % % 
% % % % % Find the first derivative of the 3D data
% % % % dataGeoDVDT = dataGeo(:,2:end)-dataGeo(:,1:end-1);
% % % % % Find the index value of the max dVdt for each signal
% % % % [~,maxGeo_i] = max(dataGeoDVDT,[],2);
% % % % % Set all noisy signals to NAN
% % % % maxGeo_i(maxGeo_i == 1) = nan;
% % % % % Remove offset by setting smallest activation time to zero
% % % % maxGeo_i = maxGeo_i-min(maxGeo_i);
% % % % actMapGeo = maxGeo_i/Fs*1000;


end




