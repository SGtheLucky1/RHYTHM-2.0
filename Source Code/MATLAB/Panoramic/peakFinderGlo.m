function [peakInds] = peakFinderGlo(data,CL,bind,eind)
% Description: A variant of peakfinder, specifically for paced data, that 
% grabs the largest peak and then looks for other peaks at cycle length
% intervals before and after the peak.
%
% Input:
%   data = data in a signals (rows) by time (columns) format
%   CL = cycle length
%   bind = number of indices to grab before activation
%   eind = number of indices to grab after activation
% 
% Output:
%   peakInds = location of dVdt peaks indicating activation
%
% Author: Christopher Gloschat
% Date: August 31, 2016
%
% Modification Log:
%
%
%% Code %%

% find first derivative
dVdt = data(2:end,:) - data(1:end-1,:);
% find max peak of derivative
[~,biggestPeak] = max(dVdt);

% find the peaks before
if biggestPeak > CL+bind
    checkB = biggestPeak:-CL:0;
else
    checkB = biggestPeak;
end
% find the peaks after
checkA = [];
if biggestPeak < length(data)-eind
    checkA = biggestPeak:CL:length(data);
end
% combine peaks
check = [flip(checkB) checkA(2:end)];
% verify existence of peaks
peakInds = zeros(1,length(check));
for n = 1:length(check)
    % window to check
    ind = check(n)-5:check(n)+5;
    % check window[
    [~,maxInd] = max(dVdt(ind));
    % update peak value
    peakInds(n) = ind(maxInd);
end
% check for peak to close to the end
if peakInds(end)+eind > length(data)
    peakInds(end) = [];
end
% check for peak to close to the front
if peakInds(1)+bind < 0
    peakInds(1) = [];
end


end