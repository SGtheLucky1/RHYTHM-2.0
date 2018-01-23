%% SNR Map %%
close all
clear all
clc
% Load the data
load('06_projectedData')
handles = a;
load('06_SNRpointsAnt')
M = a;
rawData = x;
clear a x
% Grab the signals and the pacing spike
% data = rawData;
data = handles.dataProj;
pacedSignal = handles.ecg;
nRate = handles.nRate;
centroids = handles.centroids;

noSig = find(sum(data,2) == 0);
sig = (1:size(data,1))';
sig(noSig) = [];


% Identify time when the paced data is equal to one
up = pacedSignal == 0;
down = pacedSignal ~= 0;

% Identify the indices to be used in identifying the begining of pacing
i = 2:length(pacedSignal)-1;

% Identify the indices associated with the upstroke of pacing
pacedInd = up(i).*up(i+1).*down(i-1);
pacedInd = unique(find(pacedInd));
% Convert to optical data equivalent index
pacedInd = round(pacedInd/nRate);

% Pacing rates
pacedRates = pacedInd(2:end)-pacedInd(1:end-1);
CL = round(mean(pacedRates));
sections = zeros(size(data,1),CL,length(pacedRates));
SNR = zeros(size(data,1),length(pacedRates));
% Isolate periods for analysis
start = pacedInd(1:end-1) - 50;
fin = start+CL-1;

for n = 1:length(start)
   % Find linear drift
   frontAve = mean(data(:,start(n):start(n)+19),2);
   backAve = mean(data(:,fin(n)-19:fin(n)),2);
   segment = fin(n)-start(n);
   slope = (frontAve-backAve)/segment;
   slope = repmat(0:segment,[size(slope,1) 1]).*repmat(slope,[1 segment+1]);
   sections(:,:,n) = data(:,start(n):fin(n))+slope;
   signal = max(sections(:,:,n),[],2)-min(sections(:,:,n),[],2);
   noise = max(sections(:,1:50,n),[],2)-min(sections(:,1:50,n),[],2);
   SNR(:,n) = signal./noise;
    
end

snr = mean(SNR,2);

figure
trisurf(handles.cells(sig,:),handles.pts(:,1),...
handles.pts(:,2),handles.pts(:,3),snr(sig,:),'LineStyle','none');
colormap('jet')
caxis([0 max(snr)])
hold on
% Plot triangles where no data was assigned
trisurf(handles.cells(noSig,:),handles.pts(:,1),handles.pts(:,2),...
    handles.pts(:,3),'FaceColor','none',...
    'EdgeColor',[0.5 0.5 0.5]);

% Find the points in the centroids variable
repSignals = reshape(M{5}',[1 3 3]);
repSignals = repmat(repSignals,[size(centroids,1) 1 1]);
repSignals = sum(repmat(centroids,[1 1 3])-repSignals,2);
repSignals = reshape(repSignals,[size(repSignals,1) size(repSignals,3)]);
[~,I] = min(abs(repSignals));
% Visualize
% % % colorF = [1 0 1;
% % %     128/255 128/255 128/255;
% % %     0 0 0];
colorF = [1 0 1;
    128/255 128/255 128/255;
    0 0 0];

for n = 1:3
scatter3(centroids(I(n),1),centroids(I(n),2),centroids(I(n),3),...
    'SizeData',128,'MarkerFaceColor',colorF(n,:),'MarkerEdgeColor',[1 1 1],'LineWidth',2)
end
axis equal
view(225,0)
set(gca,'Visible','off')

% Plot OAPs
figure
for n = 1:3
    subplot(3,1,n)
    plot(sections(I(n),:,round(length(pacedRates)/2)),'Color',colorF(n,:),'LineWidth',2)
    set(gca,'XLim',[0 CL-1])
end





% Calculate the SNR as the ratio of the base 