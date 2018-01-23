
up = find(a.ecg == 0);
% Each spike is actually a step with multiple values, subtract to find the
% jumps from one spike to the next
front = up(2:end)-up(1:end-1);
% Add 1 to the subtracted values to get the values of the next steps
up = up([1; find(front~=1)+1]);
% Pacing rate
% % % CL = round(mean((up(2:end)-up(1:end-1))/nRate));
% Pacing times at 1kHz timing
pacedInd = round(up/a.nRate);

% Group paced beats based on cycle length
pacedCL = pacedInd(2:end)-pacedInd(1:end-1);
pacedCL = [pacedCL(1);pacedCL];
pacedTransInd = [1; find((pacedCL(2:end)-pacedCL(1:end-1)) == -30)+1];
pacedTrans = pacedInd(pacedTransInd);


color = 'bmgkcr';
figure
hold on
for n = 1:length(data{3})
    if n == 1
        ind1 = 1;
        ind2 = pacedTrans(n+1);
    elseif n == length(data{3})
        ind2 = size(a.cmosData{3},3);
    else
        ind2 = pacedTrans(n+1);
    end
    plot((ind1:ind2)/a.Fs,squeeze(a.cmosData{3}(50,30,ind1:ind2)),'Color',color(n),'LineWidth',2)
    ind1 = ind2;
end






