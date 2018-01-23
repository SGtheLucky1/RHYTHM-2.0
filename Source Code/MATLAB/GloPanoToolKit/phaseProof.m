% Do a representative data set

n = 3;

[p,t,h] = phaseMapAll(cmosData{n},0,(size(cmosData{n},3)-1)/a.Fs,a.Fs);

% Identify time when the paced data is equal to zero
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

dhdt = real(h(2:end,:))-real(h(1:end-1,:));
act_R = nan(length(pacedInd),size(h,2));
nonan = ~isnan(real(h(1,:)));
for m = 1:length(pacedInd)
    [~,act_R(m,nonan)] = max(dhdt(pacedInd(m):pacedInd(m)+80,nonan));
    act_R(m,nonan) = act_R(m,nonan)+pacedInd(m);
end

act_C = 1:size(h,2);
act_C(~nonan) = nan;
act_C = repmat(act_C,[length(pacedInd) 1]);

act_I = sub2ind([size(h,1),size(h,2)],act_R,act_C);
act_Vx = nan(size(act_I,1),size(act_I,2));
act_Vx(~isnan(act_I)) = real(h(act_I(~isnan(act_I))));
act_Vy = nan(size(act_I,1),size(act_I,2));
act_Vy(~isnan(act_I)) = imag(h(act_I(~isnan(act_I))));


% % % h_I = sub2ind([size(h,1),size(h,2)],act_C,act_R);
% % % h_V = nan(size(h_I,1),size(h_I,2));
% % % h_V(~isnan(h_I)) = h(h_I(~isnan(h_I)));
% Visualize single signal
n = 4950;
figure,
% h2 = fliplr(h');
plot(real(h(:,n)),imag(h(:,n)))
hold on
scatter(0,0,'ro','SizeData',16)
xlabel('V_m(t)','FontSize',24)
ylabel('H(V_m(t))','FontSize',24)
set(gca,'FontSize',18)

scatter(act_Vx(:,n),act_Vy(:,n),'mo')




