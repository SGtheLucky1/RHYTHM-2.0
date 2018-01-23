% Restitution Maps

% Max slope maps
figure
map = zeros(100,100,4);
for n = 1:4
tmp = nan(100,100);
tmp(dataInd{n}(:,1)) = maxSlope{n};
map(:,:,n) = tmp;
subplot(2,2,n)
imagesc(map(:,:,n))
colormap('jet')
caxis([0 2.5])
end


% Ratio Maps
APD_change = cell(1,4);
DI_change = cell(1,4);
Rk = cell(1,4);
for n = 1:4
APD_change{n} = max(APD_ave{n},[],2)-min(APD_ave{n},[],2);
DI_change{n} = max(DI_ave{n},[],2)-min(DI_ave{n},[],2);
Rk{n} = APD_change{n}./DI_change{n};
end

figure
for n = 1:4
tmp = nan(100,100);
tmp(dataInd{n}(:,1)) = Rk{n};
map(:,:,n) = tmp;
subplot(2,2,n)
imagesc(map(:,:,n))
colormap('jet')
end