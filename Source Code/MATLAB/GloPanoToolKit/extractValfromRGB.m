%% Extract rgb data from ply to colomap values %%

% Load point cloud
ptCloud21 = pcread('HausdorffPoissonMesh_02_01.ply');

% Convert to values between 0 and 1
poisson21 = ptCloud21.Color;
poisson21 = double(poisson21)/255;

% Create dense colormap
cmap21 = colormap(jet(1000));

% Create matching data values
cmapInd21 = linspace(0,1.04,1000);
match = zeros(1,size(poisson21,1));
for n = 1:size(poisson21,1)
    if sum(poisson21(n,:) == 1) == 0
        if sum(poisson21(n,:) == 0) < 2
            tmp = find(poisson21(n,:)~= 0);
            if poisson21(tmp(1)) > poisson21(tmp(2))
                poisson21(n,tmp(1)) = 1;
            else
                poisson21(n,tmp(2)) = 1;
            end
        end
    end
    %Identify matches for each color
    Rup = cmap21(:,1) < poisson21(n,1)+0.0020;
    Rdn = cmap21(:,1) > poisson21(n,1)-0.0020;
    R = Rup.*Rdn;
    Gup = cmap21(:,2) < poisson21(n,2)+0.0020;
    Gdn = cmap21(:,2) > poisson21(n,2)-0.0020;
    G = Gup.*Gdn;
    Bup = cmap21(:,3) < poisson21(n,3)+0.0020;
    Bdn = cmap21(:,3) > poisson21(n,3)-0.0020;
    B = Bup.*Bdn;
    % Find the color that matches all three
    match(n) = find(R.*G.*B);   
end




% % % match = repmat(reshape(cmap21',[1 size(cmap21,2) size(cmap21,1)]),[size(poisson21,1) 1 1]);
% % % match = sum((match-repmat(poisson21,[1 1 size(cmap21,1)])).^2,2);
% % % [mV,mI] = min(match,[],3);
% % % 
% % % m21val = cmapInd21(mI);
