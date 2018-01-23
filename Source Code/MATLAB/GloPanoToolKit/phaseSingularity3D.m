%% Phase Singularity Calculation%%
function [nt] = phaseSingularity3D(phase)
%
% Description:This function calculates the locations of phase singularities
% in phase data.
%
% Inputs:
% phase = This input is a single time point of phase data
%
% Outputs:
% nt = topological charge as calculated by Bray et al. 2001. The values of
% nt are normalized between -1 and 1. Values of 1 indicate CCW rotation and
% values of -1 indicate CW rotation.
%
% Author: Christopher Gloschat
% Date: January 6, 2016
%
% References:
% [1] M. A. Bray, S. F. Lin, R. R. Aliev, B. J. Roth, and J. P. Wiksow Jr.,
% "Experimental and theoretical Analysis of Phase Singularity Dynamics in 
% Cardiac Tissue," J. Cardiovasc. Electrophysiol., vol. 12, pp. 716-722,
% 2001.
%
%% CODE %%

% Calculate partial derivative in the x-axis
kx = phase(:,2:end)-phase(:,1:end-1);
kx = [kx phase(:,end)];
kbx = zeros(size(kx,1)+2,size(kx,2)+2);
kbx(2:end-1,2:end-1) = kx;
% Adjust phase jump values to be between -pi and pi
tmp = kbx > pi;
kbx(tmp) = (2*pi-kbx(tmp))*-1;
tmp = kbx < -pi;
kbx(tmp) = 2*pi+kbx(tmp);
% Calculate partial derivative in the y-axis
ky = phase(2:end,:)-phase(1:end-1,:);
ky = [ky; phase(end,:)];
kby = zeros(size(ky,1)+2,size(ky,2)+2);
kby(2:end-1,2:end-1) = ky;
% Adjust phase jump values to be between -pi and pi
tmp = kby > pi;
kby(tmp) = (2*pi-kby(tmp))*-1;
tmp = kby < -pi;
kby(tmp) = 2*pi+kby(tmp);

indX = reshape(1:size(kbx,1)*size(kbx,2),[size(kbx,1) size(kbx,2)]);
indInsideX = indX(2:end-1,2:end-1);
indInsideX = repmat(indInsideX,[1 1 9]);

indY = reshape(1:size(kby,1)*size(kby,2),[size(kby,1) size(kby,2)]);
indInsideY = indY(2:end-1,2:end-1);
indInsideY = repmat(indInsideY,[1 1 9]);

kLocX = [-size(kbx,1)-1 -1 size(kbx,1)-1;-size(kbx,1) 0 size(kbx,1); -size(kbx,1)+1 1 size(kbx,1)+1];
kLocX = reshape(kLocX,[1 1 size(kLocX,1)*size(kLocX,2)]);
kLocX = repmat(kLocX,[size(kbx,1)-2 size(kbx,2)-2 1]);
kLocX = kLocX+indInsideX;

kLocY = [-size(kby,1)-1 -1 size(kby,1)-1;-size(kby,1) 0 size(kby,1); -size(kby,1)+1 1 size(kby,1)+1];
kLocY = reshape(kLocY,[1 1 size(kLocY,1)*size(kLocY,2)]);
kLocY = repmat(kLocY,[size(kby,1)-2 size(kby,2)-2 1]);
kLocY = kLocY+indInsideY;

COx = [-0.5 0 0.5; -1 0 1; -0.5 0 0.5];
COx = reshape(COx,[1 1 size(COx,1)*size(COx,2)]);
COx = repmat(COx,[size(kLocY,1) size(kLocY,2) 1]);
COy = [0.5 1 0.5; 0 0 0; -0.5 -1 -0.5];
COy = reshape(COy,[1 1 size(COy,1)*size(COy,2)]);
COy = repmat(COy,[size(kLocX,1) size(kLocX,2) 1]);

% Calculation of topological charge
nt = sum(kby(kLocY).*COx,3) + sum(kbx(kLocX).*COy,3);
if abs(max(nt(:))) > abs(min(nt(:)))
    nt = nt/max(nt(:));
else
    nt = nt/min(nt(:));
end

% % % figure(1)
% % % % subplot(1,2,1)
% % % imagesc(phase)
% % % caxis([-pi pi])
% % % colormap('jet')
% % % set(gca,'XTick',[],'YTick',[])
% % % axis square
% % % % subplot(1,2,2)
% % % figure(2)
% % % imagesc(nt)
% % % caxis([-1 1])
% % % colormap('jet')
% % % axis square
% % % set(gca,'XTick',[],'YTick',[])

end
