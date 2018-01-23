function [kernel] = makeKernel2D(im,num)
% Description:
%
% Input:
%
%
% Output:
%
%
% Author: Christopher Gloschat
% Date: January 5, 2017
%
% Modfication Log:
%
%
%% Code %%

% Get the number of rows for figuring out index values
imY = size(im,1);
% Create variable for vertical component
vert = ((num-1)/2:-1:-(num-1)/2)';
% Create empty variable for kernel
kernel = zeros(length(vert),length(vert));
% Create kernel
for n = 1:length(vert)
    if n == round(num/2)
        kernel(:,n) = vert;
    elseif n < round(num/2)
        kernel(:,n) = -imY*(n+floor(num/2)-((n-1)*2+1))-vert;
    elseif n > round(num/2)
        kernel(:,n) = imY*(n-round(num/2))-vert;
    end
end
        









end