

% Description: Interpolates data projected on the elements based on
% the connectivity information of the nodes.
%
% Inputs:
%   data = cmos data that needs to be interpolated
%   centroids = centroid information
%   level = number of iterations for interpolation
%
% Outputs:
%   interpolated_data = Interpolated data
%
% Author: Kedar Aras
% Date: September 11, 2017
%
% Modification Log:
%
[sorted_centroids,index] = sortrows(centroids,3);
for i=1:size(sorted_centroids,1)
    [r,~] = ind2sub(size(dataProj),index(i,1));
    sorted_cells(i,:) = cells(r,:);
    sorted_data(i,:) = dataProj(r,:);
end
sorted_data(14500:end,:) = 0; 

