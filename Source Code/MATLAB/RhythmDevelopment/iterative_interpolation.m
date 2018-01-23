function interpolated_data = iterative_interpolation(data, centroids, level )

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

for i=1:level
    if i == 1
        interpolated_data = interpolate_data(data,centroids,'nearest');
    else
        interpolated_data = interpolate_data(data,centroids,'linear');
    end
    zero_values = find(interpolated_data == 0);
    size(zero_values,1)
    if isempty(zero_values)
        break;
    end 
end


end

