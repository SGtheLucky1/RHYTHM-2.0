function interp_data = interpolate_data(data,centroids,method)
% Description: Interpolates data projected on the elements based on
% the connectivity information of the nodes.
%
% Inputs:
%   data = cmos data that needs to be interpolated
%   centroids = centroid information
%
% Outputs:
%   interp_data = Interpolated data
%
% Author: Kedar Aras
% Date: September 11, 2017
%
% Modification Log:
%

%Step 0: Initialize the variables as needed
interp_data = data;

%Step 1: Separate centroids with known data from those with missing data
%that needs to be interpolated
num_frames = size(data,2);
global_known_centroids = cell(1,num_frames);
global_known_data = cell(1,num_frames);
global_unknown_centroids = cell(1,num_frames);

for i=1:num_frames
    %find centroids with data available
    known_c = find(data(:,i) > 0);
    [known_c,~] = ind2sub(size(data),known_c(:));
    known_c = unique(known_c);
    known_centroids = centroids(known_c,:);
    [unique_known_centroids, rows] = unique(known_centroids, 'rows', 'stable');
    dup_rows = (setdiff(1:size(known_centroids,1), rows))';
    known_data = data(rows,i);
    
    %find centroids with data unavailable
    data(dup_rows,i) = 0;
    unknown_c = find(data(:,i) == 0);
    [unknown_c,~] = ind2sub(size(data),unknown_c(:));
    unknown_c = unique(unknown_c);
    unknown_centroids = centroids(unknown_c,:);
    
    min_candidate = size(unknown_c,1);
    if i==1
        frame_index = 1;
        min_val = min_candidate;
    else
        if min_candidate < min_val
            min_val = min_candidate;
            frame_index = i;
        end
    end
        
    
    %book keeping
    global_known_centroids{1,i} = [unique_known_centroids,rows];
    global_known_data{1,i} = [known_data,rows];
    global_unknown_centroids{1,i} = [unknown_centroids,unknown_c];
    
end


%Step 2: Interpolate the missing data
global_interpolated_data = cell(1,num_frames);
interpolation_function = [];
fn_known_centroids = global_known_centroids{1,frame_index};
fn_known_data = global_known_data{1,frame_index};
for j=1:num_frames
    known_centroids = global_known_centroids{1,j};
    known_data = global_known_data{1,j};
    unknown_centroids = global_unknown_centroids{1,j};
    if isempty(interpolation_function)
        interpolation_function = scatteredInterpolant(fn_known_centroids(:,1),fn_known_centroids(:,2),fn_known_centroids(:,3), fn_known_data(:,1),method);
    end
    interpolated_data = interpolation_function(unknown_centroids(:,1),unknown_centroids(:,2),unknown_centroids(:,3)); 
    global_interpolated_data{1,j} = interpolated_data;
end

%Step 3: Consolidate the interpolated data with original data
for k=1:num_frames
    interpolated_data = global_interpolated_data{1,k};
    unknown_centroids = global_unknown_centroids{1,k};
    num_points = size(interpolated_data,1); 
    for m=1:num_points
        index = unknown_centroids(m,4);
        interp_data(index,k) = interpolated_data(m,1);
    end
end


end

