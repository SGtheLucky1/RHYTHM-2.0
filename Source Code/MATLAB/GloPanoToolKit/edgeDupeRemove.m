function [dupes] = edgeDupeRemove(edges)
% Description: The simplification of a mesh is done by removing the edges
% of triangles that will least change the mesh geometry. When a edge is
% removed, duplicate edges are generated. This alogrithm finds and removes
% them.
% 
% Input:
% edges = the edges of a triangular mesh
%
% Outputs:
% dupes = the indices of the duplicates to be removed
%
% Author: Christopher Gloschat
% Date: April 19, 2017
%
% Modification Log:
%
%
%% CODE %%
% Sort the edges to match edges with different ordered vertices
edges = sort(edges,2);
% Create empty variable for duplicate indices
dupes = [];
% Establish counter
n = 1;
% Step through each edge finding and removing duplicates
while n < size(edges,1)
    % Create variable for checking for duplicates
    check = repmat(edges(n,:),[size(edges,1) 1]);
    % Check against current edges variable
    check = sum(check == edges,2) == 2;
    % If there are duplicates, remove all but the first one
    if sum(check)~=0
        % Grab the indices
        check = find(check);
        if size(check,2) > 1
            check = check';
        end
        % Discard all but the first one
        dupes = [dupes;check(2:end)];
% % %         edges(check,:) = [];
    end
    % Increment the counter
    n = n+1;
end

end