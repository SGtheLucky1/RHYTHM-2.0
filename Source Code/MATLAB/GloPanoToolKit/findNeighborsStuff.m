function [nbrCells,nbrNorms,nbrVerts] = findNeighborsStuff(ind,edges,...
faces,norms,nbrCells,nbrNorms,nbrVerts)
% Description: Finds the cells, their normals, and their associated 
% vertices that are adjacent to each edge. The index provided by the user
% specifices for which edges this is done.
% 
% Input:
% ind = A column vector of the indices of the edges for which this 
%       information will be calculated
% edges = A matrix with all the edges
% faces = A matrix with all faces
% norms = A matrix with the normal at each face
% nbrCells = A cell with the faces adjacent to each edge
% nbrNorms = A cell with the normals of these faces
% nbrVerts = A cell with the vertices of these faces
%
% Output:
% nbrCells = Updated version
% nbrNorms = Updated version
% nbrVerts = Updated version
%
% Author: Christopher Gloschat
% Date: April 19, 2017
%
% Modification Log:
%
%
%% Code %%

for n = 1:length(ind)
    % Find the neighboring cells
    tmp = reshape(edges(ind(n),:),[1 1 2]);
    tmp = repmat(tmp,[size(faces,1) 3]);
    nbrCells{ind(n)} = find(sum(sum(tmp == repmat(faces,[1 1 2]),2),3)~=0);
    nbrNorms{ind(n)} = norms(nbrCells{ind(n)},:);
    % Find the connected vertices
    nbrVerts{ind(n)} = unique(faces(nbrCells{ind(n)},:));
    tmp = repmat(nbrVerts{ind(n)},[1 2]) == repmat(edges(ind(n),:),...
        [size(nbrVerts{ind(n)},1) 1]);
    %     tmp = find(sum(tmp,2)~=0);
    nbrVerts{ind(n)}(sum(tmp,2)~=0) = [];
end


end