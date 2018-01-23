%% Average triangle surface area %%
function [meanTSA,stdTSA] = aveTriSurfArea(cells,pts)
% Preallocate variables for speed
edgeDists = zeros(size(cells,1),size(cells,2));
A = zeros(size(cells,1),1);
for n = 1:size(edgeDists,1)
   % Grab the vertices of a cell
   xyz =  pts(cells(n,:),:);
   % Calculate the length of each edge
   edgeDists(n,1) = sqrt(sum((xyz(2,:)-xyz(1,:)).^2));
   edgeDists(n,2) = sqrt(sum((xyz(2,:)-xyz(3,:)).^2));
   edgeDists(n,3) = sqrt(sum((xyz(1,:)-xyz(3,:)).^2));   
   % Use Heron's Formula to calculate the area using the edge lengths
   S = sum(edgeDists(n,:))/2;
   A(n) = sqrt(S*(S-edgeDists(n,1))*(S-edgeDists(n,2))*(S-edgeDists(n,3)));
end
% Save out the average triangle surface area and the standard deviation
meanTSA = mean(A);
stdTSA = std(A);
end