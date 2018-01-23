% Remove Triangle Islands
function [track] = removeIslands(cell)

% Save cells to a storage variable
check = cell;
% Create a variable for tracking indices and connected components
track = [1,1];
ind = (1:size(cell,1))';
% Create a variable for storing already checked verts
alreadyVerts = [];
% Remove checked cells
check(1,:) = [];
ind(1) = [];

iter = 1;
while ~isempty(check)
    % Create a variable for vertices who need to be check for adjoining cells
    verts = check(1,:)';
% % %     i = 0;
    while ~isempty(verts)
% % %         i = i+1
        if verts(1) == 2812
            disp(i)
        end
        % Check the vertices of the active triangle for use in other cells
        tmp = check == verts(1);
        tmp = tmp.*repmat((1:size(check,1))',[1 size(tmp,2)]);
        tmp = unique(tmp);
        tmp = tmp(2:end);
        % Add conjoined cells to assigned & remove from check
        track = [track; [ind(tmp) iter*ones(size(tmp,1),size(tmp,2))]];
        % Add their vertices to the list of vertices to check and remove
        add = check(tmp,:);
        check(tmp,:) = [];
        ind(tmp,:) = [];
        tmp = unique(add(add~=verts(1)));
        tmp = reshape(tmp,[length(tmp) 1]);
        verts = [verts; tmp];
        % Remove checked vertix
        alreadyVerts = [alreadyVerts; verts(1)];
        verts(1) = [];
        % Ensure uniqueness of all remaining vertices
        verts = unique(verts);
        if ~isempty(verts)
            % Remove any already checked verts
            tmp = repmat(alreadyVerts',[size(verts,1) 1])==repmat(verts,[1 size(alreadyVerts,1)]);
            tmp = sum(tmp,2);
            if sum(tmp) ~= 0
                verts = verts(~tmp);
            end
        end
    end
    % Increment connected component iterator
    iter = iter+1;
end

[~,newOrder] = sort(track(:,1),'ascend');
track = track(newOrder,:);




end