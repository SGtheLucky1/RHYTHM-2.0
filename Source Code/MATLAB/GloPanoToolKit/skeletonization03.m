%% Skeletonization Algorithm %%
function [skel] = skeletonization03(GRID,k)
% Make copy of image
skel = GRID;
skelCE = skel;
% Get height of image
imY = size(skel,1);
% Create 8 connectivity kernel
kernel8 = [-imY-1 -imY -imY+1 -1 0 1 imY-1 imY imY+1];
% Create 4 connectivity kernel
kernel4 = [-imY -1 0 1 +imY];
if k == 8
    kernel = kernel8;
elseif k == 4
    kernel = kernel4;
end

% Variable for tracking non-simple pixels
nonSimple = [];
% Variabe for tracking curve end pixels
curveEnd = [];
% Variable for storing border pixel indices
pixBorder = [];
% Variable for storing removed border pixel indices
pixRemove = [];
% Pixels to check
check = find(skel);
% Number of pixels to check
numCheck = length(check);
while numCheck > (length(pixRemove)+length(curveEnd)+length(nonSimple))
    % Reset the curve end and non simple variables
    curveEnd = [];
    nonSimple = [];
    % Identify the indices for the neighbors of each object pixel
    pixNbrs = repmat(check,[1 size(kernel,2)])+repmat(kernel,[size(check,1) 1]);
    % Grab the neighbors image value (0 or 1)
    findBorder = skel(pixNbrs);
    % Sum the neighbors values
    findBorder = sum(findBorder,2);
    % Border pixels will have at least one background neighbor
    findBorder = findBorder < k+1;
    % Grab image indices
    pixBorder = [pixBorder; check(findBorder)];
    % Remove these pixels from the check variable
    check(findBorder) = [];
    while ~isempty(pixBorder)
        % Grab the start value assessing a connected component
        seed = pixBorder(1);
        pixBorder(1) = [];
        while ~isempty(seed)
            % Identify the connected neighbors
            seedNbrsInd = seed(1)+kernel;
            % Identify the neighbor values
            seedNbrsVal = skelCE(seedNbrsInd)==1;
            % Grab object neighbors
            seedAdd = seedNbrsInd(seedNbrsVal);
            % Remove them from the border pixel variable
            if ~isempty(pixBorder)
                seedRemove = repmat(seedAdd,[size(pixBorder,1) 1]) == ...
                    repmat(pixBorder,[1 length(seedAdd)]);
                seedRemove = seedRemove.*repmat((1:length(pixBorder))',[1 length(seedAdd)]);
                % Sum to identify the matched index for each pixel, non matches
                % (equal to 0) mean it is not yet a border pixel
                seedRemove = sum(seedRemove);
                pixBorder(seedRemove(seedRemove~=0)) = [];
                % Add the object border neighbors to the seed queue
                seed = [seed seedAdd(seedRemove~=0)];
            end
            % Take first pixel in queue and grab its neighbor values
            CE = skel(seed(1)+kernel);
            % Check if it has a single neighbor
            if sum(CE) == 2
                % If it does place it in a tracking variable
                curveEnd = [curveEnd; seed(1)];
                % Remove from curve end skel variable
                skelCE(seed(1)) = 0;
                % Remove it from the queue
                seed = seed(2:end);
            else
                % If not, check how connectivity changes with removal
                % Grab neighbhorhood around pixel
                rm = skel(seed(1)+kernel8);
                % Check object connectivity before removal
                obj1 = gloConnCount(rm,1);
                % Check background connectivity before removal
                bkg1 = gloConnCount(rm,0);
                rm(5) = 0;
                % Check object connectivity after removal
                obj2 = gloConnCount(rm,1);
                % Check background connectivity after removal
                bkg2 = gloConnCount(rm,0);
                % If the number of foreground and background components does
                % not change with removal it is simple and can be removed
                if obj1 == obj2 && bkg1 == bkg2
                    % Remove the pixel from the image
                    skel(seed(1)) = 0;
                    % Remove from curve end check image
                    skelCE(seed(1)) = 0;
                    % Place index in a removed pixel variable
                    pixRemove = [pixRemove seed(1)];
                    % Remove the pixel from the queue
                    seed = seed(2:end);
                else
                    % Place index in a non simple pixel variable
                    nonSimple=[nonSimple; seed(1)];
                    % Remove the pixel from the queue
                    seed = seed(2:end);
                end
            end
            disp(length(seed))
        end
    end
    % Place the curve end and non simple pixels back in the check variable
    check = [check; curveEnd; nonSimple];
end

end