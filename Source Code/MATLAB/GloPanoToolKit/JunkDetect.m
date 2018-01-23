%% Automatic Junction Detection %%
% % % close all
% Read in calibration image
I = imread('Cube01.tiff');
% % % figure
% % % image(I)

%% Isolate block from background %%
% Use a threshold to grab the block
TH = I(:,:,1) > 25;
% Identify largest connected component to separate from noise in bkgrd
CC = bwconncomp(TH,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
BW = zeros(size(TH,1),size(TH,2));
BW(CC.PixelIdxList{idx}) = 1;
% Erode and dilate to eliminate spurs
se = strel('square',7);
BW2 = imerode(BW,se);
BW3 = imdilate(BW2,se);
% Dilate and erode to close any gaps
BW4 = imdilate(BW3,se);
BW5 = imerode(BW4,se);
% Isolate this portion of the image
I2 = I;
I2(~repmat(BW5,[1 1 3])) = 0;
% Convert image to a grayscale image
BW = rgb2gray(I2);
% % % BW2 = bwmorph(BW,'skel',2);
% Identify grid
BLOCK = BW > 75;
% se = strel('square',9);
BLOCK = (BLOCK == 0).*imerode(BW5,se);
CC = bwconncomp(BLOCK,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
GRID = zeros(size(TH,1),size(TH,2));
GRID(CC.PixelIdxList{idx}) = 1;
% % % figure
% % % imagesc(GRID)

% % % % Make kernel for finding neighborhood size 
% % % kSize = 9;
% % % kernel = makeSquareKernel2D(GRID,kSize);
% % % ind = find(GRID);
% % % pixNbrhd = repmat(reshape(ind,[1 1 length(ind)]),[kSize kSize])+repmat(kernel,[1 1 length(ind)]);
% % % pixNbrhd = GRID(pixNbrhd);
% % % pixNbrhdCnt = reshape(sum(sum(pixNbrhd,1),2),[size(pixNbrhd,3) 1]);
% % % junk = find(pixNbrhdCnt > 55);
% % % GRID(ind(junk)) = 0.5;
% % % figure
% % % imagesc(GRID)


%% Identify junctions %%
% % % [H,T,R] = hough(BW, 'Theta',-10:10);   %# specific theta range
% % % P  = houghpeaks(H, 5);
% % % lines = houghlines(BW, T, R, P);
% % % % Overlay detected lines over image
% % % figure
% % % imagesc(BW)
% % % hold on
% % % for n = 1:length(lines)
% % %     xy = [lines(n).point1; lines(n).point2];
% % %     plot(xy(:,1), xy(:,2), 'g.-', 'LineWidth',2);
% % % end
% % % hold off
