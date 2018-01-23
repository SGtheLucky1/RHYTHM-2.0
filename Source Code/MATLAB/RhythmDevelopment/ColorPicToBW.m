myDir = uigetdir;
cd(myDir)
myFiles = dir(fullfile(myDir, '*.tiff'));
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    a = imread(baseFileName);
    a = rgb2gray(a);
    cd('../Heart')
    imwrite(a, baseFileName);
    cd(myDir)
end