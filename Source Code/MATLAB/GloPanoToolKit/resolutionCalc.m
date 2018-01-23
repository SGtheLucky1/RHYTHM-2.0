function [camRes,camResMicro] = resolutionCalc(w)
% Description: This function is used to identify the resolution of each
% camera in the panoramic imaging setup. It does so by accepting the width
% of the base of the calibration cuboid as an input, calculating the
% diagonal of the base with this value, and then identifying the diagonal
% on each image to calculate the real world width of each pixel.
%
% Inputs:
% width - width of the cuboid, all side of the base should be equal to this
% value
%
% Outputs:
% camRes - the camera resoutions in millimeters
% 
% Author: Christopher Gloschat
% Date: April 25, 2016
%
% Change Log:
%
%
%% Code %%
% Save current directory
current_dir = pwd;
% Navigate to the directory where the calibration images are stored
disp('Navigate to the Cube directory.')
cubeDir = uigetdir;
cd(cubeDir)
% Create variables for assembling the image names for load in
camLabels = 'ABCDG';
C = zeros(2,2,5);
im = cell(5,1);
% Instruct the user on how to collect the needed information
disp('Select a point on the left edge of the cube and another on the right side at the same height.')
for n = 1:5
    f = figure;
    eval(sprintf('im{n} = imread(''Cube%s.tiff'');',camLabels(n)))
    % Get the horizontal cross section for each calibration image
    image(im{n})
    C(:,:,n) = ginput(2);
    close(f)
end
% The number of pixels the diagonal is equal to
pix = squeeze(C(2,1,:) - C(1,1,:));
% The length being measured is the diagonal of the base of the cuboid. Both
% side of the base are the length of the width (w). The diagonal (diag) is
% equal to this value squared, multiplied by two, and square rooted.
% Multiplying by 25.4 converts it to millimeters.
diag = sqrt(2*w^2)*25.4;
% The camera resolution
camRes = diag./pix;
% The camera resolutin in micrometers
camResMicro = round(camRes*1000);
% Print values to command window
for n = 1:4
    fprintf('Mapping camera %s has a resolution of approximately %3.0d micrometers.\n',camLabels(n),camResMicro(n))
end
fprintf('The geometry camera has a resolution of approximately %3.0d micrometers.\n',camResMicro(5))
% Navigate back to original directory
cd(current_dir)

end