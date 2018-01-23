function image = readtiff(filename)
%READPGM Read a raw pgm file as a matrix
%
%        IMAGE = READPGM(FILENAME) reads the binary PGM image data from
%        the file named FILENAME and returns the image as a 2-dimensional
%        array of integers IMAGE. Assumes the file is a raw PGM file
%        containing 8-bit unsigned character data to represent pixel values.
%
% Matthew Dailey, 1997
% Matthew Kay, 2003

% Open the file
A = imread(filename);
A = A(:,:,1)

% Parse and check the header information. 
A = fgets(fid);
if strcmp(A(1:2),'P5') ~= 1
  error('File is not a raw PGM');
end;
A = fgets(fid);
if strcmp(A(1),'#') == 1   % for Gimp PGM files
  A = fgets(fid);
  sizes = sscanf(A,'%d');
  w = sizes(1);
  h = sizes(2);
  A = fgets(fid);
else
  sizes = sscanf(A,'%d');
  w = sizes(1);
  h = sizes(2);
end

%A = fgets(fid);
%max = sscanf(A,'%d');
tlength = w*h;

%if max ~= 255
% error('Cannot handle anything but 8-bit graymaps');
%end; 

% Read the raw data
[v,count] = fread(fid,inf,'uint8');
if count ~= tlength
 error('File size does not agree with specified dimensions.');
end;

% Pack the column vector v into the image matrix
image = reshape(v,w,h)';
image=mat2gray(image,[0 255]);

fclose(fid);
