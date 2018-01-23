% function [] = vtk2obj()
% Description: Converts a *.vtk file format to an *.obj file format
%
% Inputs: None
% 
% Outputs: None
%
% Author: Christopher Gloschat
% Date: July 1, 2016
%
% Modification Log:
%
%% Code %%

% Identify file you want to convert
[fname,fpath,~] = uigetfile('*.vtk','Select the VTK file to convert.');
cd(fpath)

% Open binary file
fid = fopen([fpath fname]);
% Pull string
fstr = fread(fid,'int8=>char');
fstr = fstr';
% Close file
fclose(fid);

% Find end line characters
endLine = strfind(fstr,char(10));
vertStart = endLine(5);
faceStart = strfind(fstr,'POLYGON');

% Grab all end lines associated with vertices
tmp = endLine < vertStart;
vertInd = endLine;
vertInd(tmp) = [];
tmp = vertInd > faceStart;
vertInd(tmp) = [];


% end