%% Script to convert VTK surface to OBJ
% % % function vtk2obj(vtk_dir)
%% Read surface from VTK
% % % % Grab current directory path
% % % home = pwd;
% Select vtk file
[geoName,geoPath] = uigetfile('*.vtk','Select the VTK output');
geoName = geoName(1:end-8);
% Load associated geometric into handles variables
[centroids,~,~]=readdat([geoPath geoName 'centroids.dat']);
[norms,~,~]=readdat([geoPath geoName 'normals.dat']);
[neighnum,neighs,~,~]=readneighs([geoPath geoName 'neighs.dat']);   % 2/23/05 MWKay
% neighs is a structure
[pts,numparams,txtparams]=readdat([geoPath geoName 'pts.dat']);
[cells,~,~]=readdat([geoPath geoName 'cells.dat']);
cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one

%% Write surface to OBJ
% % % geoName = 'xyzv_povcyl_max08_Reduced';
% Create file
fid = fopen([geoName(1:end-1) '.obj'],'w');

% Write the header
fprintf(fid,'# OBJ file generated from VTK surface output\n');
fprintf(fid,['# Vertices: ' num2str(size(pts,1)) '\n']);
fprintf(fid,['# Faces: ' num2str(size(centroids,1)) '\n']);
fprintf(fid,'##########################################################\n');

% Write the vertices
formatVert = 'v %f %f %f\n';
fprintf(fid,formatVert,pts');
fprintf(fid,['\n# ' num2str(size(pts,1)) ' vertices\n\n']);

% Write the faces
formatFaces = 'f %d %d %d\n';
fprintf(fid,formatFaces,cells');
fprintf(fid,['\n# ' num2str(size(cells,1)) ' faces']);

% Close file
fclose(fid);


% % % end