function [points,polys]=readdat(fname)
%
% [points, polys]=readvtk('scan2_d.vtk');
%
%
% MWKay, 4/2005
% revised to also read polygons, MWKay 03/2006
%

skiplines=4;

fid=fopen(fname,'r','n');
for i=1:skiplines
  l=lower(fgets(fid));
end
[charbuff,delim]=readvtkchar(fid,32,10);
[pointstr,delim]=readvtkchar(fid,32,10);
numpoints=str2num(pointstr);
[charbuff,delim]=readvtkchar(fid,32,10);

points=zeros(numpoints,3);

for i=1:numpoints
   S=fgets(fid);
   points(i,:)=sscanf(S,'%f')';
end


l=lower(fgets(fid)); % skip 1 line
[charbuff,delim]=readvtkchar(fid,32,10);
[polystr,delim]=readvtkchar(fid,32,10);
numpolys=str2num(polystr);
[charbuff,delim]=readvtkchar(fid,32,10);

polys=zeros(numpolys,4);
for i=1:numpolys
  S=fgets(fid);
  polys(i,:)=sscanf(S,'%d')';
end

fclose(fid);
