fid = fopen('sciGeo.obj','w');
fprintf(fid,'# OBJ file generated for panoramic imaging\n');
fprintf(fid,'# Vertices: %d\n',size(pts,1));
fprintf(fid,'# Faces: %d\n',size(cells,1));
fprintf(fid,'##########################################\n');
for n = 1:size(pts,1)
    fprintf(fid,'v %.6f %.6f %0.6f\n',pts(n,1),pts(n,2),pts(n,3));
end
fprintf(fid,'# %d vertices\n\n',size(pts,1));
for n = 1:size(cells,1)
   fprintf(fid,'f %d %d %d\n',cells(n,1),cells(n,2),cells(n,3));
end
fprintf(fid,'\n# %d faces',size(cells,1));
fclose(fid);