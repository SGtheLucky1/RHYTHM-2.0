% Apply surface texture to reconstructed objects.
% Uses snapnspin snapshots
%
% 08/28/02
% MWKay

close all
clear all

basefname='scan2/xyzv_1.2mm';
Parfname='cal/Par.dat';
Rcfname='cal/Rc.dat';
bfilename=input('Basis of image filenames(ie: images/cube is basis for images/cube001.bmp ): ','s');

% Establish filenames -----------------------
ptsfname=sprintf('%s_pts.dat',basefname);
cellsfname=sprintf('%s_cells.dat',basefname);
normsfname=sprintf('%s_norms.dat',basefname);
centroidsfname=sprintf('%s_centroids.dat',basefname);
neighborsfname=sprintf('%s_neighs.dat',basefname);
txtyourfname=sprintf('%s_texture.dat',basefname);

% Read the data -----------------------
%disp(sprintf('Reading %s',ptsfname));
%[pts,numparams,txtparams]=readdat(ptsfname);
%disp(sprintf('Reading %s',cellsfname));
%[cells,numparams,txtparams]=readdat(cellsfname);
disp(sprintf('Reading %s',normsfname));
[norms,numparams,txtparams]=readdat(normsfname);
numnorms=size(norms,1);
disp(sprintf('Loaded %d normals.',numnorms));

% Read calibration data -----------------------
disp(sprintf('Reading %s',Rcfname));
[Rc,Rc_ntheta,Rc_dtheta,message]=readRccal(Rcfname);
disp(sprintf('Reading %s',Parfname));
[Par,Par_ntheta,Par_dtheta,message]=readParcal(Parfname);
disp(sprintf('\nCamera parameters loaded.\n'));
disp(sprintf('\nWas the object scanned in the same direction as the camera'));
disp(sprintf('was calibrated? ie: camera calibrated cw and object scanned ccw? or'));
disp(sprintf('                    camera calibrated ccw and object scanned cw?'));
disp(sprintf('If so, then camera calibration parameters need to be re-ordered.'));
reorder='g';
while (isempty(reorder) | (reorder~='y' & reorder~='Y' & reorder~='n' & reorder~='N'))
  reorder=input('Should camera calibration parameters be re-ordered? (y/n): ','s');
end
if (reorder=='y' | reorder=='Y')
  disp(sprintf('\n Reordering camera parameters .... \n'));
  Par(1:8,2:end)=fliplr(Par(1:8,2:end));
  Rc_new=Rc;
  for i=2:size(Rc,3)
    Rc_new(1:4,1:4,i)=Rc(1:4,1:4,size(Rc,3)-i+2);
  end
  Rc=Rc_new;
  clear Rc_new
else
  disp(sprintf('\n Ok, will not reorder camera parameters .... \n'));
end
% Plot orbit to check spin and angle consistency
disp('Check camera rotation .... ');
xyz=ones(size(Rc,3),4);
look=ones(size(Rc,3),4);
lookn=ones(size(Rc,3),3);
for i=1:size(Rc,3)
  xyz(i,:)=(inv(Rc(1:4,1:4,i))*[0 0 0 1]')';             % camera origin in world coords
  look(i,:)=(inv(Rc(1:4,1:4,i))*[0 0 1 1]')';            % camera view direction
  look(i,1:3)=look(i,1:3)-xyz(i,1:3);                    % camera view direction
  lookn(i,:)=look(i,1:3)./norm(look(i,1:3));           % normalized view direction
end
tt=figure; hold on; axis('square');
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) min(xyz(:,3)) max(xyz(:,3)) ]); 
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
for i=1:size(Rc,3)
  disp(sprintf('Angle %d',Rc(1,5,i)));
  figure(tt);
  if i==1, text(xyz(i,1),xyz(i,2),xyz(i,3),sprintf('%d',Rc(1,5,i)));
  else, text(xyz(i,1),xyz(i,2),xyz(i,3),sprintf('%d',Rc(1,5,i))); end;
  quiver3(xyz(i,1),xyz(i,2),xyz(i,3),look(i,1),look(i,2),look(i,3));
  pause(0.1);
end
clear xyz;
clear look;
close(tt);

disp('Computing angles between normals and view directions .... ');
langle=zeros(numnorms,72);
disp(sprintf('%d normals and %d views',size(langle,1),size(langle,2)));
two_comp=0;   % 1 for 2 components, 0 for three
if ~two_comp
  disp('Computing angles using x,y, and z components.');
  gg=norms*lookn';
  gg=acos(gg);
  gg=real(gg);
  langle=gg.*180/pi;  % x,y, and z components
elseif two_comp 
  disp('Computing angles using x and y components.');
  lookn=lookn(:,1:2);
  norms=norms(:,1:2);
  gg=norms*lookn';
  gg=acos(gg);
  gg=real(gg);
  langle=gg.*180/pi;  % only x and y components
end
clear gg;
clear lookn;
clear norms;


return

% Assign views to cells
disp('Assigning views to cells ...');
[Y,I]=sort(langle,2);
clear langle;
clear Y;
I=I(:,end);             % I is view number for each cell (the max langle)
disp(sprintf('Reading %s',neighborsfname));
[neighs,numparams,txtparams]=readdat(neighborsfname);
neighs=neighs+1;        % VTK starts with zero, matlab starts with one
disp(sprintf('Loaded %d neighbors.',size(neighs,1)));
disp('Assigning surface regions ...');
I2=I;
kount=0;
for i=1:numnorms
  if (I2(neighs(i,1))-I2(neighs(i,2)))==0
    I2(i)=I2(neighs(i,1));
    kount=kount+1;
  elseif (I2(neighs(i,1))-I2(neighs(i,3)))==0;
    I2(i)=I2(neighs(i,1));
    kount=kount+1;
  elseif (I2(neighs(i,2))-I2(neighs(i,3)))==0;
    I2(i)=I2(neighs(i,2));
    kount=kount+1;
  end
end

% load image filenames ----------------
n_images=size(Rc,3);
dtheta=360/n_images;
imnames=cell(n_images,1);
imnamesfid=fopen(sprintf('%s_fnames.txt',bfilename),'r');
disp(sprintf('Loading image filenames from %s_fnames.txt',bfilename));
for i=1:n_images
  imnames(i)=cellstr(fgets(imnamesfid));
end
fclose(imnamesfid);

disp(sprintf('Reading %s',centroidsfname));
[centroids,numparams,txtparams]=readdat(centroidsfname);
disp(sprintf('Loaded %d centroids.',size(centroids,1)));

% get texture from images
txtyour=zeros(numnorms,1);
for i=1:n_images
  disp(sprintf('Computing texture from image %s ...',char(imnames(i))));
  disp(sprintf('  Using a calibration angle of %d deg.',Rc(1,5,i)));
  a=rgb2gray(imread(char(imnames(i))));
  gg=find(I2==i);
  disp(sprintf('  Found %d centroids: ',length(gg)));
  disp('    Taking snapshot of model (ie: perspective projection) ...');
  [Xi,Yi]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(length(gg),1)],Rc(:,1:4,i),Par(1:8,i));
  disp('    Interpolating original image to obtain texture values ... ');
  txtyour(gg)=interp2(double(a),Xi,Yi); 
  
  figure(1);
  imagesc(a); colormap('gray'); hold on;
  plot(Xi,Yi,'w.');
  pause
  close(1);
  
end
clear a gg;
clear Xi Yi;
clear imnames txtyourcount;


% Save the texture!
mediantxtyour=median(avgtxtyour(:,1));
stdtxtyour=std(avgtxtyour(:,1));
centroids(:,4)=avgtxtyour(:,1);
centroids(centroids(:,4)>mediantxtyour+stdtxtyour,4)=mediantxtyour+stdtxtyour;
centroids(centroids(:,4)<mediantxtyour-stdtxtyour,4)=mediantxtyour-stdtxtyour;
centroids(:,4)=im_norm(centroids(:,4),0,1);
fid=fopen(txtyourfname,'w');
fwrite(fid,centroids(:,1:4)','float');
fclose(fid);
disp(sprintf('Saved texture in %s',txtyourfname));





disp(sprintf('Reading %s',neighborsfname));
[neighs,numparams,txtparams]=readdat(neighborsfname);
disp(sprintf('Loaded %d sets of neighbors.',size(neighs,1)));


