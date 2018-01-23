% 08-20-02, MWKay
%
% Revised for increased user-directed analysis with consistent
% saving of analyzed silhouettes as well as user-selected image numbers
% for analysis. 12/08/2004, MWKay
%


disp(sprintf('\n whirl.m \n version 11, 12-08-04 \n'));

clear all
close all

%dtheta=5;
%n_images=360/dtheta;
%hemi=360/(2*dtheta);
%theta=0:dtheta:(360-dtheta);

go=1;
while go
  dtheta=input('Degree step [5]: ');
  if isempty(dtheta), dtheta=5; end;
  n_images_com=sprintf('n_images=input(''Number of snapshots [%d]: '');',floor(360/dtheta));
  eval(n_images_com);
  if isempty(n_images), n_images=floor(360/dtheta); end;
  
  r=0:dtheta:(n_images*dtheta-1);
  if rem(360,dtheta)==0
  
    % Determine front/back positions
    frontback=1;
    rr=zeros(size(r,2),5).*NaN;
    r1=find(r-180<0);
    r2=find(r-180>=0);
    rr(1:length(r1),1)=r(r1)';
    rr(1:length(r2),2)=r(r2)';
    rr(1:length(r1),3)=r1'-1;
    rr(1:length(r2),4)=r2'-1;
    rnot=find(isnan(rr(:,1)) & isnan(rr(:,2)));
    rr(rnot,:)=[]; 
    sprintf('Front/Back positions:')
    rr
    clear r1 r2 rnot
  else
  
    % No front/back images
    frontback=0;
    rr=zeros(size(r,2),5).*NaN;
    rr(:,1)=r'
    sprintf('No Front/Back positions!')
  end
  disp(sprintf('\nTheta_step= %d deg \nNumber of snapshots: %d \nTheta_start= 0 deg \nTheta_final= %d deg\n',dtheta,n_images,dtheta*(n_images-1)))  
  yn=input('Correct? [Y]: ','s');
  if (isempty(yn) | yn=='Y' | yn=='y'), go=0; end;
end

%---------------------------------------------------
% Get the silhouettes

analyzeit=1;
fid=fopen('silhs1.mat');
if fid~=-1
  fclose(fid);
  disp('Found previously analyzed data in silhs1.mat!');
  disp('Loading silhs1.mat ....');
  load silhs1.mat
  notanalyzed=length(find(~analyzed));
  if notanalyzed==0
    disp('All images appear to have been analyzed.');
    goflag=1;
    while goflag
      answer=input('Continue analyzing the images (c) or use the stored analysis for reconstruction (u) [c/u]: ','s');
      if (answer=='c' | answer=='C' | answer=='u' | answer=='U')
        goflag=0;
      end
    end
    if (answer=='u' | answer=='U')
      analyzeit=0;
    end
  else
    disp(sprintf('Found %d images that have not been analyzed.',notanalyzed));
    disp('Continue with silhouette analysis...');
    clear silh1 lims1 centroids1 orientations1 areas1 thresharr analyzed
  end
else
  disp('Did not find previously analyzed data (no silhs1.mat file).');
  disp('Initiating new silhouette analysis...');
end

if analyzeit
  foundfiles=0;
  while ~foundfiles
    studylabel=input('Enter study label (ie: map24): ','s');
    scanlabel=input('Enter scan label (ie: scan4): ','s');
    defbdir=sprintf('~/%s/%s/',studylabel,scanlabel);
    bdir=input(sprintf('images are stored in this directory [%s]: ',defbdir),'s');
    if isempty(bdir); bdir=defbdir; end
    defbfilename=sprintf('%s_',scanlabel);
    bfilename=input(sprintf('basis of image filenames [%s]: ',defbfilename),'s');
    if isempty(bfilename); bfilename=defbfilename; end
    defsfilename='pgm';
    sfilename=input(sprintf('suffix of image filenames [%s]: ',defsfilename),'s');
    if isempty(sfilename); sfilename=defsfilename; end
    defndigits=3;
    ndigits=input(sprintf('number of digits in image filenames (ie: 3 digits in scan4_001.pgm) [%d]: ',defndigits));
    if isempty(ndigits); ndigits=defndigits; end
    sdigit=-1;
    while (sdigit~=0 & sdigit~=1)
      sdigit=input('start digit (must be 0 or 1): ');
      if isempty(sdigit); sdigit=-1; end
    end
    fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,sdigit,sfilename);',ndigits);
    eval(fnamecom);
    fid=fopen(fname);
    if fid==-1
      disp(sprintf('Cannot open %s!',fname));
    else
      fclose(fid);
      foundfiles=1;
    end
  end
  [silh1,lims1,centroids1,orientations1,areas1,thresharr,analyzed]=get_silhs(bdir,bfilename,sfilename,ndigits,sdigit,n_images);
end

notanalyzed=length(find(~analyzed));
if notanalyzed>0
  disp(sprintf('%d images have not been analyzed.',notanalyzed));
  disp('Exiting reconstruction ...');
  return
end

if frontback
  go=1;
  while go
    yn=input('Use largest silhouettes to collapse redundant \ninformation in front/back snapshots? [N]: ','s');
    if (isempty(yn) | yn=='N' | yn=='n');
      go=0;
      dofrontback=0; 
    elseif (yn=='Y' | yn=='y')
      go=0;
      dofrontback=1; 
    end
  end
  if dofrontback
    for i=1:size(rr,1)
      ars=[areas1(rr(i,3)) areas1(rr(i,4))];
      maxarsi=find(ars==max(ars));
      if maxarsi==1, rr(i,5)=rr(i,3); end;
      if maxarsi==2, rr(i,5)=rr(i,4); end;
    end
    sprintf('First 2 columns are snapshot angle, \nNext 2 columns are corresponding snapshot numbers, \nLast column is snapshot number of larger silhouette:')
    rr
    sprintf('Snapshots taken at these angles will be used. \nFirst column is snapshot number, \nsecond column is the angle. ')
    [rsort,isort]=sort(r(rr(:,5)));
    inumsort=rr(isort,5);
    irsort=[inumsort rsort'];
    [inumsort rsort']
  else
    sprintf('Snapshots taken at these angles will be used. \nFirst column is snapshot number, \nsecond column is the angle. ')
    rsort=r;
    inumsort=1:length(r);
    irsort=[inumsort' rsort'];
    [inumsort' rsort']
  end  
end
silh=silh1;
lims=lims1;
clear silh1 lims1 centroids1 orientations1 areas1 thresharr ars maxarsi

%---------------------------------------------------

% load camera calibration parameters
disp(sprintf('\n Load camera parameters:\n'));

defRcfname='../cal/Rc.dat';
defParfname='../cal/Par.dat';

gofname=1;
while gofname
  Rcfname=input(sprintf('Enter Rc.dat path and filename [%s]: ',defRcfname),'s');
  if isempty(Rcfname); Rcfname=defRcfname; end
  fid=fopen(Rcfname);
  if fid==-1
    disp(sprintf('Unable to open %s!',Rcfname));
  else
    fclose(fid);
    disp(sprintf('Found %s',Rcfname));
    gofname=0;
  end
end

gofname=1;
while gofname
  Parfname=input(sprintf('Enter Par.dat path and filename [%s]: ',defParfname),'s');
  if isempty(Parfname); Parfname=defParfname; end
  fid=fopen(Parfname);
  if fid==-1
    disp(sprintf('Unable to open %s!',Parfname));
  else
    fclose(fid);
    disp(sprintf('Found %s',Parfname));
    gofname=0;
  end
end

disp(sprintf('\n Loading camera parameters....\n'))

[Rc,Rc_ntheta,Rc_dtheta,message]=readRccal(Rcfname);
if dtheta~=Rc_dtheta 
  disp(sprintf('\n dtheta does not match Rc_dtheta!! \n Aborting ...\n'));
  return
end
[Par,Par_ntheta,Par_dtheta,message]=readParcal(Parfname);
if dtheta~=Par_dtheta 
  disp(sprintf('\n dtheta does not match Par_dtheta!! \n Aborting ...\n'));
  return
end

disp(sprintf('\n Camera parameters loaded.\n'));
disp(sprintf('\n Was the object scanned in the same direction as the camera'));
disp(sprintf(' was calibrated? ie: camera calibrated cw and object scanned ccw? or'));
disp(sprintf('                     camera calibrated ccw and object scanned cw?'));
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
for i=1:size(Rc,3)
  xyz(i,:)=(inv(Rc(1:4,1:4,i))*[0 0 0 1]')';
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
  pause(0.1);
end
clear xyz;

%----------------------------------------------------

% load octree basis
octree_basis

% Establish size and origin of initial carving cube 
% Keep it cubic!!!  This makes scaling volume and surface area much easier!     
xdim=6.0*25.4;  % mm
ydim=xdim;    % mm
zdim=xdim;    % mm

X0=-0.5*xdim;
Xn=0.5*xdim;
Y0=-0.5*ydim;
Yn=0.5*ydim;
Z0=0.45*zdim;
Zn=-0.55*zdim;

levs=[[1:10]' (Xn-X0)./(2.^[1:10]')];        % Level number, delt
[[1:10]' (Xn-X0)./(2.^[1:10]')]
startlevel=input('Start at level number: ');
maxlevel=input('Max level number: ');

fnametag=input('Filename tag (ie: _povcyl_1mm_corrected): ','s');
savedir=input('Save datafiles in this directory, leave off the "/" (ie: povcyl_05mm): ','s');

% verticies of voxel h-1 of octree t is (remember, voxels are from 2:9):
% vert(octree(t,find(octree(t,:,h)),1),:);

level=startlevel;
delt(1,1)=(Xn-X0)/(2^level);
delt(1,2)=(Yn-Y0)/(2^level);
delt(1,3)=(Zn-Z0)/(2^level);

%--------------------------------------------
% Record analysis parameters in textfile
parafname=sprintf('%s/whirl%s.txt',savedir,fnametag);
parafid=fopen(parafname,'w');
fprintf(parafid,'dtheta, degree step: %4.3f \n',dtheta);
fprintf(parafid,'n_images, total number of snapshots: %d \n',n_images);
fprintf(parafid,'dofrontback, 0: use all silhouettes or 1: use only largest silhouettes: %d \n',dofrontback);
fprintf(parafid,'xdim, size of clay cube in x dir (mm): %4.3f \n',xdim);
fprintf(parafid,'(X0,Xn) mm: (%4.3f,%4.3f) \n',X0,Xn);
fprintf(parafid,'ydim, size of clay cube in y dir (mm): %4.3f \n',ydim);
fprintf(parafid,'(Y0,Yn) mm: (%4.3f,%4.3f) \n',Y0,Yn);
fprintf(parafid,'zdim, size of clay cube in z dir (mm): %4.3f \n',zdim);
fprintf(parafid,'(Z0,Zn) mm: (%4.3f,%4.3f) \n',Z0,Zn);
fprintf(parafid,'Levels: \n');
for i=1:size(levs,1)
  for j=1:size(levs,2)
    fprintf(parafid,'%4.6f   ',levs(i,j));
  end
  fprintf(parafid,'\n');
end
fprintf(parafid,'startlevel, start at this spacing level: %d \n',startlevel);
fprintf(parafid,'maxlevel, stop after completion of this spacing level: %d \n',maxlevel);
fprintf(parafid,'fnametag: %s \n',fnametag);
fprintf(parafid,'savedir: %s \n',savedir);
fprintf(parafid,'Camera calibration, Position file: %s \n',Rcfname);
fprintf(parafid,'Camera calibration, Parameters file: %s \n',Parfname);
fprintf(parafid,'Angle array (r):\n');
for i=1:length(r)
  fprintf(parafid,'%4.6f \n',r(i));
end
fprintf(parafid,'Hemisphere array (rr):\n');
for i=1:size(rr,1)
  for j=1:size(rr,2)
    fprintf(parafid,'%4.6f   ',rr(i,j));
  end
  fprintf(parafid,'\n');
end
fprintf(parafid,'Sorted snapshot array (irsort):\n');
for i=1:size(irsort,1)
  for j=1:size(irsort,2)
    fprintf(parafid,'%4.6f   ',irsort(i,j));
  end
  fprintf(parafid,'\n');
end

%--------------------------------------------

goodthresh=0.95;   % threshold level of v indicating which verticies are in or out!
close(tt);

while level<=maxlevel 
%--------------------------------------------
  % Create the initial verticies matrix (verts) and octree matrix (octree)
  % vert(:,1)=x
  % vert(:,2)=y
  % vert(:,3)=z
  % vert(:,4)=v (0 for outside, 1 for inside, NaN if not checked
  % vert(:,5)=level, resolution level when this vertex was added
  % octree is 3 dimensional: (Num of octrees) x (27 indicies denoting rows in vert) 
  % x (9, where octree(:,:,1) contains vert indicies and octree(:,:,2:9) denotes voxel verticies
  % octree(octree_num,:,1)=vert indicies
  % octree(octree_num,:,2:9)=voxel verticies
  % At each level the octree matrix is rebuilt and the vert matrix is expanded.
  
  if (level==1 & level==startlevel) 
    disp(sprintf('Creating octree for level 1...\n'))
    vert=octb*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(27,3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert=vert+sshift;
    clear sshift
    vert(:,4)=NaN;
    vert(:,5)=1;
    octree=zeros(1,27);
    octree=1:27;
  elseif (level==2 & level==startlevel)
    disp(sprintf('Loading level 2 octree from datafile...\n'));
    load /home/mwk/matlab/Panoramic/level2
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==3 & level==startlevel)
    disp(sprintf('Loading level 3 octree from datafile...\n'));
    load /home/mwk/matlab/Panoramic/level3
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==4 & level==startlevel)
    disp(sprintf('Loading level 4 octree from datafile...\n'));
    load /home/mwk/matlab/Panoramic/level4
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==5 & level==startlevel)
    disp(sprintf('Loading level 5 octree from datafile...\n'));
    load /home/mwk/matlab/Panoramic/level5
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==6 & level==startlevel)
    disp(sprintf('Loading level 6 octree from datafile...\n'));
    load /home/mwk/matlab/Panoramic/level6
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  %--------------------------------------------
  else
    % which voxels were inside AND outside the volume?
    % find them, then subdivide and repeat for each octree
    % to which voxels were those that were inside belong
    % 'good' tells me which verticies were inside the volume
    % and which verticies were outside the volume of the verticies
    % which have not been checked.
    
    sprintf('Creating octree for level %d.... ',level)
    
    good=find(vert(:,4)>=goodthresh);  % Voxels to be subdivided have at least
                                 % one (but not all) vertex greater than 0.99
				 % Use 0.95 to provide good surface definition.

    % Grab the vertex index from good.
    % Find which octrees contain that index
    clear ovi;
    for j=1:size(octree,1)
      for i=1:length(good)
	[ovi1,ovi2]=find(octree(j,:)==good(i));
	if ~isempty(ovi1)
          if ~exist('ovi')
	    ovi=[j ovi2];
	  else
	    ovi=[ovi; j ovi2];
	  end
	end
      end
    end

    % ovi- octree and vertex indicies
    % ovi(:,1) contains octree numbers for verticies inside the volume
    % ovi(:,2) contains octree vertex numbers for the verticies

    % with which octrees are we dealing?
    octs=unique(ovi(:,1));

    % now divide and conquer
    for i=1:length(octs)   % loop through each octree
      % get vertex numbers for this octree that are inside the volume
      vs=ovi(find(ovi(:,1)==octs(i)),2);    
      % which voxels have these vertex numbers?
      % first built a 'mini' octbi in os
      os=octbi(vs,:);
      % next find the voxel numbers
      voxelsi=find(max(os,[],1)>0);  % these are actual voxel numbers for the respective octree!
      % is a voxel totally inside the volume?
      vflag=0;
      for j=1:length(voxelsi)
	%if length(find(ismember([1 2 3 4 5 6 7 8],os(:,voxelsi(j)))))==8     % CAN I SPEED THIS UP????
	if length(unique(os(:,voxelsi(j))))==9
          voxelsi(j)=NaN; 
	  vflag=1;   % We have dropped one or more!
	  %sprintf('Dropped one!')
	end
      end
      if vflag voxelsi=voxelsi(find(~isnan(voxelsi))); end % drop voxels totally inside the volume
      % we now have the voxel number of each voxel of the octree 
      % (stored in voxelsi) that is both inside and outside the volume
      % now subdivide each voxel, rebuild the 'octree' matrix 
      % and append to the 'vert' matrix

      for j=1:length(voxelsi)              % loop through each voxel of interest in the current octree, all voxels are the same size
	vni=find(octbi(:,voxelsi(j))>0);   % find vertex number indicies of the voxel of interest
	vvi=octree(octs(i),vni);           % actual indicies into vert for each vertex  ***REVISED 1/29/02***
	vxyz=vert(vvi,1:3);                % vertex locations are in order!
	
	% now create an octree inside the voxel
	% first fill in the verticies that already exist!
	if i==1 & j==1
          okount=1;
          new_octree=zeros(1,27).*NaN;
	  new_octree([1 4 9 11 17 19 24 27])=vvi;
	  delt(level,1)=(vxyz(2,1)-vxyz(1,1))/2;
	  delt(level,2)=(vxyz(5,2)-vxyz(1,2))/2;
	  delt(level,3)=(vxyz(3,3)-vxyz(1,3))/2;
	else
          okount=okount+1;
          new_octree(okount,[1 4 9 11 17 19 24 27])=vvi;  % zeros elsewhere
	end
        
	xn=vxyz(2,1);
	x0=vxyz(1,1);
	yn=vxyz(5,2);
	y0=vxyz(1,2);
	zn=vxyz(3,3);
	z0=vxyz(1,3);
	 
	new_vert=zeros(19,5).*NaN;
	new_vert(:,5)=level;
	new_vert(:,1:3)=octb([2 3 5 6 7 8 10 12 13 14 15 16 18 20 21 22 23 25 26],:)*[(xn-x0) 0 0; 0 (yn-y0) 0; 0 0 (zn-z0)];
        sshift=zeros(19,3);
        sshift(:,1)=x0;
        sshift(:,2)=y0;
        sshift(:,3)=z0;
        new_vert(:,1:3)=new_vert(:,1:3)+sshift;

	new_vert_ii=[1:19]+size(vert,1);

	% do any of these verticies already exist in vert? 
	% if so then delete from new_vert and update new_vert_ii
	[C,rnew_vert,rvert]=intersect(new_vert(:,1:3),vert(:,1:3),'rows');  % rows of new_vert that are already in vert
	clear C;

	if ~isempty(rnew_vert)   % if redundant verticies exist
          new_vert_ii(rnew_vert)=rvert;  % This makes sense
	  new_new_vert=new_vert;
	  new_new_vert(rnew_vert,:)=[];  % drop redundant rows in new_vert
          [C,r1,r2]=intersect(new_new_vert(:,1:3),new_vert(:,1:3),'rows');
	  clear C;
	  new_vert_ii(r2)=r1+size(vert,1);
	  new_vert=new_new_vert;
	end
	
	new_octree(okount,[2 3 5 6 7 8 10 12 13 14 15 16 18 20 21 22 23 25 26])=new_vert_ii;
	vert=[vert; new_vert];         
	
	if length(unique(new_octree(okount,:)))~=27
	  sprintf('Junk!')
	  return
	end	
      end  % End voxel loop (j)
    end % End octree loop (i), okount is combination of i and j
    octree=new_octree;
  end
  %--------------------------------------------

  % vert(:,4) of NaN: vertex not checked yet
  % vert(:,4) of 0: vertex is definately not in the volume
  
  sites2check=find(isnan(vert(:,4))); % indicies into vert
  check_num=length(sites2check);
  check=ones(check_num,2);
  check(:,1)=sites2check;
  
  % inumsort is the snapshot number to check
  % rsort is the angle corresponding to each snapshot number in inumsort
  
  for i=1:length(inumsort)     
  %for i=1:1
    good=find((check(:,2))>0);   % must shave anything greater than zero that is not inside b/c projection
    if isempty(good)             % places verticies that are outside inside, depending upon view!
      sprintf('No points inside volume. Aborting Snapshots. Restart and increase initial level.')
      return
    else
      disp(sprintf('Level %d, Snapshot #%d, %d deg.\nThis is snapshot %d of %d, checking %d verticies ... ',level,inumsort(i),rsort(i),i,length(inumsort),length(good)));
      angi=find(Rc(1,5,:)==rsort(i)); 
      disp(sprintf('  Applying camera model ...'));
      [Xi,Yi]=pred2([vert(check(good,1),1) vert(check(good,1),2) vert(check(good,1),3) ones(length(good),1)],Rc(:,1:4,angi),Par(1:8,angi));
      % [Xi,Yi]=snapcube(vert(check(good,1),1),vert(check(good,1),2),vert(check(good,1),3),pos,par,rsort(i));
      disp(sprintf('  Determining which points are inside the volume ...'));
      % Which points of Xi and Yi are inside the silhouette?
      vicinity=ones(length(good),1);
      vicinity(find(Xi<lims(inumsort(i),1) | Xi>lims(inumsort(i),2) | Yi<lims(inumsort(i),3) | Yi>lims(inumsort(i),4)))=0;
      check(good(find(~vicinity)),2)=0;
      vici=find(vicinity);
      Zi=interp2(silh(:,:,inumsort(i)),Xi(vici),Yi(vici));
      check(good(vici),2)=Zi;
      disp(sprintf('  Done with snapshot #%d',inumsort(i)));
      
      if 0
        if i==1, hh=figure; end;
        figure(hh);
        imagesc(silh(:,:,inumsort(i))); hold on;
        plot(Xi,Yi,'y.','MarkerSize',2);
        plot(Xi(vici(find(Zi)>0)),Yi(vici(find(Zi)>0)),'wo','MarkerSize',2);
        hold off;	
        pause;
        if i==length(inumsort), close(hh); end;
      end
      
    end
  end
  vert(check(:,1),4)=check(:,2);
  
  % keep up with the points that have been tested
  % vert(:,4) of ~NaN means that the vertex has been tested
  % vert(:,4) of NaN means that the vertex has not been tested
  
  level=level+1
end % End level loop

%----------------

% No need to normalize but code stays for posterity

% Normalize and save data
vertn=zeros(size(vert,1),size(vert,2)-1);
vertn(:,1)=vert(:,1)-min(vert(:,1));
vertn(:,2)=vert(:,2)-min(vert(:,2));
vertn(:,3)=vert(:,3)-min(vert(:,3));
vertn(:,4)=vert(:,4);

scale_x=max(vertn(:,1));
scale_y=max(vertn(:,2));
scale_z=max(vertn(:,3));

vertn(:,1)=vertn(:,1)./scale_x;
vertn(:,2)=vertn(:,2)./scale_y;
vertn(:,3)=vertn(:,3)./scale_z;

savecom=sprintf('save %s/octvert%s octree vert delt level xdim ydim zdim scale_x scale_y scale_z rr inumsort rsort',savedir,fnametag);
eval(savecom);

scalesfname=sprintf('%s/scales%s.dat',savedir,fnametag);
fid=fopen(scalesfname,'w');
fprintf(fid,'%4.3f\n',scale_x);
fprintf(fid,'%4.3f\n',scale_y);
fprintf(fid,'%4.3f\n',scale_z);
fclose(fid);

xyzvfname=sprintf('%s/xyzv%s.dat',savedir,fnametag);
fid=fopen(xyzvfname,'w');
fwrite(fid,vert(:,1:4)','float');   % Edited 3/12/02 from fwrite(fid,vertn','float'); 
fclose(fid);

fprintf(parafid,'Saved %d points in %s \n',size(vert,1),xyzvfname);
fclose(parafid);

return
%--------------------------------------------

sprintf('Rendering ....')

close all
clear all

load octvert

xsign=sign(delt(end,1));
ysign=sign(delt(end,2));
zsign=sign(delt(end,3));

inside=find(vert(:,4)>=0.80);  % Get the size of volume for meshgrid
minx=min(xsign.*vert(inside,1)).*xsign;
maxx=max(xsign.*vert(inside,1)).*xsign;
miny=min(ysign.*vert(inside,2)).*ysign;
maxy=max(ysign.*vert(inside,2)).*ysign;
minz=min(zsign.*vert(inside,3)).*zsign;
maxz=max(zsign.*vert(inside,3)).*zsign;

[x,y,z]=meshgrid(minx:delt(end,1)/2:maxx,miny:delt(end,2)/2:maxy,minz:delt(end,3)/2:maxz);
X=[x(:), y(:), z(:)];
v=griddatan(vert(:,1:3),vert(:,4),X);
v=reshape(v,size(x));
vs=smooth3(v);

gg=find(vert(:,4)>=0.90);
plot3(vert(gg,1),vert(gg,2),vert(gg,3),'b.')
hold on

p=patch(isosurface(x,y,z,vs,0.90)); 
isonormals(x,y,z,vs,p); 
%set(p,'FaceColor','red','EdgeColor','none'); 
set(p,'FaceColor','red');

pcaps=patch(isocaps(x,y,z,vs,0.90));
set(pcaps,'FaceColor','red');

view(3); 
camlight; 
lighting phong 





% OLD STUFF DOWN HERE -----------------
% vert(:,4)>=0.5 and vert(:,4)<0.9: vertex might be in the volume but check again
% vert(:,4)>=0.9: vertex is definately in the volume

% testa=zeros(size(vert,1),1);
% testb=testa;
% testa=isnan(vert(:,4));
% testb(find(vert(:,4)>=0.5 & vert(:,4)<0.9))=1;
% sites2check=find(testa | testb);  % indicies into vert
