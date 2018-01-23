% 05-15-02
% MWKay

disp(sprintf('\n whirl.m \n version 7, 05-14-02 \n'));

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
    rr(1:length(r1),3)=r1';
    rr(1:length(r2),4)=r2';
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
  sprintf('Theta_step= %d deg \nNumber of snapshots: %d \nTheta_start= 0 deg \nTheta_final= %d deg',dtheta,n_images,dtheta*(n_images-1))  
  yn=input('Correct? [Y]: ','s');
  if (isempty(yn) | yn=='Y' | yn=='y'), go=0; end;
end

%---------------------------------------------------
% Get the silhouettes

fid=fopen('silhs1.mat');
if fid~=-1
  sprintf('Loading silhs1.mat ....')
  fclose(fid);
  load silhs1
else 
  sprintf('Could not find silhs1.mat: Creating silhouette matrix....')
  bfilename=input('basis of image filenames(ie: cube is basis for cube001.bmp ): ','s');
  sfilename=input('suffix of image filenames(ie: bmp is suffix for cube001.bmp ): ','s');
  ndigits=input('number of digits in image filenames(ie: 3 digits in cube001.bmp ): ');
  sdigit= input('start digit (ie: 0 or 1): ');
  [silh1,lims1,centroids1,orientations1,areas1,thresharr]=get_silhs(bfilename,sfilename,ndigits,sdigit,n_images);
end
if frontback
  go=1;
  while go
    yn=input('Use largest silhouettes to collapse redundant \ninformation in front/back snapshots? [Y]: ','s');
    if (isempty(yn) | yn=='Y' | yn=='y')
      go=0;
      dofrontback=1; 
    elseif (yn=='N' | yn=='n');
      go=0;
      dofrontback=0; 
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
Rcfname=input('Rc.dat file (ie: ../pre_cal/Rc.dat): ','s');
Parfname=input('Par.dat file (ie: ../pre_cal/Par.dat): ','s');
disp(sprintf('\n Loading camera parameters....\n'))

% Load Rc
fid=fopen(Rcfname,'r');
Rc_dtheta=str2num(fgets(fid));
Rc_ntheta=str2num(fgets(fid));
Rc_message=fgets(fid);
Rc_dat=fread(fid,5*4*Rc_ntheta,'double');
fclose(fid);
Rc=reshape(Rc_dat,4,5,Rc_ntheta);
if dtheta~=Rc_dtheta 
  disp(sprintf('\n dtheta does not match Rc_dtheta!! \n Aborting ...\n'));
  return
end

% Load Par
fid=fopen(Parfname,'r');
Par_dtheta=str2num(fgets(fid));
Par_ntheta=str2num(fgets(fid));
Par_message=fgets(fid);
Par_dat=fread(fid,10*Par_ntheta,'double');
fclose(fid);
Par=reshape(Par_dat,10,Rc_ntheta);
if dtheta~=Par_dtheta 
  disp(sprintf('\n dtheta does not match Par_dtheta!! \n Aborting ...\n'));
  return
end

disp(sprintf('\n Camera parameters loaded.\n'));
disp(sprintf('\n Assuming camera parameters are oriented COUNTER-CLOCKWISE from zero!\n'));
cw=0;
while (cw~=1 & cw~=2)
  cw=input('Are images ordered cw or ccw from zero? (1:cw, 2:ccw): ');
end
if cw==1
  disp(sprintf('\n Reordering camera parameters for CLOCKWISE rotation.... \n'));
  Par=fliplr(Par);
  Rc_new=Rc;
  for i=1:size(Rc,3)
    Rc_new(:,:,i)=Rc(:,:,size(Rc,3)-i+1);
  end
  Rc=Rc_new;
  clear Rc_new
end

%----------------------------------------------------

% load octree basis
octree_basis

% Establish size and origin of initial carving cube 
% Keep it cubic!!!  This makes scaling volume and surface area much easier!     
xdim=8.0*25.4;  % mm
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
fid=fopen(parafname,'w');
fprintf(fid,'dtheta, degree step: %4.3f \n',dtheta);
fprintf(fid,'n_images, total number of snapshots: %d \n',n_images);
fprintf(fid,'dofrontback, 0: use all silhouettes or 1: use only largest silhouettes: %d \n',dofrontback);
fprintf(fid,'xdim, size of clay cube in x dir (mm): %4.3f \n',xdim);
fprintf(fid,'(X0,Xn) mm: (%4.3f,%4.3f) \n',X0,Xn);
fprintf(fid,'ydim, size of clay cube in y dir (mm): %4.3f \n',ydim);
fprintf(fid,'(Y0,Yn) mm: (%4.3f,%4.3f) \n',Y0,Yn);
fprintf(fid,'zdim, size of clay cube in z dir (mm): %4.3f \n',zdim);
fprintf(fid,'(Z0,Zn) mm: (%4.3f,%4.3f) \n',Z0,Zn);
fprintf(fid,'Levels: \n');
for i=1:size(levs,1)
  for j=1:size(levs,2)
    fprintf(fid,'%4.6f   ',levs(i,j));
  end
  fprintf(fid,'\n');
end
fprintf(fid,'startlevel, start at this spacing level: %d \n',startlevel);
fprintf(fid,'maxlevel, stop after completion of this spacing level: %d \n',maxlevel);
fprintf(fid,'fnametag: %s \n',fnametag);
fprintf(fid,'savedir: %s \n',savedir);
fprintf(fid,'Camera calibration, Position file: %s \n',Rcfname);
fprintf(fid,'Camera calibration, Parameters file: %s \n',Parfname);
fprintf(fid,'Angle array (r):\n');
for i=1:length(r)
  fprintf(fid,'%4.6f \n',r(i));
end
fprintf(fid,'Hemisphere array (rr):\n');
for i=1:size(rr,1)
  for j=1:size(rr,2)
    fprintf(fid,'%4.6f   ',rr(i,j));
  end
  fprintf(fid,'\n');
end
fprintf(fid,'Sorted snapshot array (irsort):\n');
for i=1:size(irsort,1)
  for j=1:size(irsort,2)
    fprintf(fid,'%4.6f   ',irsort(i,j));
  end
  fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------

goodthresh=0.95;   % threshold level of v indicating which verticies are in or out!

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
    sprintf('Creating octree for level 1...')
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
    sprintf('Loading level 2 octree from datafile...')
    load /home/mwk/matlab/Panoramic/level2
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==3 & level==startlevel)
    sprintf('Loading level 3 octree from datafile...')
    load /home/mwk/matlab/Panoramic/level3
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==4 & level==startlevel)
    sprintf('Loading level 4 octree from datafile...')
    load /home/mwk/matlab/Panoramic/level4
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==5 & level==startlevel)
    sprintf('Loading level 5 octree from datafile...')
    load /home/mwk/matlab/Panoramic/level5
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==6 & level==startlevel)
    sprintf('Loading level 6 octree from datafile...')
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
    good=find((check(:,2))>0);   % must shave anything greater than zero that is not inside b/c projection
    if isempty(good)             % places verticies that are outside inside, depending upon view!
      sprintf('No points inside volume. Aborting Snapshots. Restart and increase initial level.')
      return
    else
      sprintf('Level %d, Snapshot #%d, %d deg.\nThis is snapshot %d of %d, checking %d verticies ... ',level,inumsort(i),rsort(i),i,length(inumsort),length(good))
      angi=find(Rc(1,5,:)==rsort(i)); 
      [Xi,Yi]=pred2([vert(check(good,1),1) vert(check(good,1),2) vert(check(good,1),3) ones(length(good),1)],Rc(:,1:4,angi),Par(1:8,angi));
      % [Xi,Yi]=snapcube(vert(check(good,1),1),vert(check(good,1),2),vert(check(good,1),3),pos,par,rsort(i));

      % Which points of Xi and Yi are inside the silhouette?
      vicinity=ones(length(good),1);
      vicinity(find(Xi<lims(inumsort(i),1) | Xi>lims(inumsort(i),2) | Yi<lims(inumsort(i),3) | Yi>lims(inumsort(i),4)))=0;
      check(good(find(~vicinity)),2)=0;
      vici=find(vicinity);
      Zi=interp2(silh(:,:,inumsort(i)),Xi(vici),Yi(vici));
      check(good(vici),2)=Zi;
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
