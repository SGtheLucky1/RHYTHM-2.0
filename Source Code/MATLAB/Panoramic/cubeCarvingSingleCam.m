function [] = cubeCarvingSingleCam(hdir,silh,lims,dtheta,n_images,dofrontback,r,...
    rr,irsort,inumsort,rsort)
% Description: A modification of the occluding contours portion of whirl.m
% written by Dr. Matthew Kay. The silhsProcess_callback was getting
% excessively long so I decided to consolidate this code into a function.
% ATTENION: You will need to modify lines 129, 139, 149, 159, 169 to point
% to the Panoramic directory which should be downloaded with this set of
% functions. The recommended location for the Panoramic directory is in the
% Matlab root directory.
%
% Inputs:
% hdir = the home directory
% silh = silhouettes
% lims = limits of the bounding box for the silhouettes
% dtheta = degree of each step
% n_images = total number of images used
% dofrontback = (1) collapse redudant images, (0) do not collapse
% r = image angles
% rr = list of angles being used and their index
% irsort = 
% inumsort =
% rsort = 
%
% Outputs:
%
% Author: Christopher Gloschat
% Date: June 23, 2016
%
% Modification Log:
%
%
%% Code %%

% Load camera calibrations
fprintf('\n Load camera parameters:\n');
% Find Rc.dat file
Rcfname = [hdir '/Rc.dat'];
fid = fopen(Rcfname);
if fid == -1
    errordlg(['Could not find Rc.dat! Place Rc.dat in home directory '...
        'and try again.']);
    return
else
    fclose(fid);
    disp('Found Rc.dat.')
    [Rc,Rc_ntheta,Rc_dtheta,message]=readRccal(Rcfname);
    if dtheta~=Rc_dtheta
        errordlg(['Designated value for theta in WHIRL does not match '...
            'the theta value in Rc.dat. Process aborted!']);
        return
    end
end
Parfname = [hdir '/Par.dat'];
fid = fopen(Parfname);
if fid == -1
    errordlg(['Could not find Par.dat! Place Par.dat in home directory '...
        'and try again.']);
    return
else
    fclose(fid);
    disp('Found Par.dat.')
    [Par,Par_ntheta,Par_dtheta,message]=readParcal(Parfname);
    if dtheta~=Par_dtheta
        errordlg(['Designated value for theta in WHIRL does not match '... 
            'the theta value in Par.dat. Process aborted!']);
        return
    end
end
disp('Camera parameters loaded.')
% Check camera calibration directionality
calibDirection = questdlg(['IF CAMERA CAILBRATION AND OBJECT SCANNING '...
    'DIRECTIONALITY (E.G. CW OR CCW) DO NOT MATCH THE CALIBRATION MUST '...
    'BE REORDERED. IS REORDERING NEEDED? [N]'],'Reorder Camera Calibration',...
    'Yes','No','No');
switch calibDirection
    case 'Yes'
        disp('Reordering camera parameters ....')
        Par(1:8,2:end)=fliplr(Par(1:8,2:end));
        Rc_new=Rc;
        for i=2:size(Rc,3)
            Rc_new(1:4,1:4,i)=Rc(1:4,1:4,size(Rc,3)-i+2);
        end
        Rc=Rc_new;
        clear Rc_new
    case 'No'
        disp('Ok, will not reorder camera parameters...')
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
    fprintf('Angle %d',Rc(1,5,i));
    figure(tt);
    if i==1, text(xyz(i,1),xyz(i,2),xyz(i,3),sprintf('%d',Rc(1,5,i)));
    else
        text(xyz(i,1),xyz(i,2),xyz(i,3),sprintf('%d',Rc(1,5,i)));
    end
    pause(0.1);
end
clear xyz;

%% Load octree basis %%
[octb,octbi] = octree_basis;

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

levs=[(1:10)' (Xn-X0)./(2.^(1:10)')];        % Level number, delt
disp([(1:10)' (Xn-X0)./(2.^(1:10)')])
startlevel=input('Start at level number: ');
maxlevel=input('Max level number: ');

fnametag=input('Filename tag (ie: _povcyl_1mm_corrected): ','s');
% savedir=input('Save datafiles in this directory, leave off the "/" (ie: povcyl_05mm): ','s');
savedir = [hdir '/Geometry/'];

% verticies of voxel h-1 of octree t is (remember, voxels are from 2:9):
% vert(octree(t,find(octree(t,:,h)),1),:);

level=startlevel;
delt(1,1)=(Xn-X0)/(2^level);
delt(1,2)=(Yn-Y0)/(2^level);
delt(1,3)=(Zn-Z0)/(2^level);

%% Record analysis parameters in textfile %%
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

%%
goodthresh=0.95;
% Rotation matrix for rotating the cube
Rz = [cosd(dtheta) -sind(dtheta) 0; sind(dtheta) cosd(dtheta) 0; 0 0 1];
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
  
  if (level==1 && level==startlevel) 
    fprintf('Creating octree for level 1...\n')
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
  elseif (level==2 && level==startlevel)
    fprintf('Loading level 2 octree from datafile...\n');
    load /Users/Chris/Documents/MATLAB/Panoramic/level2
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==3 && level==startlevel)
    fprintf('Loading level 3 octree from datafile...\n');
    load /Users/Chris/Documents/MATLAB/Panoramic/level3
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift   
  elseif (level==4 && level==startlevel)
    fprintf('Loading level 4 octree from datafile...\n');
    load /Users/Chris/Documents/MATLAB/Panoramic/level4
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==5 && level==startlevel)
    fprintf('Loading level 5 octree from datafile...\n');
    load /Users/Chris/Documents/MATLAB/Panoramic/level5
    vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
    sshift=zeros(size(vert,1),3);
    sshift(:,1)=X0;
    sshift(:,2)=Y0;
    sshift(:,3)=Z0;
    vert(:,1:3)=vert(:,1:3)+sshift;  
    clear sshift 
  elseif (level==6 && level==startlevel)
    fprintf('Loading level 6 octree from datafile...\n');
    load /Users/Chris/Documents/MATLAB/Panoramic/level6
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
	if i==1 && j==1
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
  
  sites2check=find(isnan(vert(:,4))); % indices into vert
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
% %       fprintf('Level %d, Snapshot #%d, %d deg.\nThis is snapshot %d of %d, checking %d verticies ... ',level,inumsort(i),rsort(i),i,length(inumsort),length(good));
% % %       angi=find(Rc(1,5,:)==rsort(i));
% %       fprintf('  Applying camera model ...');
% %       [Xi,Yi]=pred2([vert(check(good,1),1) vert(check(good,1),2) vert(check(good,1),3) ones(length(good),1)],Rc(:,1:4,1),Par(1:8,1));
% %       % [Xi,Yi]=snapcube(vert(check(good,1),1),vert(check(good,1),2),vert(check(good,1),3),pos,par,rsort(i));
% %       fprintf('  Determining which points are inside the volume ...');
% %       % Which points of Xi and Yi are inside the silhouette?
% %       vicinity=ones(length(good),1);
% %       vicinity(Xi<lims(inumsort(i),1) | Xi>lims(inumsort(i),2) | Yi<lims(inumsort(i),3) | Yi>lims(inumsort(i),4))=0; % finds voxels outside the silh bounding box limits
% %       check(good(~vicinity),2)=0;
% %       vici=find(vicinity);
% %       Zi=interp2(silh(:,:,inumsort(i)),Xi(vici),Yi(vici));
% %       check(good(vici),2)=Zi;
% %       fprintf('  Done with snapshot #%d',inumsort(i));

fprintf('Level %d, Snapshot #%d, %d deg.\nThis is snapshot %d of %d, checking %d verticies ... ',level,inumsort(i),rsort(i),i,length(inumsort),length(good));
angi=find(Rc(1,5,:)==rsort(i));
fprintf('  Applying camera model ...');
[Xi,Yi]=pred2([vert(check(good,1),1) vert(check(good,1),2) vert(check(good,1),3) ones(length(good),1)],Rc(:,1:4,angi),Par(1:8,angi));
% [Xi,Yi]=snapcube(vert(check(good,1),1),vert(check(good,1),2),vert(check(good,1),3),pos,par,rsort(i));
fprintf('  Determining which points are inside the volume ...');
% Which points of Xi and Yi are inside the silhouette?
vicinity=ones(length(good),1);
vicinity(find(Xi<lims(inumsort(i),1) | Xi>lims(inumsort(i),2) | Yi<lims(inumsort(i),3) | Yi>lims(inumsort(i),4)))=0;
check(good(find(~vicinity)),2)=0;
vici=find(vicinity);
Zi=interp2(silh(:,:,inumsort(i)),Xi(vici),Yi(vici));
check(good(vici),2)=Zi;
fprintf('  Done with snapshot #%d',inumsort(i));
% Rotate cube
tmp = Rz*vert(:,1:3)';
vert(:,1:3) = tmp';

% %       if 0
% %         if i==1
% %             hh=figure; 
% %         end
% %         figure(hh);
% %         imagesc(silh(:,:,inumsort(i))); hold on;
% %         plot(Xi,Yi,'y.','MarkerSize',2);
% %         plot(Xi(vici(find(Zi)>0)),Yi(vici(find(Zi)>0)),'wo','MarkerSize',2);
% %         hold off;	
% %         pause;
% %         if i==length(inumsort) 
% %             close(hh);
% %         end
% %       end
      
    end
  end
  vert(check(:,1),4)=check(:,2);
  
  % keep up with the points that have been tested
  % vert(:,4) of ~NaN means that the vertex has been tested
  % vert(:,4) of NaN means that the vertex has not been tested
  
  level=level+1
end % End level loop

%% Normalization %%
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

%% Rendering %%
disp('Rendering ....')

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
figure
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


end