clear all
close all

% Get the silhouettes
% save silhs silh lims theta
if 0 
  theta=zeros(1,72);
  for j=1:72
    i=j+3;
    fname=sprintf('cube/cube_11_30_01/cube%03d.BMP',i);
    a=rgb2gray(imread(fname));
    level=70;
    if j==1
      silh=zeros(size(a,1),size(a,2),72);
      lims=zeros(72,4);
      [silh(:,:,1),lims(1,:)]=getsilh(a,level);
    else
      [silh(:,:,j),lims(j,:)]=getsilh(a,level);
    end
    theta(j)=(i-4)*5;
    j
  end  
else
  load silhs
end

% load camera calibration parameters
load cube/cube_11_30_01/pospar.dat
pos=pospar(1:6);
par=pospar(7:14);

% load octree basis
octree_basis

% Establish size and origin of initial carving cube 
% Keep it cubic!!!     
xdim=4.5*25.4;  % mm
ydim=xdim;      % mm
zdim=xdim;      % mm

x0=-0.5*xdim;
xn=0.5*xdim;
y0=-0.5*ydim;
yn=0.5*ydim;
z0=0.25*zdim;
zn=-0.75*zdim;

% Create the initial verticies matrix (verts) and octree matrix (octree)
% vert(:,1)=x
% vert(:,2)=y
% vert(:,3)=z
% vert(:,4)=v (0 for outside, 1 for inside, NaN if not checked
% octree is 3 dimensional: (Num of octrees) x (27 indicies denoting rows in vert) 
% x (9, where octree(:,:,1) contains vert indicies and octree(:,:,2:9) denotes voxel verticies
% octree(octree_num,:,1)=vert indicies
% octree(octree_num,:,2:9)=voxel verticies
% At each level the octree matrix is rebuilt and the vert matrix is expanded.

vert=octb*[(xn-x0) 0 0; 0 (yn-y0) 0; 0 0 (zn-z0)];
sshift=zeros(27,3);
sshift(:,1)=x0;
sshift(:,2)=y0;
sshift(:,3)=z0;
vert=vert+sshift;
vert(:,4)=NaN;

octree=zeros(1,27,9);
octree(1,:,1)=1:27;
octree(1,:,2:end)=octbi;

% verticies of voxel h-1 of octree t is (remember, voxels are from 2:9):
% vert(octree(t,find(octree(t,:,h)),1),:);

level=1;
delt(level)=(xn-x0)/2;     % assumes initial carving block is a cube
while delt(level)>=5       % mm
  check_num=length(find(isnan(vert(:,4))));
  check=ones(check_num,2);
  check(:,1)=find(isnan(vert(:,4)));
  for i=1:72
    good=find(check(:,2));
    [Xi,Yi]=snapcube(vert(good,1),vert(good,2),vert(good,3),pos,par,theta(i));
    xl1=lims(i,1);
    xl2=lims(i,2);
    yl1=lims(i,3);
    yl2=lims(i,4);

    % Which points of Xi and Yi are inside the silhouette?
    notvicinity=find(Xi<xl1 | Xi>xl2 | Yi<yl1 | Yi>yl2);
    Zi=interp2(silh(:,:,i),Xi,Yi);
    Zi(notvicinity)=0.0;
    Zi(find(Zi>=0.5))=1.0;
    Zi(find(Zi<0.5))=0.0;
    check(good(find(~Zi)),2)=0;
  end
  good=find(check(:,2));
  
  % keep up with the points that have been tested
  % vert(:,4) of 0 or 1 means that the vertex has been tested
  % vert(:,4) of NaN means that the vertex has not been tested
  vert(check(:,1),4)=check(:,2);
  
  % which voxels were inside AND outside the volume?
  % find them, then subdivide and repeat for each octree
  % to which voxels were those that were inside belong?
  % 'good' tells me which verticies were inside the volume
  % and which verticies were outside the volume of the verticies
  % which have not been checked.
  
  % Grab the vertex index from good.
  % Find which octrees contain that index
  
  for j=1:size(octree,1)   % rows of octree(:,:,1) - cycle through each octree
    for i=1:length(good)
      [ovi(i,1),ovi(i,2)]=find(octree(:,:,1)==good(i));
    end
  end
 
  % ovi- octree and vertex indicies
  % ovi(:,1) contains octree numbers for verticies inside the volume
  % ovi(:,2) contains octree vertex numbers for the verticies
  
  % with which octrees are we dealing?
  octs=unique(ovi(:,1));
  
  % now divide and conquer
  for i=1:length(octs)
    % get verticies for this octree
    vs=ovi(find(ovi(:,1)==octs(i)),2);    
    % which voxels have these verticies?
    os=reshape(octree(i,vs,2:9),length(vs),8);
    % find the voxel numbers
    voxelsi=find(max(os)>0);
    % is a voxel totally inside the volume?
    vflag=0;
    for j=1:length(voxelsi)
      if length(find(ismember([1 2 3 4 5 6 7 8],os(:,voxelsi(j)))))==8 
        voxelsi(j)=NaN; 
	vflag=1;   % We have dropped one or more!
      end
    end
    if vflag voxelsi=voxelsi(find(~isnan)); end % drop voxels totally inside the volume
    
    % we now have the voxel number of each voxel of the octree 
    % (stored in voxelsi) that is both inside and outside the volume
    % now subdivide each voxel, rebuild the 'octree' matrix 
    % and append to the 'vert' matrix
    
    for j=1:length(voxelsi)                   % loop through each voxel of interest in the current octree
      ooi=find(octree(i,:,voxelsi(j)+1)>0);   % find verticies of voxel of interest
      vviorder=octree(i,ooi,voxelsi(j)+1);    % vertex order of indicies into vert
      vvi=octree(i,ooi,1);                    % actual indicies into vert for each vertex
      vxyz=vert(vvi,1:3);
      [f,I]=sort(vviorder); 
      vxyz=vxyz(I,:);                         % sort by vertex number, easier to determine limits this way
      
      % now create an octree inside the voxel
      if i==1 & j==1
        new_octree=zeros(1,27,9).*NaN;
	new_octree(1,ooi,1)=vvi;          % these verticies already exist
	new_octree(1,:,2:end)=octbi;   IF THIS IS RIGHT THEN THE THIRD DIM OF OCTREE IS REDUNDANT!
      else
      
      end
      
      
    end
  
  end
  
  
  level=level+1;
end
delt(end)=[];




return

V(find(~v))=0;
fv=isosurface(X,Y,Z,V,0.9)
p=patch(fv);
isonormals(X,Y,Z,V,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

return

%---------------------------
% Old Stuff down here

for i=0:max(max(bw1label))
  bw1labeln(i+1)=length(find(bw1label==i));
end

biguns=(find(bw1labeln==max(bw1labeln))-1);
bw2=zeros(size(bw1));
bw2(find(bw1label==biguns))=1;
bw2=bwmorph(bw2,'majority',Inf);
bw2=bwmorph(bw2,'close',inf);
bw2=bwmorph(bw2,'majority',inf);

bw2label=bwlabel(bw2);
for i=0:max(max(bw2label))
  bw2labeln(i+1)=length(find(bw2label==i));
end


return

figure
bw2=bwmorph(bw1, 'skel', 3);
b=bwmorph(bw2, 'majority', 1);
subplot(1,3,1), subimage(a);
subplot(1,3,2), subimage(b);
edg1=bwperim(b);
edg2=bwmorph(edg1, 'skel', 3);
[ex,ey]=find(edg2);
subplot(1,3,3), subimage(b), hold on, plot(ey,ex,'.');



xsize=size(X);
ysize=size(Y);
zsize=size(Z);
x=reshape(X,xsize(1)*xsize(2)*xsize(3),1);
y=reshape(Y,ysize(1)*ysize(2)*ysize(3),1);
z=reshape(Z,zsize(1)*zsize(2)*zsize(3),1);

% move origin to center of cube
x=x-(xdim)/2;
y=y-(ydim)/2;
z=z-(zdim)/2;
% move the origin to the top of the cube
z=z+min(z);
% inches to mm
x=x.*25.4;
y=y.*25.4;
z=z.*25.4;
 
 
 
 
fv=isosurface(X,Y,Z,V,0.9)
p=patch(fv);
isonormals(X,Y,Z,V,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
