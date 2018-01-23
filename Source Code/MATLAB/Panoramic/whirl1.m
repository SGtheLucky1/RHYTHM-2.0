clear all
close all
xdim=4.5*25.4;  % mm
ydim=4.5*25.4;  % mm
zdim=4.5*25.4;  % mm

% Make a cube where [X,Y,Z] define voxel corners
delt=20;
xsize=delt+1;
ysize=delt+1;
zsize=delt+1;
X=ones(xsize,ysize,zsize);
Y=X;
Z=X;
V=X;
[X,Y,Z]=meshgrid(-0.5*xdim:(xdim)/delt:0.5*xdim,-0.5*ydim:(ydim)/delt:0.5*ydim,0.25*zdim:-(zdim)/delt:-0.75*zdim);

x=reshape(X,size(Z,1)*size(Z,2)*size(Z,3),1);
y=reshape(Y,size(Z,1)*size(Z,2)*size(Z,3),1);
z=reshape(Z,size(Z,1)*size(Z,2)*size(Z,3),1);
v=reshape(V,size(Z,1)*size(Z,2)*size(Z,3),1);

% camera calibration parameters
load cube/cube_11_30_01/pospar.dat
pos=pospar(1:6);
par=pospar(7:14);

for i=4:75
  i
  calfile=sprintf('cube/cube_11_30_01/cube%03d.BMP',i);
  a=rgb2gray(imread(calfile));

  theta=(i-4)*5;
  
  good=find(v);
  [Xi,Yi]=snapcube(x(good),y(good),z(good),pos,par,theta);

  % get silhouette
  level=70;
  [bw3,Zi]=silhouette(a,Xi,Yi,level);
  ggood=find(Zi);
  notggood=find(~Zi);
  v(good(notggood))=0;
  
  if 0
    figure
    imshow(bw3); hold on;
    plot(Xi(ggood),Yi(ggood),'r.');
    pause
    close
  end
  
  i
end

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
