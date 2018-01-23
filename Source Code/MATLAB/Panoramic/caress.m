function [silh2,lims2]=caress(silh1,dtheta,centroids1,orientations1,thresharr)
% correct for essentricity in spinny image silhouettes

theta=0:dtheta:(360-dtheta);
hemi=360/(2*dtheta);

plts=1;
if plts
  figure
  set(gcf,'Position',[700,50,400,900])
end

silh2=zeros(size(silh1));
lims2=zeros(length(theta),4);
dc=zeros(hemi,2);
dtheta=zeros(1,hemi);

for i=1:hemi
  sprintf('Correcting silhouettes #%d and #%d ... ',i,i+hemi)
  [gx1,gy1]=find(silh1(:,:,i)>0);
  [gx2,gy2]=find(silh1(:,:,i+hemi)>0);
  dc(i,:)=(centroids1(i+hemi,:)-centroids1(i,:))./2;
  dtheta(i)=((abs(orientations1(i))-abs(orientations1(i+hemi)))./2)*(pi/180);
  
  if plts
    subplot(3,1,1)
    plot(gx1,gy1,'b.'); hold on;
    plot(gx2,gy2,'y.');
    title(sprintf('Original images #%d and #%d',i,i+hemi))
  end
  
  % Translate first
  
  if plts
    subplot(3,1,2)
    plot(gx1+dc(i,2),gy1+dc(i,1),'b.'); hold on;
    plot(gx2-dc(i,2),gy2-dc(i,1),'y.');
    title('Translate')
  end
  
  % Then rotate about the centroid
  
  w1=[cos(-dtheta(i))  sin(-dtheta(i));
    -sin(-dtheta(i))  cos(-dtheta(i))];
  w2=[cos(dtheta(i))  sin(dtheta(i));
    -sin(dtheta(i))  cos(dtheta(i))];
  
  gxy1=[gx1-centroids1(i,2),gy1-centroids1(i,1)]*w1 + [ones(size(gx1))*(dc(i,2)+centroids1(i,2)) ones(size(gx1))*(dc(i,1)+centroids1(i,1))];
  gxy2=[gx2-centroids1(i+hemi,2),gy2-centroids1(i+hemi,1)]*w2 + [ones(size(gx2))*(-dc(i,2)+centroids1(i+hemi,2)) ones(size(gx2))*(-dc(i,1)+centroids1(i+hemi,1))];
  
  if plts  
    subplot(3,1,3)
    plot(gxy1(:,1),gxy1(:,2),'b.'); hold on;
    plot(gxy2(:,1),gxy2(:,2),'y.');
    title('Rotate')
  end
  
  lims2(i,1)=min(floor(gxy1(:,2)));
  lims2(i,2)=max(floor(gxy1(:,2)));
  lims2(i,3)=min(floor(gxy1(:,1)));
  lims2(i,4)=max(floor(gxy1(:,1)));
  lims2(i+hemi,1)=min(floor(gxy2(:,2)));
  lims2(i+hemi,2)=max(floor(gxy2(:,2)));
  lims2(i+hemi,3)=min(floor(gxy2(:,1)));
  lims2(i+hemi,4)=max(floor(gxy2(:,1)));
  
  ind1=sub2ind([size(silh1,1) size(silh1,2)],floor(gxy1(:,1)),floor(gxy1(:,2)));
  ind2=sub2ind([size(silh1,1) size(silh1,2)],floor(gxy2(:,1)),floor(gxy2(:,2)));
  jj1=zeros(size(silh1,1),size(silh1,2));
  jj1(ind1)=1;
  jj2=zeros(size(silh1,1),size(silh1,2));
  jj2(ind2)=1;

  se=strel('square',2);
  silh2(:,:,i)=imclose(jj1,se);
  silh2(:,:,i+hemi)=imclose(jj2,se);
  
  if plts
    figure(2)
    subplot(1,2,1)
    imagesc(silh1(:,:,i));
    subplot(1,2,2)
    imagesc(silh2(:,:,i));
    pause
    close(2)
  end
end

save silhs2 silh2 lims2 thresharr dc dtheta
