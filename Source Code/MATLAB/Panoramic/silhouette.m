function [bw3,Zi]=silhouette(a,Xi,Yi,level);

  bw1=a>level;
  bw1=bwmorph(bw1,'clean',Inf);
  bw1=bwmorph(bw1,'close',Inf);
  bw1label=bwlabel(bw1);

  stats=regionprops(bw1label,'all');
  areas=[stats.Area];
  blob=find(areas==max(areas));
  bw2=[stats(blob).FilledImage];
  lims=[stats(blob).BoundingBox];
  xl1=ceil(lims(1));
  yl1=ceil(lims(2));
  xl2=xl1+lims(3);
  yl2=yl1+lims(4);

  % Make the mask
  bw3=zeros(size(bw1));
  bw3(yl1:yl2-1,xl1:xl2-1)=bw2;

  %edg1=bwperim(bw3);
  %edg2=bwmorph(edg1, 'skel', 3);
  %[ex,ey]=find(edg2);

  %figure
  %subplot(1,3,1), subimage(a);
  %subplot(1,3,2), subimage(bw3);
  %subplot(1,3,3), subimage(a), hold on, plot(ey,ex,'.');

  % Which points of Xi and Yi are inside the silhouette?
  notvicinity=find(Xi<xl1 | Xi>xl2 | Yi<yl1 | Yi>yl2);
  Zi=interp2(bw3,Xi,Yi);
  Zi(notvicinity)=0.0;
  Zi(find(Zi>=0.5))=1.0;
  Zi(find(Zi<0.5))=0.0;
