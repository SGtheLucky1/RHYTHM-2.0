function [bw3,lims,centroid,orientation,area]=getsilh(a,level,aabb)

  warning off 
  
  if aabb 
    bw1=a>level;
  else
    bw1=a<level;
  end
  
  bw1=bwmorph(bw1,'clean',Inf);
  bw1=bwmorph(bw1,'close',Inf);
  bw1label=bwlabel(bw1);

  stats=regionprops(bw1label,'all');  % regionprops also gives a measure of essentricity
  areas=[stats.Area];
  blob=find(areas==max(areas));
  bw2=[stats(blob).FilledImage];
  lims=[stats(blob).BoundingBox];
  centroid=[stats(blob).Centroid];
  area=[stats(blob).Area];
  orientation=[stats(blob).Orientation];
  xl1=ceil(lims(1));
  yl1=ceil(lims(2));
  xl2=xl1+lims(3);
  yl2=yl1+lims(4);

  % Make the silhouette
  bw3=zeros(size(bw1));
  bw3(yl1:yl2-1,xl1:xl2-1)=bw2;

  % limits
  lims=[xl1 xl2 yl1 yl2];

  warning on 
  
