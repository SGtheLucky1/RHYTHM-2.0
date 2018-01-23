ee=1;
dd=1;
plot3(vert(octree(ee,:),1),vert(octree(ee,:),2),vert(octree(ee,:),3),'o')
hold on
plot3(vert(new_octree(dd,:),1),vert(new_octree(dd,:),2),vert(new_octree(dd,:),3),'k*')

hh=[2 10 13 17 19 24 25 27];
plot3(vert(new_octree(dd,hh),1),vert(new_octree(dd,hh),2),vert(new_octree(dd,hh),3),'k*')


plot3(vert(:,1),vert(:,2),vert(:,3),'k*','MarkerSize',10);
gg=find(vert(:,4)==1);
hold on
plot3(vert(gg,1),vert(gg,2),vert(gg,3),'rs','MarkerSize',11);

plot3(vert(octree(:,:),1),vert(octree(:,:),2),vert(octree(:,:),3),'o')

for ii=1:size(ovi,1)
plot3(vert(octree(ovi(ii,1),ovi(ii,2)),1),vert(octree(ovi(ii,1),ovi(ii,2)),2),vert(octree(ovi(ii,1),ovi(ii,2)),3),'gs')
end

plot3(vert(octree(i,vs),1),vert(octree(i,vs),2),vert(octree(i,vs),3),'ks','MarkerSize',10)

plot3(vxyz(:,1),vxyz(:,2),vxyz(:,3),'rs','MarkerSize',12)

plot3(vert(octree(1,1),1),vert(octree(1,1),2),vert(octree(1,1),3),'ko','MarkerSize',12)
