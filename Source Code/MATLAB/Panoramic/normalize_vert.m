% normalizes the vert matrix
vert(:,1)=vert(:,1)-min(vert(:,1));
vert(:,2)=vert(:,2)-min(vert(:,2));
vert(:,3)=vert(:,3)-min(vert(:,3));
vert(:,1)=vert(:,1)./max(vert(:,1));
vert(:,2)=vert(:,2)./max(vert(:,2));
vert(:,3)=vert(:,3)./max(vert(:,3));
