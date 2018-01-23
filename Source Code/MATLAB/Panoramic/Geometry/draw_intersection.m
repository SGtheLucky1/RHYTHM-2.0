function draw_intersection(k,p_i,X_1,X_2,faces)

%DRAW_INTERSECTION Visualization of tetrahedron intersection.
%   DRAW_INTERSECTION(K,P_I,X_1,X_2,faces) returns three figures.
%   Figure 1 shows the opaque view of the two tetrahedra.
%   Figure 2 shows the transparent view of the two tetrahedra together with
%            their intersection polytope.
%   Figure 3 shows the transparent intersection polytope seperately.
%
%   NOTE: To ensure, that the 3D view is shown correctly use the command
%         'opengl neverselect', before using Figure 1 and use the command
%         'opengl autoselect', before using Figure 2 and Figure 3.
%

maxx = max( [X_1(:,1);X_2(:,1)] );
maxy = max( [X_1(:,2);X_2(:,2)] );
maxz = max( [X_1(:,3);X_2(:,3)] );
minx = min( [X_1(:,1);X_2(:,1)] );
miny = min( [X_1(:,2);X_2(:,2)] );
minz = min( [X_1(:,3);X_2(:,3)] );

figure(1), opengl neverselect, hold on
fill3(X_1(faces)',X_1(faces+4)',X_1(faces+8)','g');
fill3(X_2(faces)',X_2(faces+4)',X_2(faces+8)','b'); 
alpha(0.999), view(60,20), axis([minx,maxx,miny,maxy,minz,maxz]), axis off

figure(2), opengl autoselect, hold on
fill3(X_1(faces)',X_1(faces+4)',X_1(faces+8)','g'); 
fill3(X_2(faces)',X_2(faces+4)',X_2(faces+8)','b');
alpha(0.5), view(60,20), axis([minx,maxx,miny,maxy,minz,maxz]), axis off

[num,tmp] = size(k);
if num > 3
    for i = 1:num
        a = p_i(k(i,2),:)-p_i(k(i,1),:);
        b = p_i(k(i,3),:)-p_i(k(i,1),:);
        n(i,:) = cross(a,b); n(i,:) = n(i,:)/norm(n(i,:));
    end
end
i=1;
while i <= num
    
    cr = cross(ones(num,1)*n(i,:),n,2);
    parallel = find(max(abs(cr),[],2) < 1e-12);
    coord = [0,0; 1,0; 0,1]; in = [0 0];
    if length(parallel) > 1
        for  j = 2:length(parallel)
            for jj = 1:3
                [in(j),coeff] = in_plane(p_i(k(i,1),:)',p_i(k(i,2),:)',p_i(k(i,3),:)',p_i(k(parallel(j),jj),:)');
                coord = [coord;coeff'];
            end
        end
        order = convhull(coord(:,1),coord(:,2));
        vert = [p_i(k(i,:),:);p_i(k(parallel(find(in==1)),:)',:)];
        vert = vert(order,:);
    else
        vert = p_i(k(i,:),:);
    end
    
    figure(2)
    fill3(vert(:,1),vert(:,2),vert(:,3),'r'); alpha(0.3),view(60,20);
    axis([minx,maxx,miny,maxy,minz,maxz]), axis off
    figure(3), hold on
    fill3(vert(:,1),vert(:,2),vert(:,3),'r'); alpha(0.6),view(60,20);
    axis([minx,maxx,miny,maxy,minz,maxz]), axis off
    k(parallel(in==1),:) = []; n(parallel(in==1),:)=[];
    [num,tmp] = size(k);
    i = i+1; clear in;
    
end
