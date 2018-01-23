function draw_poly_intersection(k_i,p_i,X_1,X_2,k_1,k_2)

%DRAW_INTERSECTION Visualization of tetrahedron intersection.
%   DRAW_INTERSECTION(K,P_I,X_1,X_2,faces) returns three figures.
%   Figure 1 shows the opaque view of the two tetrahedrons.
%   Figure 2 shows the transparent view of the two tetrahedrons together with
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

for f = 1:3
    
    if f == 1
        k = k_1; p = X_1;
    elseif f == 2
        k = k_2; p = X_2;
    elseif f == 3
        k = k_i; p = p_i;
    end

    [num,tmp] = size(k);
    if num > 3
        for i = 1:num
            a = p(k(i,2),:)-p(k(i,1),:);
            b = p(k(i,3),:)-p(k(i,1),:);
            n(i,:) = cross(a,b); n(i,:) = n(i,:)/norm(n(i,:));
        end
    end

    i=1;
    while i <= num
   
        cr = cross(ones(num,1)*n(i,:),n,2);
        parallel = find(max(abs(cr),[],2) < 1e-10);
        parallel(parallel==i) = [];
        coord = [0,0; 1,0; 0,1]; same = 0;
        vert = p(k(i,:),:);
        if length(parallel) >= 1
            for  j = 1:length(parallel)
                for jj = 1:3
                    [in(jj),coeff] = in_plane(p(k(i,1),:)',p(k(i,2),:)',p(k(i,3),:)',p(k(parallel(j),jj),:)');
                    coord = [coord;coeff'];  
                    if in(jj) == 1
                        vert = [vert;p(k(parallel(j),jj),:)];
                    end
                end
                if any(in)
                    same(j) = 1;
                else
                    same(j) = 0;
                end
            end
            order = convhull(coord(:,1),coord(:,2));
            vert = vert(order,:);
        end
    
        if ( f == 1 )
            figure(1), hold on, axis equal
            fill3(vert(:,1),vert(:,2),vert(:,3),'b'); alpha(0.999),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
            figure(2), hold on
            fill3(vert(:,1),vert(:,2),vert(:,3),'b'); alpha(0.4),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
        elseif ( f == 2 )
            figure(1), hold on, axis equal
            fill3(vert(:,1),vert(:,2),vert(:,3),'g'); alpha(0.999),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
            figure(2), hold on
            fill3(vert(:,1),vert(:,2),vert(:,3),'g'); alpha(0.4),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
        elseif ( f == 3 )
            figure(2), hold on, axis equal
            fill3(vert(:,1),vert(:,2),vert(:,3),'r'); alpha(0.4),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
            figure(3), hold on
            fill3(vert(:,1),vert(:,2),vert(:,3),'r'); alpha(0.7),view(155,44);
            axis([minx,maxx,miny,maxy,minz,maxz]), axis off
        end
        
        k(parallel(same==1),:) = []; n(parallel(same==1),:)=[];
        [num,tmp] = size(k);
        i = i+1; clear in;
    
    end
    
    clear n;
    
end

opengl autoselect;