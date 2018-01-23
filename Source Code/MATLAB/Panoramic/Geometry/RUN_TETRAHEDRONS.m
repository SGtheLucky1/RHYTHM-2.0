disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%   DRIVER PROGAMME FOR TETRAHEDRON INTERSECTION      %%')
disp('    %%                                                     %%')
disp('    %%   (c) Martin Kaeser, TU Muenchen, 11.11.2002        %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

clear, close all;

data_tetrahedrons

% Coordinates of Tetraeders
disp(sprintf('\n COORDINATES OF TETRAHEDRON 1:')), disp(X_1) 
disp(sprintf(' COORDINATES OF TETRAHEDRON 2:')), disp(X_2) 

[k_1,V_1] = convhulln(X_1);
[k_2,V_2] = convhulln(X_2);

edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
faces = [1 2 3; 2 3 4; 3 4 1; 4 1 2];

% Computation of intersection points
p_i = [];

for i=1:6
    e = X_1(edges(i,:),:)';
    for j=1:4
        f = X_2(faces(j,:),:)';
        p_i = [ p_i; intersect_plane_line(f(:,1),f(:,2),f(:,3),e(:,1),e(:,2)) ];
    end
end

for i=1:6
    e = X_2(edges(i,:),:)';
    for j=1:4
        f = X_1(faces(j,:),:)';
        p_i = [ p_i; intersect_plane_line(f(:,1),f(:,2),f(:,3),e(:,1),e(:,2)) ];
    end
end
        
p_i = [ p_i; X_2(find(tsearchn(X_1,[1 2 3 4],X_2)==1),:) ];  
p_i = [ p_i; X_1(find(tsearchn(X_2,[1 2 3 4],X_1)==1),:) ];  

p_i = sortrows(p_i); p_i(find(max(abs(diff(p_i)),[],2)<1e-10),:) = [];

[pl,tmp] = size(p_i);
if pl > 3
    [k_i,V_i] = convhulln(p_i);
else
    k_i = []; V_i = 0;
end

disp(sprintf('\n\t VOLUME OF TETRAHEDRON 1:\t\t%8.6f',V_1));
disp(sprintf('\n\t VOLUME OF TETRAHEDRON 2:\t\t%8.6f',V_2));
disp(sprintf('\n\t VOLUME OF INTERSECTION POLYTOPE:\t%8.6f\n\n',V_i));
disp(sprintf('\n COORDINATES OF INTERSECTION POLYTOPE:')), disp(p_i)
draw_intersection(k_i,p_i,X_1,X_2,faces);
