disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                      %%')
disp('    %%   DRIVER PROGAMME FOR CONVEX POLYTOPE INTERSECTION   %%')
disp('    %%                                                      %%')
disp('    %%   (c) Martin Kaeser, TU Muenchen, 18.11.2002         %%')
disp('    %%                                                      %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

clear, close all;

data_polytopes

% Coordinates of Polytopes
disp(sprintf('\n COORDINATES OF POLYTOPE 1:')), disp(X_1) 
disp(sprintf(' COORDINATES OF POLYTOPE 2:')), disp(X_2) 

[k_1,V_1] = convhulln(X_1);  Tes_1 = delaunayn(X_1);
[k_2,V_2] = convhulln(X_2);  Tes_2 = delaunayn(X_2);

edges_1 = unique(sort([k_1(:,1:2);k_1(:,2:3)],2),'rows');
edges_2 = unique(sort([k_2(:,1:2);k_2(:,2:3)],2),'rows');

% Computation of intersection points
p_i = [];

for i=1:size(edges_1,1)
    e = X_1(edges_1(i,:),:)';
    for j=1:size(k_2,1)
        f = X_2(k_2(j,:),:)';
        p_i = [ p_i; intersect_plane_line(f(:,1),f(:,2),f(:,3),e(:,1),e(:,2)) ];
    end
end

for i=1:size(edges_2,1)
    e = X_2(edges_2(i,:),:)';
    for j=1:size(k_1,1)
        f = X_1(k_1(j,:),:)';
        p_i = [ p_i; intersect_plane_line(f(:,1),f(:,2),f(:,3),e(:,1),e(:,2)) ];
    end
end
        
t = tsearchn(X_1,Tes_1,X_2); p_i = [ p_i; X_2(find(t./t==1),:) ];  
t = tsearchn(X_2,Tes_2,X_1); p_i = [ p_i; X_1(find(t./t==1),:) ];  

p_i = sortrows(p_i); p_i(find(max(abs(diff(p_i)),[],2)<1e-10),:) = [];

[pl,tmp] = size(p_i);
if pl > 3
    [k_i,V_i] = convhulln(p_i);
else
    k_i = []; V_i = 0;
end

disp(sprintf('\n\t VOLUME OF POLYTOPE 1:\t\t%8.6f',V_1));
disp(sprintf('\n\t VOLUME OF POLYTOPE 2:\t\t%8.6f',V_2));
disp(sprintf('\n\t VOLUME OF INTERSECTION POLYTOPE:\t%8.6f\n\n',V_i));
disp(sprintf('\n COORDINATES OF INTERSECTION POLYTOPE:')), disp(p_i)

draw_poly_intersection(k_i,p_i,X_1,X_2,k_1,k_2);