function [SA]=cellSA3d(V1,V2,V3);
% [SA]=3DcellSA(V1,V2,V3);
%
% Computes the surface area of 3D triangles defined by arrays V1, V2, and V3.
% V1, V2, V3 contain the coordinates of the verticies in 3 column format.
% Triangle centroids are computed and subtracted from the verticies to
% improve numerical percision.
%
% 2004, MWKay

% for better numerical precision shift origin to centroid of each cell
centroids=[(V1(:,1)+V2(:,1)+V3(:,1))./3 (V1(:,2)+V2(:,2)+V3(:,2))./3 (V1(:,3)+V2(:,3)+V3(:,3))./3]; 
V1n=V1-centroids;
V2n=V2-centroids;
V3n=V3-centroids;
x1=V1n(:,1); y1=V1n(:,2); z1=V1n(:,3);
x2=V2n(:,1); y2=V2n(:,2); z2=V2n(:,3);
x3=V3n(:,1); y3=V3n(:,2); z3=V3n(:,3);
t1=y1.*(z2-z3)-z1.*(y2-y3)+y2.*z3-y3.*z2;
t2=z1.*(x2-x3)-x1.*(z2-z3)+z2.*x3-z3.*x2;
t3=x1.*(y2-y3)-y1.*(x2-x3)+x2.*y3-x3.*y2;
SA=0.5.*sqrt(t1.*t1+t2.*t2+t3.*t3);
