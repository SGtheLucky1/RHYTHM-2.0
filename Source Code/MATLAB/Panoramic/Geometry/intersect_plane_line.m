function p_i = intersect_plane_line(a,b,c,d,e)

%INTERSECT_PLANE_LINE Intersection of a plane and a line segment in 3D.
%   P_I = INTERSECT_PLANE_LINE(A,B,C,D,E) returns the coordinates of the 
%   point of intersection of a triangular section of a plane specified 
%   by the three corner points A, B, C and a line segment specified by the 
%   two end points D, E.
%
%   The points A, B, C, D, E are given as 3-by-1 vectors.

p_i = [];

n = cross(b-a,c-a);
if ( abs(dot(e-d,n)) > eps ) 

    M = [b-a,c-a,d-e];
    f = d-a;

    s = M\f;

    if ( (s(3)>-eps) & (s(3)<1+eps) )
        if ( (s(1)>-eps) & (s(2)>-eps) & (s(1)+s(2)<1+eps) )
            p_i = d' + s(3)*(e-d)';
        end
    end

end