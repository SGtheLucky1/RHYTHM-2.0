function [in,coeff] = in_plane(a,b,c,p)

%IN_PLANE Check if a point is element of a plane.
%   [IN,COEFF] = IN_PLANE(A,B,C,P) returns a flag IN with value 1, if
%   the point P is an element of the plane specified by the points A,B, C.
%   The value of the flag IN is 0 otherwise.
%   COEFF are the coordinates of the point P in the coordinate system
%   spanned by the vectors B-A and C-A.
%
%   The points A, B, C, P are given as 3-by-1 vectors.

in = 0; coeff = [];

if ( cond([b(1:2)-a(1:2),c(1:2)-a(1:2)]) < 1e6 )
    
    coeff = [b(1:2)-a(1:2),c(1:2)-a(1:2)]\[p(1:2)-a(1:2)];
    if ( abs(a(3)+coeff(1)*(b(3)-a(3))+coeff(2)*(c(3)-a(3))-p(3)) < 1e-10 )
        in = 1;
    else
        in = 0;
        coeff = [];
    end
    
elseif ( cond([b(2:3)-a(2:3),c(2:3)-a(2:3)]) < 1e6 )
    
    coeff = [b(2:3)-a(2:3),c(2:3)-a(2:3)]\[p(2:3)-a(2:3)];
    if ( abs(a(1)+coeff(1)*(b(1)-a(1))+coeff(2)*(c(1)-a(1))-p(1)) < 1e-10 )
        in = 1;
    else
        in = 0;
        coeff = [];
    end
    
elseif ( cond([b([1,3])-a([1,3]),c([1,3])-a([1,3])]) < 1e6 )
    
    coeff = [b([1,3])-a([1,3]),c([1,3])-a([1,3])]\[p([1,3])-a([1,3])];
    if ( abs(a(2)+coeff(1)*(b(2)-a(2))+coeff(2)*(c(2)-a(2))-p(2)) < 1e-10 )
        in = 1;
    else
        in = 0;
        coeff = [];
    end
    
end
