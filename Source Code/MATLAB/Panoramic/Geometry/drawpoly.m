function [p,q,n,m] = drawpoly(dummy)

%DRAWPOLY Draw two polygons.
%   [P,Q,N,M] = DRAWPOLY(dummy) returns an N-by-2 matrix P and 
%   and M-by-2 matrix Q of coordinates drspecifying two polygons.
%
%   Verticies of each polygon can be added with the left
%   mouse-button. To complete a polygon, the last vertex 
%   has to be added with the right mouse-button.
%
%   Each polygon has to consist of at least 3 non-collinear vertices.
%

figure,axis([-1 1 -1 1]),axis square
title('DRAWPOLY: Draw polygon P')

% get first polygon
i=1; fine=1;
while fine ~= 3
  [p(i,1),p(i,2),b(i)]=ginput(1);
  if i==1
      plot(p(i,1),p(i,2),'g*'),axis([-1 1 -1 1]),axis square,hold on
      title('DRAWPOLY: Draw polygon P')
  else
      plot(p(i,1),p(i,2),'g*')
      plot([p(i-1,1);p(i,1)],[p(i-1,2);p(i,2)],'g-');    
  end
  b(1:2)=1; fine=b(i);
  i=i+1;
end
n = length(p); 
clf,fill(p(:,1),p(:,2),'g'),axis([-1 1 -1 1]),axis square,hold on
title('DRAWPOLY: Draw polygon Q')


%get second polygon
i=1; fine=1;
while fine ~= 3
  [q(i,1),q(i,2),b(i)]=ginput(1);
  if i==1
      plot(q(i,1),q(i,2),'b*');
  else
      plot(q(i,1),q(i,2),'b*')
      plot([q(i-1,1);q(i,1)],[q(i-1,2);q(i,2)],'b-');    
  end
  b(1:2)=1; fine=b(i);
  i=i+1;
end
m = length(q);

clf
fill(p(:,1),p(:,2),'g'),hold on
fill(q(:,1),q(:,2),'b')
plot([p(:,1);p(1,1)],[p(:,2);p(1,2)],'k--'),axis([-1 1 -1 1])
grid on, axis square
title('DRAWPOLY: SPECIFIED POLYGONS')



