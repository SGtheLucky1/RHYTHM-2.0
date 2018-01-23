function [Xi,Yi]=snapcube(x,y,z,pos,par,theta)

  xp=x.*cos(theta*pi/180)-y.*sin(theta*pi/180);
  yp=x.*sin(theta*pi/180)+y.*cos(theta*pi/180);

  % take a snapshot of the cube
  [Xi,Yi]=pred([xp yp z],par,pos);




