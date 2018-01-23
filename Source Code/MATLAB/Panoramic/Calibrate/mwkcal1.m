% This evolved into mwkcal2.m

clear all
close all

%----------------------
% User defined parameters here
numpts=10;

%----------------------

done=0;
while ~done
  calfile=input('Calibration image filename: [cal.bmp] ','s');
  if isempty(calfile), calfile='cal.bmp'; end;
  aid=fopen(calfile);
  if aid~=-1
    fclose(aid);
    done=1;
  end
end
a=rgb2gray(imread(calfile));
fh=figure;
imshow(a);

% Identify 2 sides of the cube for easy determination of normal vectors
rioh=figure;
sprintf('Select side A')
sideA=roipoly(a);
[Ai]=find(sideA);
nA=input('Unit normal for side A [x,y,z]: ');

sprintf('Select side B')
sideB=roipoly(a);
[Bi]=find(sideB);
nB=input('Unit normal for side B [x,y,z]: ');
close(rioh);

xi=zeros(numpts,1);
yi=xi; X=xi; Y=xi; Z=xi; NX=xi; NY=xi; NZ=xi;

for i=1:numpts 
  sprintf('Click on calibration point #%d.',i)  
  [xi(i),yi(i),P]=impixel(a);
  text(xi(i),yi(i),sprintf('%d',i));
  xyi=sub2ind(size(a),xi(i),yi(i));
  isA=find(Ai==xyi);
  if isempty(isA)              % Could be on Side B
    isB=find(Bi==xyi);
    if isempty(isB)            % then neither Side A nor B
      sprintf('Point not on side A or B')
    end
  else                         % then on Side A
  
  end
  done=0;
  while ~done
    
    ttext=sprintf('Xi(%d): ',i);
    X(i)=input(ttext);
    ttext=sprintf('Yi(%d): ',i);
    Y(i)=input(ttext);
    ttext=sprintf('Zi(%d): ',i);
    Z(i)=input(ttext);
    ttext=sprintf('NX(%d): ',i);
    NX(i)=input(ttext);
    ttext=sprintf('NY(%d): ',i);
    NY(i)=input(ttext);
    ttext=sprintf('NZ(%d): ',i);
    NZ(i)=input(ttext);
    sstext=sprintf('Is [%2.2f, %2.2f, %2.2f], [%2.2f, %2.2f, %2.2f] correct? [Y]: ',X(i),Y(i),Z(i),NX(i),NY(i),NZ(i));
    ss=input(sstext,'s');
    if (ss=='' | ss=='Y' | ss=='y') done=1; end;
  end
  X(i)=X(i)*25.4/2;
  Y(i)=Y(i)*25.4/2;
  Z(i)=Z(i)*25.4/2;
end

data3d=[X,Y,Z,xi,yi,NX,NY,NZ];
[par,pos,iter,res,er,C]=cacal('watec',data3d);








return

% OLD STUFF DOWN HERE

% First point is world origin
% Just mark it for now
sprintf('Mark world origin (usually a corner of the cube)')
[c,r,P]=impixel(a);   % c and r are x,y coordinates; P is color triplet
figure(fh); hold on; plot(c,r,'ro');
