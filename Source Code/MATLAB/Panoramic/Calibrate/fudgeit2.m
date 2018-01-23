clear all
close all

% Calibration points defined in file.
% 3D positions in reference to a fixed right-handed frame!
%

%----------------------
% User defined values here
zthresh=90;
zrange=60;
method=2;            % 1 for click center, 2 for line intersect
%----------------------

cw=0;
while (cw~=1 & cw~=2)
  cw=input('Calibrated cw or ccw? [1:cw, 2:ccw]: ');
end
if cw==1, cwtxt='cw'; end;
if cw==2, cwtxt='ccw'; end;
ang=pi;
while (ang~=0 & ang~=90 & ang~=180 & ang~=270)
  ang=input('Calibration angle (ie: 0, 90, 180, or 270): ');
end
defaultfile=sprintf('calpts1_%d%s.txt',ang,cwtxt);
eval(sprintf('calptsfile=input(''Calibration points filename: [%s] '',''s'');',defaultfile));
if isempty(calptsfile), calptsfile=defaultfile; end;
calptsfile2=sprintf('/home/mwk/matlab/Panoramic/Calibrate/calpts/%s',calptsfile);
fid=fopen(calptsfile2);
if fid~=-1
  fclose(fid);
  done=1;
  disp(sprintf('Loading %s',calptsfile2));
else
  disp(sprintf('Could not load %s',calptsfile2));
  return
end

calpts_nondim=load(calptsfile2);       % nondimensional calibration point positions in cube basis

done=0;
while ~done
  if ang==0, anum=0;
  elseif ang==90, anum=2;
  elseif ang==180, anum=4;
  elseif ang==270, anum=6; end;
  defaultfile=sprintf('cal00%d.BMP',anum);
  eval(sprintf('calfile=input(''Calibration image filename: [%s] '',''s'');',defaultfile));
  if isempty(calfile), calfile=defaultfile; end;
  fid=fopen(calfile);
  if fid~=-1
    fclose(fid);
    done=1;
  end
end
disp(sprintf('Loading %s',calfile));
a=imread(calfile);
if size(a,3)>1, a=rgb2gray(a); end;
data=crosseyed(a,calpts_nondim,calfile,method,zrange,zthresh);  % THE MEAT OF THE FRUIT
data(:,1:3)=data(:,1:3).*(25.4/2);       % convert positions to mm
datacom=sprintf('data%d=data;',aa);
eval(datacom);

%-----------------------------------------
% Calibrate

disp(sprintf('\n Enter number for camera: \n'));
disp(sprintf('    1) watec_and_snappy \n'));
disp(sprintf('    2) interpolated_dalsa \n'));
done=0;
while ~done
  camnum=input('Camera number [1]: ');
  if isempty(camnum), camnum=1; end;
  switch camnum
    case 1
      camera='watec_and_snappy';
      done=1;
    case 2
      camera='interpolated_dalsa';
      done=1;
  end
end

[par,pos,iter,res,er,C,success]=cacal(camera,data1);
savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 camera',calfile); 
[Xi,Yi]=pred(data1(:,1:3),par,pos,camera); 
eval(savecommand);

% Save calibration parameters
posparfname=sprintf('%s.pospar',calfile);
posparfid=fopen(posparfname,'w');
for i=1:6
  fprintf(posparfid,'%15e\n',pos(i));
end
for i=1:8
  fprintf(posparfid,'%15e\n',par(i));
end
fclose(posparfid);

close all
imshow(a)
image(a)
hold on
xi=data1(:,4);
yi=data1(:,5);
plot(xi+0.5,yi+0.5,'go','MarkerFaceColor','g');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
plot(Xi+0.5,Yi+0.5,'ro','MarkerFaceColor','r');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
title('Predicted calibration points','FontSize',14);
jpegcommand=sprintf('print -djpeg %s.jpg',calfile);
eval(jpegcommand);

