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

imgnum=input('Number of calibration images (4 max) [1]: ');
if isempty(imgnum), imgnum=1; end;
if imgnum>4 
  sprintf('Can only process a max of 4 images!')
  return
end

for aa=1:imgnum
  done=0;
  while ~done
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
    defaultfile=sprintf('calpts%d%s.txt',ang,cwtxt);
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
    end
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
end
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

if imgnum==1 
  [par,pos,iter,res,er,C,success]=cacal(camera,data1);
  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 camera',calfile); 
  [Xi,Yi]=pred(data1(:,1:3),par,pos,camera); 
       
elseif imgnum==2 
  [par,pos,iter,res,er,C,success]=cacal(camera,data1,data2);
  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 data2 camera',calfile); 

elseif imgnum==3 
  [par1,pos1,iter,res,er,C,success]=cacal(camera,data1,data2,data3);
  Rerr1=NaN;
  if success
    [Xi11,Yi11]=pred(data1(:,1:3),par1,pos1(:,1),camera);
    [Xi12,Yi12]=pred(data2(:,1:3),par1,pos1(:,2),camera);
    [Xi13,Yi13]=pred(data3(:,1:3),par1,pos1(:,3),camera);
    Xerr1=[Xi11-data1(:,4);Xi12-data2(:,4);Xi13-data3(:,4)];
    Yerr1=[Yi11-data1(:,5);Yi12-data2(:,5);Yi13-data3(:,5)];
    Rerr1=sqrt(mean(Xerr1.^2 + Yerr1.^2));
  end
  
  [par2,pos2,iter,res,er,C,success]=cacal(camera,data1,data2);
  Rerr2=NaN;
  if success
    [Xi21,Yi21]=pred(data1(:,1:3),par2,pos2(:,1),camera);
    [Xi22,Yi22]=pred(data2(:,1:3),par2,pos2(:,2),camera);
    Xerr2=[Xi21-data1(:,4);Xi22-data2(:,4)];
    Yerr2=[Yi21-data1(:,5);Yi22-data2(:,5)];
    Rerr2=sqrt(mean(Xerr2.^2 + Yerr2.^2));
  end
  
  [par3,pos3,iter,res,er,C,success]=cacal(camera,data1,data3);
  Rerr3=NaN;
  if success
    [Xi31,Yi31]=pred(data1(:,1:3),par3,pos3(:,1),camera);
    [Xi32,Yi32]=pred(data3(:,1:3),par3,pos3(:,2),camera);
    Xerr3=[Xi31-data1(:,4);Xi32-data3(:,4)];
    Yerr3=[Yi31-data1(:,5);Yi32-data3(:,5)];
    Rerr3=sqrt(mean(Xerr3.^2 + Yerr3.^2));
  end
  
  [par4,pos4,iter,res,er,C,success]=cacal(camera,data2,data3);
  Rerr4=NaN;
  if success
    [Xi41,Yi41]=pred(data2(:,1:3),par4,pos4(:,1),camera);
    [Xi42,Yi42]=pred(data3(:,1:3),par4,pos4(:,2),camera);
    Xerr4=[Xi41-data2(:,4);Xi42-data3(:,4)];
    Yerr4=[Yi41-data2(:,5);Yi42-data3(:,5)];
    Rerr4=sqrt(mean(Xerr4.^2 + Yerr4.^2));
  end
  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 data2 data3 camera',calfile); 
  
  Rerr=[Rerr1,Rerr2,Rerr3,Rerr4];
  minRerr=find(Rerr==min(Rerr));
  sprintf('Min Rerr achieved with case %d',minRerr)
  eval(sprintf('pos=pos%s; par=par%d;',minRerr,minRerr));
  
elseif imgnum==4 
  [par1,pos1,iter,res,er,C,success]=cacal(camera,data1,data2,data3);
  Rerr1=NaN;
  if success
    [Xi11,Yi11]=pred(data1(:,1:3),par1,pos1(:,1),camera);
    [Xi12,Yi12]=pred(data2(:,1:3),par1,pos1(:,2),camera);
    [Xi13,Yi13]=pred(data3(:,1:3),par1,pos1(:,3),camera);
    Xerr1=[Xi11-data1(:,4);Xi12-data2(:,4);Xi13-data3(:,4)];
    Yerr1=[Yi11-data1(:,5);Yi12-data2(:,5);Yi13-data3(:,5)];
    Rerr1=sqrt(mean(Xerr1.^2 + Yerr1.^2));
  end
  
  [par2,pos2,iter,res,er,C,success]=cacal(camera,data1,data2);
  Rerr2=NaN;
  if success
    [Xi21,Yi21]=pred(data1(:,1:3),par2,pos2(:,1),camera);
    [Xi22,Yi22]=pred(data2(:,1:3),par2,pos2(:,2),camera);
    Xerr2=[Xi21-data1(:,4);Xi22-data2(:,4)];
    Yerr2=[Yi21-data1(:,5);Yi22-data2(:,5)];
    Rerr2=sqrt(mean(Xerr2.^2 + Yerr2.^2));
  end
  
  [par3,pos3,iter,res,er,C,success]=cacal(camera,data1,data3);
  Rerr3=NaN;
  if success
    [Xi31,Yi31]=pred(data1(:,1:3),par3,pos3(:,1),camera);
    [Xi32,Yi32]=pred(data3(:,1:3),par3,pos3(:,2),camera);
    Xerr3=[Xi31-data1(:,4);Xi32-data3(:,4)];
    Yerr3=[Yi31-data1(:,5);Yi32-data3(:,5)];
    Rerr3=sqrt(mean(Xerr3.^2 + Yerr3.^2));
  end
  
  [par4,pos4,iter,res,er,C,success]=cacal(camera,data2,data3);
  Rerr4=NaN;
  if success
    [Xi41,Yi41]=pred(data2(:,1:3),par4,pos4(:,1),camera);
    [Xi42,Yi42]=pred(data3(:,1:3),par4,pos4(:,2),camera);
    Xerr4=[Xi41-data2(:,4);Xi42-data3(:,4)];
    Yerr4=[Yi41-data2(:,5);Yi42-data3(:,5)];
    Rerr4=sqrt(mean(Xerr4.^2 + Yerr4.^2));
  end
  
  [par5,pos5,iter,res,er,C,success]=cacal(camera,data1,data3);
  Rerr5=NaN;
  if success
    [Xi51,Yi51]=pred(data1(:,1:3),par5,pos5(:,1),camera);
    [Xi52,Yi52]=pred(data3(:,1:3),par5,pos5(:,2),camera);
    Xerr5=[Xi51-data1(:,4);Xi52-data3(:,4)];
    Yerr5=[Yi51-data1(:,5);Yi52-data3(:,5)];
    Rerr5=sqrt(mean(Xerr5.^2 + Yerr5.^2));
  end
  
  [par6,pos6,iter,res,er,C,success]=cacal(camera,data2,data4);
  Rerr6=NaN;
  if success
    [Xi61,Yi61]=pred(data2(:,1:3),par6,pos6(:,1),camera);
    [Xi62,Yi62]=pred(data4(:,1:3),par6,pos6(:,2),camera);
    Xerr6=[Xi61-data2(:,4);Xi62-data4(:,4)];
    Yerr6=[Yi61-data2(:,5);Yi62-data4(:,5)];
    Rerr6=sqrt(mean(Xerr6.^2 + Yerr6.^2));
  end
  
  [par7,pos7,iter,res,er,C,success]=cacal(camera,data3,data4);
  Rerr7=NaN;
  if success
    [Xi71,Yi71]=pred(data3(:,1:3),par7,pos7(:,1),camera);
    [Xi72,Yi72]=pred(data4(:,1:3),par7,pos7(:,2),camera);
    Xerr7=[Xi71-data3(:,4);Xi72-data4(:,4)];
    Yerr7=[Yi71-data3(:,5);Yi72-data4(:,5)];
    Rerr7=sqrt(mean(Xerr7.^2 + Yerr7.^2));
  end
  
  [par8,pos8,iter,res,er,C,success]=cacal(camera,data2,data3,data4);
  Rerr8=NaN;
  if success
    [Xi81,Yi81]=pred(data2(:,1:3),par8,pos8(:,1),camera);
    [Xi82,Yi82]=pred(data3(:,1:3),par8,pos8(:,2),camera);
    [Xi83,Yi83]=pred(data4(:,1:3),par8,pos8(:,3),camera);
    Xerr8=[Xi81-data2(:,4);Xi82-data3(:,4);Xi83-data4(:,4)];
    Yerr8=[Yi81-data2(:,5);Yi82-data3(:,5);Yi83-data4(:,5)];
    Rerr8=sqrt(mean(Xerr8.^2 + Yerr8.^2));
  end
  
  [par9,pos9,iter,res,er,C,success]=cacal(camera,data1,data2,data4);
  Rerr9=NaN;
  if success
    [Xi91,Yi91]=pred(data1(:,1:3),par9,pos9(:,1),camera);
    [Xi92,Yi92]=pred(data2(:,1:3),par9,pos9(:,2),camera);
    [Xi93,Yi93]=pred(data4(:,1:3),par9,pos9(:,3),camera);
    Xerr9=[Xi91-data1(:,4);Xi92-data2(:,4);Xi93-data4(:,4)];
    Yerr9=[Yi91-data1(:,5);Yi92-data2(:,5);Yi93-data4(:,5)];
    Rerr9=sqrt(mean(Xerr9.^2 + Yerr9.^2));
  end
  
  [par10,pos10,iter,res,er,C,success]=cacal(camera,data1,data3,data4);
  Rerr10=NaN;
  if success
    [Xi101,Yi101]=pred(data1(:,1:3),par10,pos10(:,1),camera);
    [Xi102,Yi102]=pred(data3(:,1:3),par10,pos10(:,2),camera);
    [Xi103,Yi103]=pred(data4(:,1:3),par10,pos10(:,3),camera);
    Xerr10=[Xi101-data1(:,4);Xi102-data3(:,4);Xi103-data4(:,4)];
    Yerr10=[Yi101-data1(:,5);Yi102-data3(:,5);Yi103-data4(:,5)];
    Rerr10=sqrt(mean(Xerr10.^2 + Yerr10.^2));
  end
  
  [par11,pos11,iter,res,er,C,success]=cacal(camera,data1,data2,data3,data4);
  Rerr11=NaN;
  if success
    [Xi111,Yi111]=pred(data1(:,1:3),par11,pos11(:,1),camera);
    [Xi112,Yi112]=pred(data2(:,1:3),par11,pos11(:,2),camera);
    [Xi113,Yi113]=pred(data3(:,1:3),par11,pos11(:,3),camera);
    [Xi114,Yi114]=pred(data4(:,1:3),par11,pos11(:,4),camera);
    Xerr11=[Xi111-data1(:,4);Xi112-data2(:,4);Xi113-data3(:,4);Xi114-data4(:,4)];
    Yerr11=[Yi111-data1(:,5);Yi112-data2(:,5);Yi113-data3(:,5);Yi114-data4(:,5)];
    Rerr11=sqrt(mean(Xerr11.^2 + Yerr11.^2));
  end
  
  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 data2 data3 data4 camera',calfile); 
  
  Rerr=[Rerr1,Rerr2,Rerr3,Rerr4,Rerr5,Rerr6,Rerr7,Rerr8,Rerr9,Rerr10,Rerr11];
  minRerr=find(Rerr==min(Rerr));
  sprintf('Min Rerr achieved with case %d with Rerr=%5.2f',minRerr,min(Rerr))
  eval(sprintf('pos=pos%d; par=par%d;',minRerr,minRerr));
  
%elseif imgnum==5 
%  [par,pos,iter,res,er,C,success]=cacal(camera,data1,data2,data3,data4,data5);
%  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 data2 data3 data4 data5',calfile); 
%elseif imgnum==6 
%  [par,pos,iter,res,er,C,success]=cacal(camera,data1,data2,data3,data4,data5,data6);
%  savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data1 data2 data3 data4 data5 data6',calfile); 
end

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

return
%-----------------------------------------
%Old stuff down here

Mcr=[1,     0,      0, 0;
     0,     1,      0, 0;
     0,     0,      1, 0
 43.21,-38.20,-272.81, 1];                         % From cube basis to reference basis, translate  
 
 
 
Trc=[-43.21,38.20,272.81];  % origin of cube basis in reference frame
%Trw=[-43,38,33.75];         % origin of world basis in reference frame
Trw=Trc;
Twc=Trc-Trw;
theta=0.606  % radians
Mcr=[cos(theta), sin(theta),       0, 0;    
  -1*sin(theta), cos(theta),       0, 0;
         0     ,     0     ,       1, 0;
	 0     ,     0     ,       0, 1];   % Align cube basis with reference basis, rotate by theta
Mcw=[1,      0,      0, 0;
     0,      1,      0, 0;
     0,      0,      1, 0
Twc(1), Twc(2), Twc(3), 1];                 % Transform cube basis to world basis, vector addition
    
xyz1c=xyz1c*Mcr; % First rotate cube basis to align with reference basis (ie: xc // xr, yc // yr)
xyz1w=xyz1c*Mcw; % Represent cube coords in terms of world coords

vec1c=vec1c*Mcr;     % First rotate cube basis to align with reference basis (ie: xc // xr, yc // yr)         
vec1w=vec1c*Mcw;     % Represent cube coords in terms of world coords          
vec1Nw=vec1w-xyz1w;  % Re-compute normals, no need to renormalize since nondim vector decomposition assumed in mm

%xi2=xi+1.*randn(1,length(xi));
%yi2=yi+1.*randn(1,length(xi));

