clear all
close all

% Calibration points defined in file.
% 3D positions in reference to a fixed right-handed frame!
%

%----------------------
% User defined values here
zthresh=180;
zrange=300;
method=2;            % 1 for click center, 2 for line intersect
%----------------------

expname=input('Experiment name [ie, map05_07_15_03]: ','s');

speciesID=0;
disp('Animal Species:');
disp('  1: Pig');
disp('  2: Rabbit');
while (speciesID~=1 & speciesID~=2)
  speciesID=input('Enter animal species [1 or 2]: ');
end
if speciesID==1
  disp('Using PIG calibration data.');
  species='pig';
elseif speciesID==2
  disp('Using RABBIT calibration data.');
  species='rabbit';
end

disp(sprintf(' posX \n negX \n posY \n negY'));
plane1=input('Enter first plane to calibrate (the one on the right): ','s');
plane1fname=sprintf('calpts2_%s.txt',plane1);
fid=fopen(sprintf('/local/home/mwk/matlab/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname));
if fid~=-1
  fclose(fid);
  disp(sprintf('Loading %s ...',plane1fname));
  calpts1=load(sprintf('/local/home/mwk/matlab/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname));
else
  disp(sprintf('Could not load %s!',plane1fname));
  return
end
plane2=input('Enter second plane to calibrate (the one on the left): ','s');
plane2fname=sprintf('calpts2_%s.txt',plane2);
fid=fopen(sprintf('/local/home/mwk/matlab/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname));
if fid~=-1
  fclose(fid);
  disp(sprintf('Loading %s ...',plane2fname));
  calpts2=load(sprintf('/local/home/mwk/matlab/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname));
else
  disp(sprintf('Could not load %s!',plane2fname));
  return
end
calpts_nondim=[calpts1;calpts2];       % nondimensional calibration point positions in cube basis

done=0;
while ~done
  defaultdir=sprintf('/disk2/raw_data/%s/fiber_pics_cal/',expname);
  defaultfile=sprintf('fiber_cal.jpg');
  eval(sprintf('caldir=input(''Calibration directory: [%s] '',''s'');',defaultdir));
  if isempty(caldir), caldir=defaultdir; end;
  eval(sprintf('calfile=input(''Calibration image filename: [%s] '',''s'');',defaultfile));
  if isempty(calfile), calfile=defaultfile; end;
  calfilename=sprintf('%s%s',caldir,calfile);
  fid=fopen(calfilename);
  if fid~=-1
    fclose(fid);
    done=1;
  end
end

disp(sprintf('Loading %s',calfilename));
a=imread(calfilename);
%a=imread(calfilename);
if size(a,3)>1, a=rgb2gray(a); end;

fid=fopen(sprintf('%s.mat',calfile));
if fid==-1
  data=crosseyed(a,calpts_nondim,calfile,method,zrange,zthresh,0);  % THE MEAT OF THE FRUIT
else
  fclose(fid);
  disp(sprintf('%s.mat exists!',calfile));
  useit='j';
  while (useit~='Y' & useit~='y' & useit~='N' & useit~='n')
    useit=input(sprintf('Want to use the pixel locations from %s.mat? [y/n]: ',calfile),'s');
  end
  if (useit=='N' | useit=='n')
    data=crosseyed(a,calpts_nondim,calfile,method,zrange,zthresh,0);  % THE MEAT OF THE FRUIT
  else
    disp(sprintf('Loading data matrix from %s.mat ...',calfile));
    eval(sprintf('load %s.mat data',calfile));
    disp('Creating extracting pixel locations and creating new data matrix ...');
    data=[calpts_nondim(:,4:6), data(:,4), data(:,5), calpts_nondim(:,1:3)];
  end
end

disp('Enter calibration target grid spacing (inches):');
disp('   1: 3/8 inch for the pig calibration target');
disp('   2: 1/4 inch for the rabbit calibration target');
calselect=0;
while (calselect~=1 & calselect~=2)
  calselect=input('Enter grid spacing [1 or 2]: ');
end
if calselect==1
  disp('Using 3/8 inch spacing (pig).');
  data(:,1:3)=data(:,1:3).*(25.4*3/8);       % convert positions to mm
elseif calselect==2
  disp('Using 1/4 inch spacing (rabbit)');
  data(:,1:3)=data(:,1:3).*(25.4*1/4);       % convert positions to mm
end

%-----------------------------------------
% Calibrate

disp(sprintf('\n Enter number for camera: \n'));
disp(sprintf('    1) watec_with_f8.5 (pig) \n'));
disp(sprintf('    2) watec_with_f12.5 (rabbit) \n'));
disp(sprintf('    3) Nikon E4300 (rabbit fibers) \n'));
disp(sprintf('    4) interpolated_dalsa \n'));
done=0;
while ~done
  camnum=input('Camera number [1]: ');
  if isempty(camnum), camnum=1; end;
  switch camnum
    case 1
      camera='watec_with_f8.5';
      done=1;
    case 2
      camera='watec_with_f12.5';
      done=1;
    case 3
      camera='E4300';
      done=1;   
    case 4
      camera='interpolated_dalsa';
      done=1;
  end
end

[par,pos,iter,res,er,C,success]=cacal(camera,data);
savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data camera',calfile); 
[Xi,Yi]=pred(data(:,1:3),par,pos,camera); 
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
imagesc(a); colormap('gray');
hold on
xi=data(:,4);
yi=data(:,5);
plot(xi+0.5,yi+0.5,'go','MarkerFaceColor','g');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
plot(Xi+0.5,Yi+0.5,'ro','MarkerFaceColor','r');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
title('Predicted calibration points','FontSize',14);
jpegcommand=sprintf('print -djpeg %s.jpg',calfile);
eval(jpegcommand);

