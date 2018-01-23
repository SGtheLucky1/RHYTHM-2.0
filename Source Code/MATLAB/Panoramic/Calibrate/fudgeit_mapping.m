clear all
close all

% Calibration points defined in file.
% 3D positions in reference to a fixed right-handed frame!
%

%----------------------
% User defined values here
zthresh=90;
zrange=30;
method=2;            % 1 for click center, 2 for line intersect
%----------------------

speciesID=0;
disp('Animal Species:');
disp('  1: Pig (BrainVision Cameras)');
disp('  2: Rabbit');
disp('  3: Pig (Andor Cameras)');
while (speciesID~=1 && speciesID~=2 && speciesID~=3)
  speciesID=input('Enter type of study [1, 2, or 3]: ');
end

if speciesID==1
  disp('Using PIG calibration data for BrainVision cameras.');
  species='pig_brainvision';
elseif speciesID==2
  disp('Using RABBIT calibration data.');
  species='rabbit';
elseif speciesID==3
  disp('Using PIG calibration data for Andor cameras.');
  species='pig_andor';
end

disp(sprintf(' posX \n negX \n posY \n negY'));
plane1=input('Enter first plane to calibrate (the one on the right): ','s');
plane1fname=sprintf('calpts2_%s.txt',plane1);
fid=fopen(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname));
if fid~=-1
  fclose(fid);
  disp(sprintf('Loading %s ...',plane1fname));
  calpts1=load(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname));
else
  disp(sprintf('Could not load %s!',plane1fname));
  return
end
plane2=input('Enter second plane to calibrate (the one on the left): ','s');
plane2fname=sprintf('calpts2_%s.txt',plane2);
fid=fopen(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname));
if fid~=-1
  fclose(fid);
  disp(sprintf('Loading %s ...',plane2fname));
  calpts2=load(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname));
else
  disp(sprintf('Could not load %s!',plane2fname));
  return
end
calpts_nondim=[calpts1;calpts2];       % nondimensional calibration point positions in cube basis

calfile=input('Enter calibration image filename: ','s');
fid=fopen(calfile);
if fid==-1
  disp(sprintf('Could not open %s!',calfile));
  return
end
fclose(fid);
disp(sprintf('Loading %s ... ',calfile));
a=imread(calfile);
if size(a,3)>1, a=rgb2gray(a); end;

fid=fopen(sprintf('%s.mat',calfile));
if fid==-1
  data=crosseyed(a,calpts_nondim,calfile,method,zrange,zthresh);  % THE MEAT OF THE FRUIT
else
  fclose(fid);
  disp(sprintf('%s.mat exists!',calfile));
  useit='j';
  while (useit~='Y' && useit~='y' && useit~='N' && useit~='n')
    useit=input(sprintf('Want to use the pixel locations from %s.mat? [y/n]: ',calfile),'s');
  end
  if (useit=='N' || useit=='n')
    data=crosseyed(a,calpts_nondim,calfile,method,zrange,zthresh);  % THE MEAT OF THE FRUIT
  else
    disp(sprintf('Loading data matrix from %s.mat ...',calfile));
    eval(sprintf('load %s.mat data',calfile));
    disp('Extracting pixel locations and creating new data matrix ...');
    data=[calpts_nondim(:,4:6), data(:,4), data(:,5), calpts_nondim(:,1:3)];
  end
end

disp('Enter calibration target grid spacing (inches):');
disp('   1: 3/8 inch for the pig with BrainVision calibration target');
disp('   2: 1/4 inch for the rabbit calibration target');
disp('   3: 0.6 inch for the pig with Andor calibration target');
calselect=0;
while (calselect~=1 && calselect~=2 && calselect~=3)
  calselect=input('Enter grid spacing [1, 2, or 3]: ');
end
if calselect==1
  disp('Using 3/8 inch spacing (pig with BrainVision).');
  data(:,1:3)=data(:,1:3).*(25.4*3/8);       % convert positions to mm
elseif calselect==2
  disp('Using 1/4 inch spacing (rabbit)');
  data(:,1:3)=data(:,1:3).*(25.4*1/4);       % convert positions to mm
elseif calselect==3
  disp('Using 0.6 inch spacing (pig with Andor)');
  data(:,1:3)=data(:,1:3).*(25.4*0.6);       % convert positions to mm
end

%-----------------------------------------
% Calibrate

disp(sprintf('\n Enter number for camera: \n'));
disp(sprintf('    1) watec_and_snappy \n'));
disp(sprintf('    2) interpolated_dalsa \n'));
disp(sprintf('    3) brainvision_high \n'));
disp(sprintf('    4) brainvision_low \n'));
disp(sprintf('    5) andor_ixon_6mm \n'));
disp(sprintf('    6) andor_ixon_2.2mm \n'));
done=0;
while ~done
  camnum=input('Camera number [5]: ');
  if isempty(camnum), camnum=5; end;
  switch camnum
    case 1
      camera='watec_and_snappy';
      done=1;
    case 2
      camera='interpolated_dalsa';
      done=1;
    case 3
      camera='brainvision_high';
      done=1;
    case 4
      camera='brainvision_low';
      done=1;  
    case 5
      camera='andor_ixon_6mm';
      PIXfname=input('Enter filename with pixel offsets [camA.pix]: ','s');
      if isempty(PIXfname); PIXfname='camA.pix'; end;
      fid=fopen(PIXfname,'r');
      if fid==-1
        disp(sprintf('Could not open %s!',PIXfname));
      else
        fclose(fid);
        disp(sprintf('Loading pixel offsets from %s ... ',PIXfname));
	p=load(PIXfname);
	hoffset=p(1)-1;
	voffset=p(4)-1;
	disp(sprintf('applying vertical pixel offset: %d\napplying horizontal pixel offset: %d\n',voffset,hoffset));
	data(:,4)=data(:,4)+hoffset;
	data(:,5)=data(:,5)+voffset;
	done=1;  
      end
    case 6
      camera='andor_ixon_2.2mm';
      PIXfname=input('Enter filename with pixel offsets [camA.pix]: ','s');
      if isempty(PIXfname); PIXfname='camA.pix'; end;
      fid=fopen(PIXfname,'r');
      if fid==-1
        disp(sprintf('Could not open %s!',PIXfname));
      else
        fclose(fid);
        disp(sprintf('Loading pixel offsets from %s ... ',PIXfname));
	p=load(PIXfname);
	hoffset=p(1)-1;
	voffset=p(4)-1;
	disp(sprintf('applying vertical pixel offset: %d\napplying horizontal pixel offset: %d\n',voffset,hoffset));
	data(:,4)=data(:,4)+hoffset;
	data(:,5)=data(:,5)+voffset;
	done=1;  
      end 
  end
end

[par,pos,iter,res,er,C,success]=cacal(camera,data);
savecommand=sprintf('save %s.mat a calpts_nondim par pos iter res er C data camera',calfile); 
[Xi,Yi]=pred(data(:,1:3),par,pos,camera); 
eval(savecommand);

if (camnum==5 || camnum==6)
  Xi=Xi-hoffset;
  Yi=Yi-voffset;
end

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
if (camnum==5 || camnum==6)
  xi=data(:,4)-hoffset;
  yi=data(:,5)-voffset;
else
  xi=data(:,4);
  yi=data(:,5);
end
plot(xi+0.5,yi+0.5,'go','MarkerFaceColor','g');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
plot(Xi+0.5,Yi+0.5,'ro','MarkerFaceColor','r');    % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
axis('equal');
title('Predicted calibration points','FontSize',14);
jpegcommand=sprintf('print -djpeg %s.jpg',calfile);
eval(jpegcommand);

done=0;
while ~done
  yn=input('Calibrate for low resolution mode? [y/n]: ','s');
  if ~isempty(yn) && (yn=='y' || yn=='Y' || yn=='n' || yn=='N')
    done=1;
  end
end

if (yn=='y' || yn=='Y')
  real_xshift=3.5;
  mirror_xshift=-2.5;
  xshift=input(sprintf('Enter horizontal shift (%4.3f pixels for real images, %4.3f pixels for mirror images): ',real_xshift,mirror_xshift));
  disp(sprintf('Shifting by %4.3f pixels.',xshift));
  data_low=data;
  data_low(:,4)=data_low(:,4)./2 + xshift;
  data_low(:,5)=data_low(:,5)./2; 
  camera='brainvision_low';
  disp(sprintf('camera: %s',camera));
  [par_low,pos_low,iter,res,er,C,success]=cacal(camera,data_low);
  [xi,yi]=pred(data_low(:,1:3),par_low,pos_low,camera);
 
  % Save calibration parameters
  posparfname=sprintf('%s_lowres.pospar',calfile);
  posparfid=fopen(posparfname,'w');
  for i=1:6
    fprintf(posparfid,'%15e\n',pos_low(i));
  end
  for i=1:8
    fprintf(posparfid,'%15e\n',par_low(i));
  end
  fclose(posparfid);

end
