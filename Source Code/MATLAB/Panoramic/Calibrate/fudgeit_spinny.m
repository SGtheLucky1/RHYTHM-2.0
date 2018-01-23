clear all
close all

% Calibration points defined in file.
% 3D positions in reference to a fixed right-handed frame!
%

%----------------------
% User defined values here
zthresh=200;
zrange=60;
method=2;            % 1 for click center, 2 for line intersect
%----------------------

expname=input('Experiment name [ie, map05_07_15_03]: ','s');

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

cw=0;
while (cw~=1 && cw~=2)
  cw=input('Calibrated cw or ccw? [1:cw, 2:ccw]: ');
end
if cw==1, cwtxt='cw'; end;
if cw==2, cwtxt='ccw'; end;
% ang=input('Calibration angle (ie: 0, 45, 90, 135, 180, 225, 270, or 315): ');
angle = [45 135 225 315];
fprintf('Calibrating for angles %d, %d, %d, and %d.\n',angle(1),angle(2),angle(3),angle(4)); 

plane1 = cell(length(angle),1);
plane1fname = cell(length(angle),1);
plane2 = cell(length(angle),1);
plane2fname = cell(length(angle),1);
calpts1 = cell(length(angle),1);
calpts2 = cell(length(angle),1);
calpts_nd = cell(length(angle),1);
for n = 1:length(angle)
    fprintf('%d degree angle:\n',angle(n))
    plane1{n}=input('Enter first plane (right hand) to calibrate (e.g. posX, negX, posY, negY) for angle: ','s');
    plane1fname{n}=sprintf('calpts2_%s.txt',plane1{n});
    fid=fopen(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname{n}));
    if fid~=-1
        fclose(fid);
        fprintf('Loading %s ...\n',plane1fname{n});
        calpts1{n}=load(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane1fname{n}));
    else
        fprintf('Could not load %s!',plane1fname{n});
        return
    end
    plane2{n}=input('Enter second plane (the one on the left) to calibrate (e.g. posX, negX, posY, negY): ','s');
    plane2fname{n}=sprintf('calpts2_%s.txt',plane2{n});
    fid=fopen(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname{n}));
    if fid~=-1
        fclose(fid);
        fprintf('Loading %s ...\n',plane2fname{n});
        calpts2{n}=load(sprintf('/Users/Chris/Documents/MATLAB/Panoramic/Calibrate/calpts/%s/%s',species,plane2fname{n}));
    else
        fprintf('Could not load %s\n!',plane2fname{n});
        return
    end
    calpts_nd{n}=[calpts1{n};calpts2{n}];       % nondimensional calibration point positions in cube basis
end

%identify image positions that correspond to chosen angles
anum = cell(length(angle),1);
for n = 1:length(angle)
    if angle(n)==0
        anum{n}='000';
    elseif angle(n)==45
        anum{n}='009';
    elseif angle(n)==90
        anum{n}='018';
    elseif angle(n)==135
        anum{n}='027';
    elseif angle(n)==180
        anum{n}='036';
    elseif angle(n)==225
        anum{n}='045';
    elseif angle(n)==270
        anum{n}='054';
    elseif angle(n)==315
        anum{n}='063';
    end
end

%identify calibration directory and calibration image
done=0;
while done == 0
  defaultdir=sprintf('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/%s/Geometry/cal/',expname);
  defaultfile=sprintf('cal_%s.pgm',anum{1});
  eval(sprintf('caldir=input(''Calibration directory: [%s] '',''s'');',defaultdir));
  if isempty(caldir)
      caldir=defaultdir; 
  end
  eval(sprintf('calfile=input(''Calibration image filename: [%s] '',''s'');',defaultfile));
  if isempty(calfile)
      calfile=defaultfile;
  end
  calfilename=sprintf('%s%s',caldir,calfile);
  fid=fopen(calfilename);
  if fid~=-1
    fclose(fid);
    done=1;
  end
end
fprintf('Loading %s',calfilename);
a=imread(calfilename);
a=a(:,:,1);
%a=imread(calfilename);
% % if size(a,3)>1
% %     a=rgb2gray(a);
% % end
data = [];
for n = 2:length(angle)
% % for n = 1
    fid=fopen(sprintf('%s_%s.mat',calfile(1:end-5),anum{n}));
    if fid==-1 && isempty(data)
        data=crosseyed(a,calpts_nd{n},calfile,method,zrange,zthresh);  % THE MEAT OF THE FRUIT
        tmp = sprintf('data%s=data;',anum{n});
        eval(tmp)
    elseif fid==-1 && ~isempty(data)
        tmp = sprintf('data%s=[calpts_nd{n}(pt_id,4:6) data(:,4:5) calpts_nd{n}(pt_id,1:3)];',anum{n});
        eval(tmp)
    else
        fclose(fid);
        fprintf('%s.mat exists!',calfile);
        useit='j';
        while (useit~='Y' && useit~='y' && useit~='N' && useit~='n')
            useit=input(sprintf('Want to use the pixel locations from %s.mat? [y/n]: ',calfile),'s');
        end
        if (useit=='N' || useit=='n')
            data=crosseyed(a,calpts_nd{n},calfile,method,zrange,zthresh);  % THE MEAT OF THE FRUIT
        else
            fprintf('Loading data matrix from %s.mat ...',calfile);
            eval(sprintf('load %s.mat data',calfile));
            disp('Extracting pixel locations and creating new data matrix ...');
            data=[calpts_nd{n}(:,4:6), data(:,4), data(:,5), calpts_nd{n}(:,1:3)];
        end
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
for n = 1:length(angle)
    if calselect==1
        disp('Using 3/8 inch spacing (pig with BrainVision).');
        tmp = sprintf('data%s(:,1:3)=data%s(:,1:3).*(25.4*3/8);',anum{n},anum{n});   
        eval(tmp)   % convert positions to mm
    elseif calselect==2
        disp('Using 1/4 inch spacing (rabbit)');
        tmp = sprintf('data%s(:,1:3)=data%s(:,1:3).*(25.4*1/4);',anum{n},anum{n});
        eval(tmp)   % convert positions to mm
    elseif calselect==3
        disp('Using 0.6 inch spacing (pig with Andor)');
        tmp = sprintf('data%s(:,1:3)=data%s(:,1:3).*(25.4*0.6)',anum{n},anum{n});
        eval(tmp)   % convert positions to mm
    end
end

% % if calselect==1
% %   disp('Using 3/8 inch spacing (pig with BrainVision).');
% %   data(:,1:3)=data(:,1:3).*(25.4*3/8);       % convert positions to mm
% % elseif calselect==2
% %   disp('Using 1/4 inch spacing (rabbit)');
% %   data(:,1:3)=data(:,1:3).*(25.4*1/4);       % convert positions to mm
% % elseif calselect==3
% %   disp('Using 0.6 inch spacing (pig with Andor)');
% %   data(:,1:3)=data(:,1:3).*(25.4*0.6);       % convert positions to mm
% % end

%-----------------------------------------
% Calibrate

fprintf('\n Enter number for camera: \n');
fprintf('    1) watec_with_f8.5 (pig) \n');
fprintf('    2) watec_with_f12.5 (rabbit) \n');
fprintf('    3) Nikon E4300 (rabbit fibers) \n');
fprintf('    4) interpolated_dalsa \n');
fprintf('    5) iDS_UI_3220CP-M-GL_with_f1.2 \n');
done=0;
while done == 0
    camnum=input('Camera number [5]: ');
    if isempty(camnum)
        camnum=5;
    end
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
        case 5
            camera='iDS_UI_3220CP-M-GL_with_f1.2';
            done=1;
    end
end

for n = 1:length(angle)
    calibcommand=sprintf('[par,pos,iter,res,er,C,success]=cacal(camera,data%s);',anum{n});
    eval(calibcommand)
%     [par,pos,iter,res,er,C,success]=cacal(camera,data);
    calpts_nondim = calpts_nd{n};
    ang = angle(n);
    datacommand = sprintf('data=data%s;',anum{n});
    eval(datacommand)
    savecommand=sprintf('save calOld_%s.mat a calpts_nondim par pos iter res er C data camera ang cwtxt',anum{n});
    eval(savecommand)
    if n == 1
        predcommand = sprintf('[Xi,Yi]=pred(data%s(:,1:3),par,pos,camera);',anum{n});
        eval(predcommand)
    end
    
    % Save calibration parameters
    posparfname=sprintf('calOld_%s.pospar',anum{n});
    posparfid=fopen(posparfname,'w');
    for i=1:6
        fprintf(posparfid,'%15e\n',pos(i));
    end
    for i=1:8
        fprintf(posparfid,'%15e\n',par(i));
    end
    fclose(posparfid);
end

close all
imshow(a)
imagesc(a); colormap('gray');
hold on
xi=data(:,4);
yi=data(:,5);
plot(xi+0.5,yi+0.5,'go','MarkerFaceColor','g');         % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
plot(Xi+0.5,Yi+0.5,'ro','MarkerFaceColor','r');         % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
title('Predicted calibration points','FontSize',14);
jpegcommand=sprintf('print -djpeg %s.jpg',calfile);
eval(jpegcommand);

