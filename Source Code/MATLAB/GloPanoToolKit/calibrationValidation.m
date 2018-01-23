%% Calibration Validation %%

angles = [45 315 225 135 45];
labels{5} = 'Geo';
labels{1} = 'A';
labels{2} = 'B';
labels{3} = 'C';
labels{4} = 'D';

A = cell(length(angles),1);
ANG = zeros(length(angles),1);
c = zeros(14,14,length(angles));
CAMERA = cell(length(angles),1);
CWTXT = cell(length(angles),1);
DATA = cell(length(angles),1);
ER = cell(length(angles),1);
ITER = zeros(length(angles),1);
PAR = zeros(8,length(angles));
POS = zeros(6,length(angles));
RES = zeros(length(angles),1);

dirContents = dir;
check = zeros(size(dir,1),5);

% % % for n = 1:length(angles)
for n = 1:5
    % Calibration fileame
    tmp = sprintf('cal%s_%0.3d.mat',labels{n},angles(n));
    % load calibration data
    for m = 1:size(dirContents,1)
       if strcmp(getfield(dirContents,{m},'name'),tmp)
        check(m,n) =  1;
       end
    end
    
    if sum(check(:,n)) > 0
        tmp = sprintf('load(''cal%s_%0.3d.mat'')',labels{n},angles(n));
        eval(tmp)
        A{n} = a;
        ANG(n) = ang;
        c(:,:,n) = C;
        CAMERA{n} = camera;
        DATA{n} = data;
        ER{n} = er;
        ITER(n) = iter;
        PAR(:,n) = par;
        POS(:,n) = pos;
        RES(n) = res;
        clear a ang C camera cwtxt data er iter par pos res
    end
end

check = sum(check).*(1:size(check,2));
if sum(check == 0) ~= 0
   check = unique(check);
   check = check(2:end);
end

% % % axesCheck = [0 10 0 0; 0 0 10 0; 0 0 0 10;];
% % % 
% % % [Xi,Yi,pc] = predMod(axesCheck',PAR(:,1),POS(:,1),CAMERA{1});
% % % 
% % % ptColor = 'bcmgk';
% % % ptMarker = 'oooo+';
% % % 
% % % figure
% % % hold on
% % % if rem(n,2) == 0
% % %     % Calibration points
% % %     tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''MarkerEdgeColor'',''none'',''MarkerFaceColor'',''%s'')',ptColor(n),ptMarker(n),ptColor(n));
% % %     eval(tmp)
% % % else
% % %     tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''LineWidth'',2)',ptColor(n),ptMarker(n));
% % %     eval(tmp)
% % % end
% % % for n = 2:4    
% % %     plot3([axesCheck(1,1) axesCheck(1,n)],[axesCheck(2,1) axesCheck(2,n)],[axesCheck(3,1) axesCheck(3,n)],'b','LineWidth',2)
% % %     plot3([pc(1,1) pc(n,1)],[pc(1,2) pc(n,2)],[pc(1,3) pc(n,3)],'r','LineWidth',2)
% % % end

% Calculation of camera location
xyz = zeros(4,length(angles));
Rc = zeros(4,4,length(angles));
for j = check
    % Convert values from degrees in to radians
    wa=POS(4,j)*pi/180;  % omega, Wx rotation (world basis, x axis)
    pa=POS(5,j)*pi/180;  % psi, y rotation  (world basis, y axis)
    ra=POS(6,j)*pi/180;  % kappa, z rotation (world basis, z axis)
    cw=cos(wa); sw=sin(wa);
    cp=cos(pa); sp=sin(pa);
    cr=cos(ra); sr=sin(ra);
    R=zeros(3,3);
    % Construct rotation portion of tranformation matrix
    R(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]';             % build the 4x4 R matrix
    R(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
    R(:,3)=[-sp cp*sw cp*cw]';
% % % R(:,1) = [cp*cr,-cp*sr,sp];
% % % R(:,2) = [cr*sw*sp+cw*sr, cw*cr-sw*sp*sr, -cp*sw];
% % % R(:,3) = [sw*sr-cw*cr*sp, cr*sw+cw*sp*sr, cw*cp];
    % Add translation portion of transformation matrix
% %     Rc(:,:,j)=[[R'; [0 0 0]] [POS(1,j) POS(2,j) POS(3,j) 1]'];  % store each R in a 3d Rc matrix for easy referencing later
Rc(:,:,j)=[[R'; [0 0 0]] [POS(1,j) POS(2,j) POS(3,j) 1]'];  % store each R in a 3d Rc matrix for easy referencing later
    % Use transformation matrix on origin to find camera location use of \
    % indicates multiplication by the inverse of Rc
    xyz(:,j) = Rc(:,:,j)\[0 0 0 1]';

end

% Grab camera view directions
lookmap=ones(size(Rc,3),4,3);
looknmap=ones(size(Rc,3),3,3);
for i=1:size(Rc,3)
    lookmap(i,:,3)=(Rc(:,:,i)\[0 0 10 1]')';       % camera view direction
    lookmap(i,:,1)=(Rc(:,:,i)\[10 0 0 1]')';       % x-axis 
    lookmap(i,:,2)=(Rc(:,:,i)\[0 10 0 1]')';       % y-axis
    lookmap(i,1:3,3)=lookmap(i,1:3,3)-xyz(1:3,i)';     % camera view direction
    lookmap(i,1:3,1)=lookmap(i,1:3,1)-xyz(1:3,i)';
    lookmap(i,1:3,2)=lookmap(i,1:3,2)-xyz(1:3,i)';
    looknmap(i,:,3)=lookmap(i,1:3,3)./norm(lookmap(i,1:3,3));     % normalized view direction
end


ptColor = 'bcmgk';
ptMarker = 'oooo+';

figure
% scatter3(0,0,0,'rx','SizeData',96)
hold on
plot3([0 10],[0 0],[0 0],'k','LineWidth',2)
plot3([0 0],[0 10],[0 0],'k','LineWidth',2)
plot3([0 0],[0 0],[0 10],'k','LineWidth',2)
% for n = 1:length(angles)
for n = 1:5
    if rem(n,2) == 0
        % Calibration points
        tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''MarkerEdgeColor'',''none'',''MarkerFaceColor'',''%s'')',ptColor(n),ptMarker(n),ptColor(n));
        eval(tmp)
    else
        tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''LineWidth'',2)',ptColor(n),ptMarker(n));
        eval(tmp)
    end
   % Camera Positions
%    tmp = sprintf('scatter3(xyz(1,n),xyz(2,n),xyz(3,n),''%s%s'',''SizeData'',48)',ptColor(n),ptMarker(n));
   eval(tmp)
   % Camera Orientation
   tmp = sprintf('plot3([xyz(1,n) xyz(1,n)+lookmap(n,1,1)],[xyz(2,n) xyz(2,n)+lookmap(n,2,1)],[xyz(3,n) xyz(3,n)+lookmap(n,3,1)],''%s'',''LineWidth'',2)','r');
   eval(tmp)
   tmp = sprintf('plot3([xyz(1,n) xyz(1,n)+lookmap(n,1,2)],[xyz(2,n) xyz(2,n)+lookmap(n,2,2)],[xyz(3,n) xyz(3,n)+lookmap(n,3,2)],''%s'',''LineWidth'',2)','r');
   eval(tmp)
   tmp = sprintf('plot3([xyz(1,n) xyz(1,n)+lookmap(n,1,3)],[xyz(2,n) xyz(2,n)+lookmap(n,2,3)],[xyz(3,n) xyz(3,n)+lookmap(n,3,3)],''%s'',''LineWidth'',2)','r');
   eval(tmp)
end
axis equal

