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


for n = 1:length(angles)
    % load calibration data
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

% Calculation of camera location
xyz = zeros(4,length(angles));
xyzNorm = zeros(4,length(angles));
Rc = zeros(4,4,length(angles));
for j = 1:length(angles)
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
    % Add translation portion of transformation matrix
    Rc(:,:,j)=[[R'; [0 0 0]] [POS(1,j) POS(2,j) POS(3,j) 1]'];  % store each R in a 3d Rc matrix for easy referencing later
    % Use transformation matrix on origin to find camera location use of \
    % indicates multiplication by the inverse of Rc
    xyz(:,j) = Rc(:,:,j)\[0 0 0 1]';
    xyzNorm(:,j) = Rc(:,:,j)\[0 0 300 1]';
end


ptColor = 'mckgb';
ptMarker = 'oooo+';

% % % figure
% % % scatter3(0,0,0,'rx','SizeData',96)
% % % hold on
for n = 1:length(angles)
    figure
    scatter3(0,0,0,'rx','SizeData',96)
    axis equal
    hold on
    if rem(n,2) == 0
        % Calibration points
        tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''MarkerEdgeColor'',''none'',''MarkerFaceColor'',''%s'')',ptColor(n),ptMarker(n),ptColor(n));
        eval(tmp)
    else
        tmp = sprintf('scatter3(DATA{n}(:,1),DATA{n}(:,2),DATA{n}(:,3),''%s%s'',''SizeData'',48,''LineWidth'',2)',ptColor(n),ptMarker(n));
        eval(tmp)
    end
   % Camera Positions
   tmp = sprintf('scatter3(xyz(1,n),xyz(2,n),xyz(3,n),''%s%s'',''SizeData'',48)',ptColor(n),ptMarker(n));
   eval(tmp)
   % Camera Orientation
   tmp = sprintf('plot3([xyz(1,n) xyzNorm(1,n)],[xyz(2,n) xyzNorm(2,n)],[xyz(3,n) xyzNorm(3,n)],''%s'',''LineWidth'',2)',ptColor(n));
   eval(tmp)
end
axis equal