close all
clear all
disp(sprintf('\npan_coloranatomy.m \nVersion: 03/08/05 \n\n'));
disp('Press Return ...'); pause;

setup_pan_coloranatomy

%------------------------------------
% Load spinny camera calibration stuff

[Rcspin,Rc_ntheta,Rc_dtheta,message]=readRccal(char(spincalfnames(1)));
[Parspin,Par_ntheta,Par_dtheta,message]=readParcal(char(spincalfnames(2)));

if Rc_dtheta==Par_dtheta
  dtheta=Rc_dtheta;
  ntheta=Rc_ntheta;
  disp(sprintf('%s, dtheta: %d, ntheta: %d',char(spincalfnames(1)),Rc_dtheta,Rc_ntheta));
  disp(sprintf('%s, dtheta: %d, ntheta: %d',char(spincalfnames(2)),Par_dtheta,Par_ntheta));
  dtheta=Rc_dtheta;
else
  disp(sprintf('dthetha mismatch in %s (%d) and %s (%d)! aborting...',char(spincalfnames(1)),Rc_dtheta,char(spincalfnames(2)),Par_dtheta));
  return
end

disp(sprintf('\nCamera parameters loaded.\n'));
disp(sprintf('\nWas the object scanned in the same direction as the camera'));
disp(sprintf('was calibrated? ie: camera calibrated cw and object scanned ccw? or'));
disp(sprintf('                    camera calibrated ccw and object scanned cw?'));
disp(sprintf('If so, then camera calibration parameters need to be re-ordered.'));
if isempty(reorder), reorder='g'; end
while (isempty(reorder) | (reorder~='y' & reorder~='Y' & reorder~='n' & reorder~='N'))
  reorder=input('Should camera calibration parameters be re-ordered? (y/n): ','s');
end
if (reorder=='y' | reorder=='Y')
  disp(sprintf('\n Reordering camera parameters .... \n'));
  Parspin(1:8,2:end)=fliplr(Parspin(1:8,2:end));
  Rcspin_new=Rcspin;
  for i=2:size(Rcspin,3)
    Rcspin_new(1:4,1:4,i)=Rcspin(1:4,1:4,size(Rcspin,3)-i+2);
  end
  Rcspin=Rcspin_new;
  clear Rcspin_new
else
  disp(sprintf('\n Ok, will not reorder camera parameters .... \n'));
end

angles=[0:dtheta:360-dtheta];
anglesi=imnums+1;               % imnums are image numbers defined by user
nviews=length(imnums);
image_angles=angles(anglesi);

% check spinny camera calibration stuff
xyzspin=ones(size(Rcspin,3),4);
lookspin=ones(size(Rcspin,3),4);
looknspin=ones(size(Rcspin,3),3);
for i=1:size(Rcspin,3)
  xyzspin(i,:)=(inv(Rcspin(1:4,1:4,i))*[0 0 0 1]')';         % camera origin in world coords
  lookspin(i,:)=(inv(Rcspin(1:4,1:4,i))*[0 0 1 1]')';        % camera view direction
  lookspin(i,1:3)=lookspin(i,1:3)-xyzspin(i,1:3);            % camera view direction
  looknspin(i,:)=lookspin(i,1:3)./norm(lookspin(i,1:3));     % normalized view direction
end

if 0
tt=figure; hold on; axis('square');
axis([min(xyzspin(:,1)) max(xyzspin(:,1)) min(xyzspin(:,2)) max(xyzspin(:,2)) min(xyzspin(:,3)) max(xyzspin(:,3)) ]); 
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
for i=1:size(Rcspin,3)
  disp(sprintf('Angle %d',Rcspin(1,5,i)));
  figure(tt);
  if i==1, text(xyzspin(i,1),xyzspin(i,2),xyzspin(i,3),sprintf('%d',Rcspin(1,5,i)));
  else, text(xyzspin(i,1),xyzspin(i,2),xyzspin(i,3),sprintf('%d',Rcspin(1,5,i))); end;
  quiver3(xyzspin(i,1),xyzspin(i,2),xyzspin(i,3),lookspin(i,1),lookspin(i,2),lookspin(i,3));
  pause(0.1);
end
disp('Press Return to continue');
pause
close(tt)
end

%-------------------------------------------
% load geometry data
disp('Loading geometry data ...');
[centroids,numparams,txtparams]=readdat(char(geomfnames(1)));
[norms,numparams,txtparams]=readdat(char(geomfnames(4)));
[neighnum,neighs,numparams,txtparams]=readneighs(char(geomfnames(5)));   % 2/23/05 MWKay
% neighs is a structure
[cells,numparams,txtparams]=readdat(char(geomfnames(3)));

cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
for i=1:size(cells,1)            % VTK starts with zero, matlab starts with one
  for j=1:neighnum(i)
    neighs{i}(j)=neighs{i}(j)+1;
  end
end

%-------------------------------------------
% compute angles b/w normals and view directions
disp('Computing angles between normals and view directions .... ');
langle=zeros(size(norms,1),nviews+1);
disp(sprintf('%d normals and %d views',size(langle,1),nviews));
two_comp=0;   % 1 for 2 components, 0 for three
if ~two_comp
  disp('Computing angles using x,y, and z components.');
  gg=norms*looknspin(anglesi,:)';
  gg=acos(gg);
  gg=real(gg);
  langle(:,1:nviews)=gg.*180/pi;  % x,y, and z components
elseif two_comp 
  disp('Computing angles using x and y components.');
  looknspin2=looknspin(anglesi,1:2);
  norms2=norms(:,1:2);
  gg=norms2*looknspin2';
  gg=acos(gg);
  gg=real(gg);
  langle(:,1:nviews)=gg.*180/pi;  % only x and y components
end

%-------------------------------------------
disp('Finding max angles ....');
for i=1:size(langle,1)
  candidatemax=find(langle(i,1:nviews)==max(langle(i,1:nviews)));
  if length(candidatemax)>1
    langle(i,end)=candidatemax(1);
  elseif length(candidatemax)==1
    langle(i,end)=candidatemax;
  end
end

%--------------------------------------------
% Assign cells to each view
viewi=zeros(size(langle,1),nviews+1);
for i=1:nviews
  ii=anglesi(i);
  disp(sprintf('Finding cells in view %d, %d deg ...',ii,image_angles(i)));
  gg=find(langle(:,i)>minangle);  % find the front face
  if isempty(gg)
    disp(sprintf('Found NO cells in view %d, %d deg ...',ii,image_angles(i)));
  else 
    viewi(gg,i)=1;
    disp(sprintf('Found %d cells in view %d, %d deg ...',length(gg),ii,image_angles(i)));
    disp(sprintf('Smoothing edges of view %d...',ii));
    [viewi(:,i)]=smoothregionbdr(viewi(:,i),neighnum,neighs);
  end
end
clear pvect;

viewi(:,end)=sum(viewi(:,1:nviews),2);

%--------------------------------------------
% projection and mapping
disp('Projection and mapping ...');
txtyour=zeros(size(norms,1),nviews+1).*NaN;
for i=1:nviews
  ii=anglesi(i);
  disp(sprintf('Assigning anatomy from image %s ...',char(imnames(i))));
  disp(sprintf('Using a calibration angle of %d deg.',Rcspin(1,5,ii)));
  gg=find(viewi(:,i));
  if ~isempty(gg);  
    a=imread(sprintf('%s%s/%s',datadir,scanlabel,char(imnames(i))));
    nr=size(a,1); nc=size(a,2);
    
    b=rgb2gray(a);
    b(find(b>=250))=0;
    colornum=mean(b(find(b)));
    if colornum<=50   % then blue and LV is blue
      b(find(b))=1;
      disp('Found region of LV');
    elseif colornum>50 & colornum<=100  % then red and posterior is red
      b(find(b))=2;
      disp('Found posterior region');
    elseif colornum>=125  % then green and anterior is green
      b(find(b))=4;
      disp('Found anterior region');
    end
    
    [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),camera);
    
    out=find(Xigg<1 | Xigg>nc);
    if ~isempty(out)
      Xigg(out)=[];
      Yigg(out)=[];
      gg(out)=[];
    end
    out=find(Yigg<1 | Yigg>nr);
    if ~isempty(out)
      Xigg(out)=[];
      Yigg(out)=[];
      gg(out)=[];
    end
    txtyour(gg,i)=interp2(double(b),Xigg,Yigg);
  end
  disp(sprintf('Completed view %d',ii));
end

%% These lines are for multiple anatomic regions
%for i=1:size(txtyour,1)
%  anatnums=unique(txtyour(i,1:nviews));
%  anatnums=anatnums(find(~isnan(anatnums)));
%  txtyour(i,end)=sum(unique(anatnums));
%end

% This is used to just identify one anatomic region
for i=1:nviews
 gg=find(~isnan(txtyour(:,i)));
 bb=find(txtyour(gg,i));
 txtyour(gg(bb),end)=1;
end

% Smooth the region boarder
disp(sprintf('Smoothing edges of anatomy boarder'));
[txtyour(:,end)]=smoothregionbdr(txtyour(:,end),neighnum,neighs);

% Reassign LV/RV cells according to user specification
if vcellsflag
  eval(sprintf('load %s',vcellsfname));
  [path,dataname,ext]=fileparts(vcellsfname);
  eval(sprintf('thecells=%s(:,1)+1;',dataname));
  thecells=unique(thecells);
  txtyour(thecells,end)=~txtyour(thecells,end);
end


%-------------------------
% Roughly identify anterior and posterior cells from the
% image numbers given in setup_pan_coloranatomy:
% anteriorimnum posteriorimnum

antpostimnum=[anteriorimnum posteriorimnum];
napviews=length(antpostimnum);
apanglesi=antpostimnum+1;
apimage_angles=angles(apanglesi);

% compute angles b/w normals and anterior/posterior views
disp('Computing angles between anterior and posterior view directions .... ');
aplangle=zeros(size(norms,1),napviews+1);
disp(sprintf('%d normals and %d views',size(aplangle,1),napviews));
two_comp=0;   % 1 for 2 components, 0 for three
if ~two_comp
  disp('Computing angles using x,y, and z components.');
  gg=norms*looknspin(apanglesi,:)';
  gg=acos(gg);
  gg=real(gg);
  aplangle(:,1:napviews)=gg.*180/pi;  % x,y, and z components
elseif two_comp 
  disp('Computing angles using x and y components.');
  looknspin2=looknspin(apanglesi,1:2);
  norms2=norms(:,1:2);
  gg=norms2*looknspin2';
  gg=acos(gg);
  gg=real(gg);
  aplangle(:,1:napviews)=gg.*180/pi;  % only x and y components
end

disp('Finding max angles for anterior/posterior views ....');
for i=1:size(aplangle,1)
  candidatemax=find(aplangle(i,1:napviews)==max(aplangle(i,1:napviews)));
  if length(candidatemax)>1
    aplangle(i,end)=candidatemax(1);
  elseif length(candidatemax)==1
    aplangle(i,end)=candidatemax;
  end
end

% Assign cells to either anterior or posterior views
apviewi=zeros(size(aplangle,1),napviews+1);
apviewi(:,end)=apviewi(:,end).*NaN;
for i=1:napviews
  ii=apanglesi(i);
  disp(sprintf('Finding cells in view %d, %d deg ...',ii,apimage_angles(i)));
  gg=find(aplangle(:,i)>minangle);  % find the front face
  if isempty(gg)
    disp(sprintf('Found NO cells in view %d, %d deg ...',ii,apimage_angles(i)));
  else 
    apviewi(gg,i)=1;
    disp(sprintf('Found %d cells in view %d, %d deg ...',length(gg),ii,apimage_angles(i)));
    disp(sprintf('Smoothing edges of view %d...',ii));
    [apviewi(:,i)]=smoothregionbdr(apviewi(:,i),neighnum,neighs);
  end
end
clear pvect;

% Assign anterior cells
gg=find(apviewi(:,1)==1);
apviewi(gg,end)=1;                  % 1 assigned to anterior cells

% Posterior cells are zero
gg=find(apviewi(:,2)==1);
apviewi(gg,end)=0;                  % 0 assigned to posterior cells

% Assign cells in both views to view with max angle
ggg=find(apviewi(gg,1)==1);
if ~isempty(ggg)
  apviewi(gg(ggg),end)=aplangle(gg(ggg),end)-1;       
end

% Assign cells not yet assigned
% Assign them to view with max angle
gg=find(isnan(apviewi(:,end)));
if ~isempty(gg)
  apviewi(gg,end)=aplangle(gg,end)-1;
end

% Reassign ant/post cells according to user specification
if acellsflag
  eval(sprintf('load %s',acellsfname));
  [path,dataname,ext]=fileparts(acellsfname);
  eval(sprintf('thecells=%s(:,1)+1;',dataname));
  thecells=unique(thecells);
  apviewi(thecells,end)=~apviewi(thecells,end);
end

% Smooth the region boarders
disp(sprintf('Smoothing edges of ant/post boarder'));
[apviewi(:,end)]=smoothregionbdr(apviewi(:,end),neighnum,neighs);

% Identify basal and apical cells
apexbase=zeros(size(txtyour,1),1);
postmidpt=median(centroids(find(apviewi(:,end)==0),3));
gg=find(centroids(:,3)>postmidpt);
apexbase(gg)=1;

% Smooth the apex/base border
disp(sprintf('Smoothing edges of apex/base boarder'));
[apexbase]=smoothregionbdr(apexbase,neighnum,neighs);

% Assign an anatomy code
anatcode=zeros(size(centroids,1),1).*NaN;
for i=1:size(centroids,1)
  binnum=[txtyour(i,end) apviewi(i,end) apexbase(i)];
  anatcode(i)=dot(binnum,[1 2 4]);
end


%--------------------------------------------
% identify the edge cells for the LV colored images
disp('Finding points along edges ...');
edges=int8(zeros(size(viewi,1),1));
for j=1:size(edges,1)
  if neighnum(j)==3
    if (viewi(neighs{j}(1),end)-viewi(neighs{j}(2),end))~=0
      edges(j)=1;
    elseif (viewi(neighs{j}(1),end)-viewi(neighs{j}(3),end))~=0
      edges(j)=1;
    elseif (viewi(neighs{j}(2),end)-viewi(neighs{j}(3),end))~=0
      edges(j)=1;
    end
  elseif neighnum(j)<3 % if cell has less than 3 neighbors then on an edge
      edges(j)=1;
  end
end

% find points along the edges (find the indicies)
edgesi=find(edges);
edgepts3col=zeros(length(edgesi),3).*NaN;
for i=1:length(edgesi)
  % lazy approach:
  % simply check all the edges of each previously identified cell
  % then 'unique' the list
  if neighnum(edgesi(i))==3
    if (viewi(edgesi(i),end)-viewi(neighs{edgesi(i)}(1),end))~=0
      edgepts3col(i,:)=[cells(edgesi(i),1), cells(edgesi(i),2), NaN];
    elseif (viewi(edgesi(i),end)-viewi(neighs{edgesi(i)}(2),end))~=0
      edgepts3col(i,:)=[NaN, cells(edgesi(i),2), cells(edgesi(i),3)];
    elseif (viewi(edgesi(i),end)-viewi(neighs{edgesi(i)}(3),end))~=0
      edgepts3col(i,:)=[cells(edgesi(i),1), NaN, cells(edgesi(i),3)];
    end
  elseif neighnum(edgesi(i))==2  % Take the 2 verticies the neighbors do not have in common
      oneedge=intersect(cells(neighs{edgesi(i)}(1),:),cells(edgesi(i),:));
      twoedge=intersect(cells(neighs{edgesi(i)}(2),:),cells(edgesi(i),:));
      edgepts3col(i,:)=[setdiff(oneedge,twoedge) setdiff(twoedge,oneedge) NaN];
  elseif neighnum(edgesi(i))==1  % Take all verticies
      edgepts3col(i,:)=cells(edgesi(i),:);
  end
end
edgepts=[edgepts3col(:,1); edgepts3col(:,2); edgepts3col(:,3)];
edgepts=edgepts(find(~isnan(edgepts)));
edgepts=unique(edgepts);   % these are indicies into the pts array!

%---------------------------------------------
% save texture for each cell to a cdf file -- Implemented on 7/29/03
txtyourfnamef=sprintf('%s.%s',char(txtyourfname{1}),char(txtyourfname{2}));
fid=fopen(txtyourfnamef,'w','b');
% write the header
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',size(centroids,1));
fprintf(fid,'delta-t=1.0\n');
fprintf(fid,'n_samples=4\n');
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=%s\n',date);
fprintf(fid,'data_source=%s,%s\n',study,scanlabel);
dat_precision=sprintf('uint%d',8*4);
fwrite(fid,[0 txtyour(:,end)'].*1000,dat_precision);   % this is the LV data
fwrite(fid,[0 apviewi(:,end)'].*1000,dat_precision);   % this is the ant/post data
fwrite(fid,[0 apexbase'].*1000,dat_precision);         % this is the apex/base data
fwrite(fid,[0 anatcode'].*1000,dat_precision);         % this is the 8 region anatomy code
fclose(fid);
disp(sprintf('Saved texture in %s',txtyourfnamef));

%---------------------------------------------
% save anatcode for each cell to a cdf file -- 2/24/05
txtyourfnamef=sprintf('%s.%s',char(txtyourfname2{1}),char(txtyourfname2{2}));
fid=fopen(txtyourfnamef,'w','b');
% write the header
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',size(centroids,1));
fprintf(fid,'delta-t=1.0\n');
fprintf(fid,'n_samples=1\n');
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=%s\n',date);
fprintf(fid,'data_source=%s,%s\n',study,scanlabel);
dat_precision=sprintf('uint%d',8*4);
fwrite(fid,[0 anatcode'].*1000,dat_precision);         % this is the 8 region anatomy code
fclose(fid);
disp(sprintf('Saved texture in %s',txtyourfnamef));

%---------------------------------------------
% save indicies of points on the edge. THESE ARE POINT INDICIES, NOT CELL INDICIES!
edgesfnamef=sprintf('%s.%s',char(edgesfname{1}),char(edgesfname{2}));
fid=fopen(edgesfnamef,'w');
fprintf(fid,'3\n');
fprintf(fid,'Edge point indicies, created by pan_textureIM.m\n');
fprintf(fid,'int\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(edgepts,1));
fprintf(fid,'%d\n',size(edgepts,2));
fwrite(fid,edgepts-1,'int');
fclose(fid);
disp(sprintf('Saved edge point indicies in %s',edgesfnamef));

%---------------------------------------------
% save view for each cell
viewsfnamef=sprintf('%s.%s',char(viewsfname{1}),char(viewsfname{2}));
fid=fopen(viewsfnamef,'w');
fprintf(fid,'3\n');
fprintf(fid,'View data, created by pan_textureIM.m\n');
fprintf(fid,'float\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(viewi(:,end),1));
fprintf(fid,'%d\n',size(viewi(:,end),2));
fwrite(fid,viewi(:,end),'float');
fclose(fid);
disp(sprintf('Saved view numbers in %s',viewsfnamef));

