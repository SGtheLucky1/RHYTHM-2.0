close all
clear all
disp(sprintf('\npan_textureIM7.m \nVersion: 11/18/02 \n\n'));
disp('Press Return ...'); pause;

expname=input('Experiment name [ie, map05_07_15_03]: ','s');
scanname=input('Scan label [ie, scan11]: ','s');

%-----------------------------------
% user defined parameters

% INPUT FILES
spincalfnames={'cal/Rc.dat';       % Rc first
               'cal/Par.dat';};    % Par second
	       
image_fnames=sprintf('%s/%s__fnames.txt',scanname,scanname);

geomfnames={sprintf('%s/scan4_d_centroids.dat',scanname);
            sprintf('%s/scan4_d_pts.dat',scanname);
            sprintf('%s/scan4_d_cells.dat',scanname);
	    sprintf('%s/scan4_d_smoothnorms.dat',scanname);
	    sprintf('%s/scan4_d_neighs.dat',scanname);};

silhfname=sprintf('%s/silhs1.mat',scanname);

% OUTPUT FILES	    
txtyourfname_prefix=sprintf('%s/%s_anatomy',scanname,scanname);
txtyourfname_postfix='cdf';
edgesfname_prefix=sprintf('%s/%s_anatomyedges',scanname,scanname);
edgesfname_postfix='dat';
viewsfname_prefix=sprintf('%s/%s_anatomyviews',scanname,scanname);
viewsfname_postfix='dat';

%image_angles=[0 45 90 135 180 225 270 315];
image_angles=[0 90 180 270];
%image_angles=[0:20:340];
%image_angles=[0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340];

minangle=95;   % front face angles must be greater than this, use 180 to force max angle assignment
minedge=0.90;
edgefiltval=100;

%------------------------------------
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
      cname='watec_with_f8.5';
      done=1;
    case 2
      cname='watec_with_f12.5';
      done=1;
    case 3
      cname='E4300';
      done=1;   
    case 4
      cname='interpolated_dalsa';
      done=1;
  end
end

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
reorder='g';
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
anglesi=image_angles./dtheta+1;   % image_angles defined by user
nviews=length(image_angles);

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

%-------------------------------------------
% load geometry data
disp('Loading geometry data ...');
[centroids,numparams,txtparams]=readdat(char(geomfnames(1)));
[norms,numparams,txtparams]=readdat(char(geomfnames(4)));
[neighs,numparams,txtparams]=readdat(char(geomfnames(5)));
[cells,numparams,txtparams]=readdat(char(geomfnames(3)));

cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
neighs=neighs+1;                 % VTK starts with zero, matlab starts with one

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
% load image silhouettes and compute edge weights
disp(sprintf('Loading image silhouettes from %s',silhfname));
load(silhfname);
nr=size(silh1,1);
nc=size(silh1,2);
disp('Computing edge weights ...');
ledge=zeros(size(norms,1),nviews+1);
for i=1:nviews
  ii=anglesi(i);
  disp(sprintf('Computing weights for view %d, %d deg ...',ii,image_angles(i)));
  gg=find(langle(:,i)>90);  % find the front face
  [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),cname);
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
  [conesilh,kernel]=conefilt(silh1(:,:,ii),edgefiltval,1);
  ledge(gg,i)=interp2(conesilh,Xigg,Yigg); 
  clear conesilh;
end
clear silh1;

%-------------------------------------------
disp('Finding max edge weights ....');
for i=1:size(ledge,1)
  gg=find(ledge(i,1:nviews)==max(ledge(i,1:nviews)));
  if (~isempty(gg) & length(gg)==1)
    ledge(i,end)=gg;            % max ledge
  elseif length(unique(langle(i,gg)))==1
    ledge(i,end)=gg(1);
  else
    ledge(i,end)=find(langle(i,1:nviews)==max(langle(i,gg)));  % max angle at max ledge
  end
end


%--------------------------------------------
% Assign cells to each view
viewi=zeros(size(langle,1),nviews+1);
for i=1:nviews
  ii=anglesi(i);
  disp(sprintf('Finding cells in view %d, %d deg ...',ii,image_angles(i)));
  gg=find(langle(:,i)>minangle & ledge(:,i)>minedge);  % find the front face
  if isempty(gg)
    disp(sprintf('Found NO cells in view %d, %d deg ...',ii,image_angles(i)));
  else 
    viewi(gg,i)=1;
    disp(sprintf('Found %d cells in view %d, %d deg ...',length(gg),ii,image_angles(i)));
    disp(sprintf('Smoothing edges of view %d...',ii));
    pvect=zeros(size(viewi,1),2);
    pvect(:,1)=viewi(:,i);
    pvect(:,2)=viewi(:,i);
    pkount=0;
    pdiffnum=1;
    while (pdiffnum>0)
      pkount=pkount+1;
      p_now=mod(pkount,2)+1;      % 2 first, then 1
      p_prev=mod(pkount+1,2)+1;
      for jj=1:size(pvect,1)
	if (pvect(neighs(jj,1),p_now)-pvect(neighs(jj,2),p_now))==0
	  pvect(jj,p_now)=pvect(neighs(jj,1),p_now);
	elseif (pvect(neighs(jj,1),p_now)-pvect(neighs(jj,3),p_now))==0
	  pvect(jj,p_now)=pvect(neighs(jj,1),p_now);
	elseif (pvect(neighs(jj,2),p_now)-pvect(neighs(jj,3),p_now))==0
	  pvect(jj,p_now)=pvect(neighs(jj,2),p_now);
	end
      end
      pdiffnum=length(find(pvect(:,1)-pvect(:,2)~=0));
      disp(sprintf('Pass %d, Reassigned %d cells',pkount,pdiffnum)); 
      pvect(:,p_prev)=pvect(:,p_now);
    end
    viewi(:,i)=pvect(:,1);   % smoothed edges
  end
end
clear pvect;

%--------------------------------------------
% load image filenames
disp(sprintf('Loading image filenames from %s',image_fnames));
imnames=cell(ntheta,1);
imnamesfid=fopen(image_fnames,'r');
for i=1:ntheta
  imnames(i)=cellstr(fgets(imnamesfid));
end
fclose(imnamesfid);

%--------------------------------------------
% projection and mapping
disp('Projection and mapping ...');
txtyour=zeros(size(norms,1),nviews+1).*NaN;
for i=1:nviews
  ii=anglesi(i);
  disp(sprintf('Computing texture from image %s ...',char(imnames(ii))));
  disp(sprintf('Using a calibration angle of %d deg.',Rcspin(1,5,ii)));
  gg=find(viewi(:,i));
  if ~isempty(gg);    
    a=gray2ind(readpgm(char(imnames(ii))),256);
    b=zeros(size(a));
    c=nan.*b;
    [Xi,Yi]=pred2([centroids(:,1) centroids(:,2) centroids(:,3) ones(size(centroids(:,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),cname);
    out=find(Xi<1 | Xi>nc);
    if ~isempty(out)
      Xi(out)=[];
      Yi(out)=[];
    end
    out=find(Yi<1 | Yi>nr);
    if ~isempty(out)
      Xi(out)=[];
      Yi(out)=[];
    end
    b(sub2ind(size(a),round(Yi),round(Xi)))=1;
    for j=1:size(b,2)
      maxc=max(find(b(:,j)));
      minc=min(find(b(:,j)));
      c(minc:maxc,j)=a(minc:maxc,j);
    end
    
    if 0
      % testing stuff here
      ci=find(~isnan(c));
      cvalues=c(ci);
      minc1=min(cvalues);
      maxc1=max(cvalues);
      midc1=(maxc1-minc1)/2;
      %w=abs(op_norm((cvalues-midc1),0.5,1));
      w=op_norm(-(cvalues-midc1).^3,0.5,1.5);
      c2=c;
      newcvalues=imlinearshift2(cvalues.*w,min(cvalues.*w),max(cvalues.*w),1,256);
      c2(ci)=newcvalues;
      minc2=min(newcvalues);
      maxc2=max(newcvalues);
      midc2=(maxc2-minc2)/2;
      c2(find(isnan(c2)))=midc2;
     %c=imlinearshift(double(histeq(uint8(c),[10:110])),10,160,1,256);
     %c=imlinearshift(double(a),100,250,1,256);  % was 160,245,1,256
     %c=imlinearshift(double(histeq(uint8(c),256)),0,120,1,256);
    end
     
     %c2=imlinearshift2(double(a),50,260,1,256);
     c2=double(a);
   
    [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),cname);
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
    txtyour(gg,i)=interp2(c2,Xigg,Yigg);   
  end
  disp(sprintf('Completed view %d',ii));
end
clear a b c c2 Xi Yi Xigg Yigg;

% Assign one texture for each cell
disp('Assigning texture ...');
for i=1:size(txtyour,1)
  gg=find(~isnan(txtyour(i,1:nviews)));
  if ~isempty(gg)&(length(gg)==1)
      txtyour(i,end)=txtyour(i,gg);
      viewi(i,end)=gg;
  elseif ~isempty(gg)
      txtyour(i,end)=sum(txtyour(i,gg).*ledge(i,gg))/sum(ledge(i,gg)); % Weighted average
      viewi(i,end)=nviews+1;
    % txtyour(i,end)=mean(txtyour(i,gg));
    % viewi(i,end)=mean(gg);
    % txtyour(i,end)=min(txtyour(i,gg));   % the darkest?
    % txtyour(i,end)=median(txtyour(i,gg));
  end  
end

%--------------------------------------------------
% assign texture to NaN sites
disp('Assigning texture to empty sites ...');
for i=1:nviews 
  ii=anglesi(i);
  gg=find(isnan(txtyour(:,end))&langle(:,i)>90&ledge(:,i)>0.25);  % data is bad if it doesn't meet this criteria
  disp(sprintf('View %d, found %d empty sites',ii,length(gg)));
  if ~isempty(gg) 
    a=gray2ind(readpgm(char(imnames(ii))),256);  
    b=zeros(size(a));
    c=b+256; 
    [Xi,Yi]=pred2([centroids(:,1) centroids(:,2) centroids(:,3) ones(size(centroids(:,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),cname);
    out=find(Xi<1 | Xi>nc);
    if ~isempty(out)
      Xi(out)=[];
      Yi(out)=[];
    end
    out=find(Yi<1 | Yi>nr);
    if ~isempty(out)
      Xi(out)=[];
      Yi(out)=[];
    end
    b(sub2ind(size(a),round(Yi),round(Xi)))=1;
    for j=1:size(b,2)
      maxc=max(find(b(:,j)));
      minc=min(find(b(:,j)));
      c(minc:maxc,j)=a(minc:maxc,j);
    end
   
    % use same command here as above
    %c2=imlinearshift2(double(a),50,260,1,256);
    c2=double(a);
   
    [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcspin(:,1:4,ii),Parspin(1:8,ii),cname);
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
    txtyour(gg,i)=interp2(c2,Xigg,Yigg); 
  end  
  disp(sprintf('Completed view %d',ii));
end  
clear a b c Xi Yi Xigg Yigg;
  
% Now combine textures
gg=find(isnan(txtyour(:,end)));
for i=1:length(gg)
  ss=find(~isnan(txtyour(gg(i),1:nviews)));
  if ~isempty(ss)&(length(ss)==1)
      txtyour(gg(i),end)=txtyour(gg(i),ss);
      viewi(gg(i),end)=ss;
  elseif ~isempty(ss)
    % txtyour(gg(i),end)=mean(txtyour(gg(i),ss));
      txtyour(gg(i),end)=sum(txtyour(gg(i),ss).*ledge(gg(i),ss))/sum(ledge(gg(i),ss)); % Weighted average
      viewi(gg(i),end)=nviews+2;
    % viewi(gg(i),end)=mean(ss);
  end  
end


% Improve contrast for nviews+1 cells
%gg=find(viewi(:,end)==nviews+1);
%txtyour(gg,end)=dat_norm(txtyour(gg,end),min(txtyour(:,end)),max(txtyour(:,end)));

% Improve contrast for nviews+2 cells
%gg=find(viewi(:,end)==nviews+2);
%txtyour(gg,end)=dat_norm(txtyour(gg,end),min(txtyour(:,end)),max(txtyour(:,end)));

% Improve contrast
gg=find(viewi(:,end));
txtyour(gg,end)=op_norm(txtyour(gg,end),1,256);

%--------------------------------------------
% identify the edge cells
disp('Finding points along edges ...');
edges=int8(zeros(size(viewi,1),1));
for j=1:size(edges,1)
  if (viewi(neighs(j,1),end)-viewi(neighs(j,2),end))~=0
    edges(j)=1;
  elseif (viewi(neighs(j,1),end)-viewi(neighs(j,3),end))~=0
    edges(j)=1;
  elseif (viewi(neighs(j,2),end)-viewi(neighs(j,3),end))~=0
    edges(j)=1;
  end
end

% find points along the edges (find the indicies)
edgesi=find(edges);
edgepts2col=zeros(length(edgesi),2);
for i=1:length(edgesi)
  if (viewi(edgesi(i),end)-viewi(neighs(edgesi(i),1),end))~=0
    edgepts2col(i,:)=cells(edgesi(i),[1 2]);
  elseif (viewi(edgesi(i),end)-viewi(neighs(edgesi(i),2),end))~=0
    edgepts2col(i,:)=cells(edgesi(i),[2 3]);
  elseif (viewi(edgesi(i),end)-viewi(neighs(edgesi(i),3),end))~=0
    edgepts2col(i,:)=cells(edgesi(i),[3 1]);
  end
end
edgepts=[edgepts2col(:,1); edgepts2col(:,2)];
edgepts=unique(edgepts);   % these are indicies into the pts array!

    
%--------------------------------------------------
% Get rid of sites on the top 
top_thresh1=25; % angle, degrees from z axis
top_thresh2=175; % angle, degrees from z axis
thetaz=acos(norms(:,3)).*180/pi;
thetaz=real(thetaz);
gg=find(thetaz<=top_thresh1 | thetaz>=top_thresh2);
txtyour(gg,end)=NaN;
viewi(gg,end)=0;

%--------------------------------------------------
% set low sites to 1
%mapthresh=0;
%gg=find(txtyour(:,end)<=mapthresh);
%txtyour(gg,end)=1;
 
%---------------------------------------------
% save texture for each cell to a cdf file -- Implemented on 7/29/03
txtyourfname=sprintf('%s.%s',txtyourfname_prefix,txtyourfname_postfix);
fid=fopen(txtyourfname,'w','b');
% write the header
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',size(centroids,1));
fprintf(fid,'delta-t=0.00\n');
fprintf(fid,'n_samples=1\n');
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=%s\n',date);
fprintf(fid,'data_source=%s,%s\n',expname,image_fnames);
dat_precision=sprintf('uint%d',8*4);
fwrite(fid,[0 txtyour(:,end)'].*1000,dat_precision);
fclose(fid);
disp(sprintf('Saved texture in %s',txtyourfname));

%---------------------------------------------
% save indicies of points on the edge. THESE ARE POINT INDICIES, NOT CELL INDICIES!
edgesfname=sprintf('%s.%s',edgesfname_prefix,edgesfname_postfix);
fid=fopen(edgesfname,'w');
fprintf(fid,'3\n');
fprintf(fid,'Edge point indicies, created by pan_textureIM.m\n');
fprintf(fid,'int\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(edgepts,1));
fprintf(fid,'%d\n',size(edgepts,2));
fwrite(fid,edgepts-1,'int');
fclose(fid);
disp(sprintf('Saved edge point indicies in %s',edgesfname));

%---------------------------------------------
% save view for each cell
viewsfname=sprintf('%s.%s',viewsfname_prefix,viewsfname_postfix);
fid=fopen(viewsfname,'w');
fprintf(fid,'3\n');
fprintf(fid,'View data, created by pan_textureIM.m\n');
fprintf(fid,'float\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(viewi(:,end),1));
fprintf(fid,'%d\n',size(viewi(:,end),2));
fwrite(fid,viewi(:,end),'float');
fclose(fid);
disp(sprintf('Saved view numbers in %s',viewsfname));

