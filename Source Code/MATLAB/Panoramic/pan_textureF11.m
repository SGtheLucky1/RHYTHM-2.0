close all
clear all

%% - added fluorescence averaging within surface mesh cells. 6/4/04, MWKay
%% - removed signal normalization for subsequent signal quality measurement,
%%   filtering and normalization in pan_postprocess.m. 11/09/04, MWKay
%% - added setup_pan_texture.m. 11/15/04, MWKay
%% - added output of texturing information. 11/23/04, MWKay
%% - added pixel fine-tuning of mapping data positions on geometry. 12/03/04, MWKay

disp(sprintf('\npan_textureF11.m \nVersion: 12/03/04 \n\n'));

disp('Press Return to read setup_pan_texture.m');
disp('(located in current dir) and begin processing ...'); pause;

% Read user defined parameters and filenames
setup_pan_texture

%-----------------------------------
% Load mapping camera calibration stuff

nviews=size(mapcalfnames,1);

Rcmap=NaN.*ones(4,4,nviews);
Parmap=NaN.*ones(8,nviews);
xyzmap=NaN.*ones(nviews,4);
for i=1:nviews
  pospar=load(char(mapcalfnames(i)));
  wa=pospar(4)*pi/180;  % omega, Wx rotation (world basis, x axis) 
  pa=pospar(5)*pi/180;  % psi, y rotation  (world basis, y axis)  
  ra=pospar(6)*pi/180;  % kappa, z rotation (world basis, z axis) 
  cw=cos(wa); sw=sin(wa);
  cp=cos(pa); sp=sin(pa);
  cr=cos(ra); sr=sin(ra);
  R=zeros(3,3);

  R(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]';             % build the 4x4 R matrix
  R(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
  R(:,3)=[-sp cp*sw cp*cw]';
  Rcmap(:,:,i)=[[R'; [0 0 0]] [pospar(1) pospar(2) pospar(3) 1]'];  % store each R in a 3d Rc matrix for easy referencing later

  xyzmap(i,:)=(inv(Rcmap(:,:,i))*[0 0 0 1]')';   % compute xyz positions of orbit from calibration matricies
  Parmap(:,i)=pospar(7:14);
end

lookmap=ones(size(Rcmap,3),4);
looknmap=ones(size(Rcmap,3),3);
for i=1:size(Rcmap,3)
  lookmap(i,:)=(inv(Rcmap(1:4,1:4,i))*[0 0 1 1]')';       % camera view direction
  lookmap(i,1:3)=lookmap(i,1:3)-xyzmap(i,1:3);            % camera view direction
  looknmap(i,:)=lookmap(i,1:3)./norm(lookmap(i,1:3));     % normalized view direction
end

%-------------------------------------------
% load geometry data
disp('Loading geometry data ...');
[centroids,numparams,txtparams]=readdat(char(geomfnames(1)));
[norms,numparams,txtparams]=readdat(char(geomfnames(4)));
[neighs,numparams,txtparams]=readdat(char(geomfnames(5)));
[cells,numparams,txtparams]=readdat(char(geomfnames(3)));
[pts,numparams,txtparams]=readdat(char(geomfnames(2)));
cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
neighs=neighs+1;                 % VTK starts with zero, matlab starts with one

%--------------------------------------------
% load geometry masks (these should all be 128x128)
disp('Loading geometry masks ...');
disp(sprintf('Loading geommask1: %s',char(maskfnames(1))));
geommask1=imread(char(maskfnames(1)));
disp(sprintf('Loading geommask2: %s',char(maskfnames(2))));
geommask2=imread(char(maskfnames(2)));
disp(sprintf('Loading geommask3: %s',char(maskfnames(3))));
geommask3=imread(char(maskfnames(3)));
disp(sprintf('Loading geommask4: %s',char(maskfnames(4))));
geommask4=imread(char(maskfnames(4)));

%---------------------------------------
% load mapping data
disp('Loading the mapping data ...');
pix=zeros(4,4);
% data1
disp(sprintf('reading %s',char(datafnames(1))));
[data1,bkgrnd1,pix(1,:)]=loadprocsif(char(datafnames(1)),1);
% data2
disp(sprintf('reading %s',char(datafnames(2))));
[data2,bkgrnd2,pix(2,:)]=loadprocsif(char(datafnames(2)),1);
% data3
disp(sprintf('reading %s',char(datafnames(3))));
[data3,bkgrnd3,pix(3,:)]=loadprocsif(char(datafnames(3)),1);
% data4
disp(sprintf('reading %s',char(datafnames(4))));
[data4,bkgrnd4,pix(4,:)]=loadprocsif(char(datafnames(4)),1);

nr=unique(pix(:,3)-pix(:,4)+1);
nc=unique(pix(:,2)-pix(:,1)+1);
if length(nr)>1 | length(nc)>1
  disp('Number of rows or number of columns does not match for each dataset.');
  disp('This is going to be a problem!');
  disp('Bailing out...');
  return
else
  disp('For each of the 4 datasets,');
  disp(sprintf('Num of rows=%d, Num of columns=%d',nr,nc));
end

%---------------------------------------
% Position and crop the masks according to user input

disp('Position mask for view 1');
[ri,ci]=find(geommask1);
ri=ri-pix(1,4)+1; ci=ci-pix(1,1)+1;
[newri,newci]=usershiftind(bkgrnd1,ri,ci);
mask1=zeros(size(bkgrnd1));
mask1(sub2ind(size(mask1),newri,newci))=1;

disp('Position mask for view 2');
[ri,ci]=find(geommask2);
ri=ri-pix(2,4)+1; ci=ci-pix(2,1)+1;
[newri,newci]=usershiftind(bkgrnd2,ri,ci);
mask2=zeros(size(bkgrnd2));
mask2(sub2ind(size(mask2),newri,newci))=1;

disp('Position mask for view 3');
[ri,ci]=find(geommask3);
ri=ri-pix(3,4)+1; ci=ci-pix(3,1)+1;
[newri,newci]=usershiftind(bkgrnd3,ri,ci);
mask3=zeros(size(bkgrnd3));
mask3(sub2ind(size(mask3),newri,newci))=1;

disp('Position mask for view 4');
[ri,ci]=find(geommask4);
ri=ri-pix(4,4)+1; ci=ci-pix(4,1)+1;
[newri,newci]=usershiftind(bkgrnd4,ri,ci);
mask4=zeros(size(bkgrnd4));
mask4(sub2ind(size(mask4),newri,newci))=1;

%---------------------------------------
% Detrend the data
disp('detrending dataset 1 ...');
[gr,gc]=find(mask1);
data1=op_detrend(data1,gr,gc);
disp('detrending dataset 2 ...');
[gr,gc]=find(mask2);
data2=op_detrend(data2,gr,gc);
disp('detrending dataset 3 ...');
[gr,gc]=find(mask3);
data3=op_detrend(data3,gr,gc);
disp('detrending dataset 4 ...');
[gr,gc]=find(mask4);
data4=op_detrend(data4,gr,gc);


%-------------------------------------------
% compute angles b/w normals and view directions
disp('Computing angles between normals and view directions .... ');
langle=zeros(size(norms,1),nviews+1); 
disp(sprintf('%d normals and %d views',size(langle,1),nviews));
two_comp=0;   % 1 for 2 components, 0 for three
if ~two_comp
  disp('Computing angles using x,y, and z components.');
  gg=norms*looknmap';
  gg=acos(gg);
  gg=real(gg);
  langle(:,1:nviews)=gg.*180/pi;  % x,y, and z components
elseif two_comp 
  disp('Computing angles using x and y components.');
  looknmap2=looknmap(:,1:2);
  norms2=norms(:,1:2);
  gg=norms2*looknmap2';
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
% Compute edge weights
disp('Computing edge weights ...');
ledge=zeros(size(norms,1),nviews+1);
for i=1:nviews
  disp(sprintf('Computing weights for view %d ...',i));
  gg=find(langle(:,i)>90);  % find the front face
  [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
  Xigg=Xigg-pix(i,1)+1; Yigg=Yigg-pix(i,4)+1;
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
  com=sprintf('[conesilh,kernel]=conefilt(mask%d,edgefiltval,1);',i);
  eval(com);
  ledge(gg,i)=interp2(conesilh,Xigg,Yigg); 
  clear conesilh;
end

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
% Assign cells to views: DIRECT, OVERLAP, and EDGE

% This sets up the DIRECT views and smooths the edges of those views
% The next step finds the OVERLAP of the DIRECT views

% viewi(:,1-nviews)=0 or 1, indicating whether the cell was mapped by a particular view
% viewi(:,nviews+1)=integer b/w 1 and nviews+2. This number indicates the assigned view. 
%      nviews+1 is overlap, nviews+2 is edge

viewi=zeros(size(langle,1),nviews+1);
for i=1:nviews
  disp(sprintf('Finding cells in DIRECT view %d ...',i));
  gg=find(langle(:,i)>=minangleDIRECT & ledge(:,i)>=minedgeDIRECT);
  if isempty(gg)
    disp(sprintf('Found NO cells in view %d ...',i));
  else 
    viewi(gg,i)=1;
    disp(sprintf('Found %d cells in DIRECT view %d ...',length(gg),i));
    disp(sprintf('Smoothing edges of DIRECT view %d...',i));
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

disp(sprintf('Assigning OVERLAP and EDGE views...'));
for i=1:size(viewi,1)
  gg=find(viewi(i,1:nviews));
  if ~isempty(gg);
    if length(gg)==1
      viewi(i,nviews+1)=gg;          % DIRECT
    else
      viewi(i,nviews+1)=nviews+1;    % OVERLAP
    end
  else
    gg=find(langle(i,1:nviews)>=90&ledge(i,1:nviews)>=minedgeEDGE);   
    if ~isempty(gg) 
      viewi(i,nviews+1)=nviews+2;    % EDGE
      viewi(i,gg)=1;  % average these views
    end
  end
end

disp(sprintf('Finding indicies of cells with multiple views...'));
mults=sum(viewi(:,1:nviews),2);
multsover1=find(mults>1);
multsii=zeros(length(find(viewi(multsover1,1:nviews))),1);
multsjj=multsii;
kount=0;
for i=1:size(viewi,1)
  gg=find(viewi(i,1:nviews));
  if ~isempty(gg) & length(gg)>1
    for j=1:length(gg)
      kount=kount+1;
      multsii(kount)=i;
      multsjj(kount)=gg(j);
    end
  end
end


%--------------------------------------------
% projection and mapping

% txtkernel stores pixel locations and weights
% txtkernel{:,:,1} is array of pixel indicies
% txtkernel{:,:,2} is array of pixel weights
txtkernel=cell(size(viewi,1),nviews,2); 
disp('Determining pixel kernels for each cell in each view ...');
for i=1:nviews
  disp(sprintf('Doing this for view %d ...',i));
  gg=find(viewi(:,i));
  if ~isempty(gg);
    [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
    Xigg=Xigg-pix(i,1)+1; Yigg=Yigg-pix(i,4)+1;
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
    V1=zeros(size(gg,1),3); V2=V1; V3=V1;
    V1=[pts(cells(gg,1),1) pts(cells(gg,1),2) pts(cells(gg,1),3)];
    V2=[pts(cells(gg,2),1) pts(cells(gg,2),2) pts(cells(gg,2),3)];
    V3=[pts(cells(gg,3),1) pts(cells(gg,3),2) pts(cells(gg,3),3)];
    [XV1gg,YV1gg]=pred2([V1(:,1) V1(:,2) V1(:,3) ones(size(V1,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
     XV1gg=XV1gg-pix(i,1)+1;YV1gg=YV1gg-pix(i,4)+1;
    [XV2gg,YV2gg]=pred2([V2(:,1) V2(:,2) V2(:,3) ones(size(V2,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
     XV2gg=XV2gg-pix(i,1)+1;YV2gg=YV2gg-pix(i,4)+1;
    [XV3gg,YV3gg]=pred2([V3(:,1) V3(:,2) V3(:,3) ones(size(V3,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
     XV3gg=XV3gg-pix(i,1)+1;YV3gg=YV3gg-pix(i,4)+1;
    for j=1:size(gg,1)
      % Define the triangle's bounding square and round up to pixel edges
      % To do this, add or subtract 0.5 to grid the integer-valued pixels  
      VXmax=ceil(max([XV1gg(j) XV2gg(j) XV3gg(j)])-0.5)+0.5; 
      VXmin=floor(min([XV1gg(j) XV2gg(j) XV3gg(j)])+0.5)-0.5;
      VYmax=ceil(max([YV1gg(j) YV2gg(j) YV3gg(j)])-0.5)+0.5;
      VYmin=floor(min([YV1gg(j) YV2gg(j) YV3gg(j)])+0.5)-0.5;
      if VXmax>nc+0.5; VXmax=nc+0.5; end;
      if VXmin<0.5; VXmin=0.5; end;
      if VYmax>nr+0.5; VYmax=nr+0.5; end;
      if VYmin<0.5; VYmin=0.5; end;
      % Define the pixel grid
      [pvX,pvY]=meshgrid([VXmin:1:VXmax],[VYmin:1:VYmax]);
      intarea=zeros(size(pvX,1)-1,size(pvX,2)-1);
      Q=[XV1gg(j),YV1gg(j);
	 XV2gg(j),YV2gg(j);
	 XV3gg(j),YV3gg(j)];
      Q=checkpoly(Q*100)/100;  % multiply and divide by 100 to avoid numerical imprecision with checkpoly.m
      txtyour(gg(j),i)=0;
    % plot(Q(:,1),Q(:,2),'k*');
      for pi=1:size(pvX,1)-1
	for pj=1:size(pvY,2)-1
	  % assign P for cw rotation (yes, cw! it works..)
	  P=[pvX(pi,pj+1) pvY(pi,pj+1);
	     pvX(pi+1,pj+1) pvY(pi+1,pj+1);
	     pvX(pi+1,pj) pvY(pi+1,pj);
	     pvX(pi,pj) pvY(pi,pj)];
	% for pp=1:size(P,1)
	%   plot(P(pp,1),P(pp,2),'o');
	%   pause
	% end
	  [pint,kk,intarea(pi,pj)]=convex_intersect(P,Q);
	% fill(P(:,1),P(:,2),'y');
	% if intarea(pi,pj)>0; fill(pint(:,1),pint(:,2),'g'); end;
	% pause
	end
      end 
      txtkernel{gg(j),i,1}=sub2ind([nr nc],pvY(1:end-1,1:end-1)+0.5,pvX(1:end-1,1:end-1)+0.5);
      txtkernel{gg(j),i,2}=intarea; 
    end
  end
end

% Open file for storing final texture values
if isempty(tii)
  txtfname=sprintf('%s.%s',char(txtyourfname(1)),char(txtyourfname(2)));
  tii=[1:size(data1,3)];
else
  txtfname=sprintf('%s_%d-%d.%s',char(txtyourfname(1)),tii(1),tii(end),char(txtyourfname(2)));
end
fid=fopen(txtfname,'w','b');
% write the header
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',size(centroids,1));
fprintf(fid,sprintf('delta-t=%2.6f\n',1/fps));
fprintf(fid,'n_samples=%d\n',length(tii));
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=%s\n',date);
fprintf(fid,'data_source=%s,%s,%s\n',study,runlabel,scanlabel);
dat_precision=sprintf('uint%d',8*4);

disp('Projection and mapping ...');
mapthresh=25;  % F value   THIS VALUE IS CURRENTLY NOT USED, MWK 5/31/04
top_thresh=25; % angle, degrees from z axis
for k=1:length(tii)
  ti=tii(k);
  dat=zeros(nr,nc,4);
  dat(:,:,1)=data1(:,:,ti);
  dat(:,:,2)=data2(:,:,ti);
  dat(:,:,3)=data3(:,:,ti);
  dat(:,:,4)=data4(:,:,ti);
  txtyour=zeros(size(norms,1),nviews+1).*NaN;

  %imagesc(dat(:,:,1)); hold on;
  
  for i=1:nviews
    disp(sprintf('Computing candidate DIRECT and EDGE textures for view %d at ti=%d...',i,ti));
    gg=find(viewi(:,i));
    thisdat=dat(:,:,i);
    if ~isempty(gg);
      for j=1:size(gg,1)
        txtyour(gg(j),i)=sum(sum(thisdat(txtkernel{gg(j),i,1}).*txtkernel{gg(j),i,2}))/sum(sum(txtkernel{gg(j),i,2}));
      end
    end
  end
%%%% added fluorescence averaging within surface mesh cells. 6/4/04, MWKay

  % Assign one texture for each cell
  disp(sprintf('Assigning textures for ti=%d...',ti));
  for i=1:size(txtyour,1)
    if (viewi(i,nviews+1)~=0 & viewi(i,nviews+1)<=nviews)   % DIRECT
      txtyour(i,end)=txtyour(i,viewi(i,nviews+1));
    elseif viewi(i,nviews+1)>nviews                         % OVERLAP or EDGE
      gg=find(viewi(i,1:nviews));
      txtyour(i,end)=sum(txtyour(i,gg).*ledge(i,gg))/sum(ledge(i,gg));    % Weighted average
    end  
  end
  
  %--------------------------------------------------
  % Get rid of sites on the top 
  % (this is not needed for decimated models
  if nuketopdata
    thetaz=acos(norms(:,3)).*180/pi;
    thetaz=real(thetaz);
    gg=find(thetaz<=top_thresh);
    txtyour(gg,end)=NaN;
  end
  
  %--------------------------------------------------
  % set low sites to 1
  %gg=find(txtyour(:,end)<=mapthresh);
  %txtyour(gg,end)=1;
  
  %--------------------------------------------------
  % set NaN sites to 0
  gg=find(isnan(txtyour(:,end)));
  txtyour(gg,end)=0;
  
  %---------------------------------------------
  % save texture
  
  %fwrite(fid,txtyour(:,end),'float');
  % remember to write null data to electrode 1
  fwrite(fid,[0 txtyour(:,end)'].*1000,dat_precision);
  disp(sprintf('Saved texture in %s',txtfname));

end
fclose(fid);

%--------------------------------------------
% identify the edge cells
disp('Finding points along edges ...');
edges=int8(zeros(size(viewi,1),1));
for j=1:size(edges,1)
  if (viewi(neighs(j,1),nviews+1)-viewi(neighs(j,2),nviews+1))~=0
    edges(j)=1;
  elseif (viewi(neighs(j,1),nviews+1)-viewi(neighs(j,3),nviews+1))~=0
    edges(j)=1;
  elseif (viewi(neighs(j,2),nviews+1)-viewi(neighs(j,3),nviews+1))~=0
    edges(j)=1;
  end
end

% find points along the edges (find the indicies)
edgesi=find(edges);
edgepts2col=zeros(length(edgesi),2);
for i=1:length(edgesi)
  if (viewi(edgesi(i),nviews+1)-viewi(neighs(edgesi(i),1),nviews+1))~=0
    edgepts2col(i,:)=cells(edgesi(i),[1 2]);
  elseif (viewi(edgesi(i),nviews+1)-viewi(neighs(edgesi(i),2),nviews+1))~=0
    edgepts2col(i,:)=cells(edgesi(i),[2 3]);
  elseif (viewi(edgesi(i),nviews+1)-viewi(neighs(edgesi(i),3),nviews+1))~=0
    edgepts2col(i,:)=cells(edgesi(i),[3 1]);
  end
end
edgepts=[edgepts2col(:,1); edgepts2col(:,2)];
edgepts=unique(edgepts);   % these are indicies into the pts array!

%---------------------------------------------
% save indicies of points on the edge. THESE ARE POINT INDICIES, NOT CELL INDICIES!
edgsfname=sprintf('%s.%s',char(edgesfname(1)),char(edgesfname(2)));
fid=fopen(edgsfname,'w');
fprintf(fid,'3\n');
fprintf(fid,'Edge point indicies, created by pan_textureF10.m\n');
fprintf(fid,'int\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(edgepts,1));
fprintf(fid,'%d\n',size(edgepts,2));
fwrite(fid,edgepts-1,'int');
fclose(fid);
disp(sprintf('Saved edge point indicies in %s',edgsfname));

%---------------------------------------------
% save view for each cell
vwsfname=sprintf('%s.%s',char(viewsfname(1)),char(viewsfname(2)));
fid=fopen(vwsfname,'w');
fprintf(fid,'3\n');
fprintf(fid,'View data, created by pan_textureF10.m\n');
fprintf(fid,'float\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'%d\n',size(viewi(:,end),1));
fprintf(fid,'%d\n',size(viewi(:,end),2));
fwrite(fid,viewi(:,nviews+1),'float');
fclose(fid);
disp(sprintf('Saved view numbers in %s',vwsfname));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output some useful statistics

npixels=zeros(nviews,1);
viewarea=zeros(nviews+2,1);      % mm
meancellarea=zeros(nviews+2,1);  % mm
stdcellarea=zeros(nviews+2,1);   % mm
averes=zeros(nviews,1);          % mm
meanpixelusedpercell=zeros(nviews,1);
stdpixelusedpercell=zeros(nviews,1);
meanpixelpercell=zeros(nviews,1);
stdpixelpercell=zeros(nviews,1);

% compute total surface area
V1=zeros(size(pts)); V2=V1; V3=V1;
V1=[pts(cells(:,1),1) pts(cells(:,1),2) pts(cells(:,1),3)];
V2=[pts(cells(:,2),1) pts(cells(:,2),2) pts(cells(:,2),3)];
V3=[pts(cells(:,3),1) pts(cells(:,3),2) pts(cells(:,3),3)];
totalSA=sum(cellSA3d(V1,V2,V3));
clear V1;
clear V2;
clear V3;

% compute surface areas of each view
% viewi(:,1-nviews)=0 or 1, indicating whether the cell was mapped by a particular view
% viewi(:,nviews+1)=integer b/w 1 and nviews+2. This number indicates the assigned view. 
%      nviews+1 is overlap, nviews+2 is edge
for i=1:nviews+2
  gg=find(viewi(:,nviews+1)==i);
  if ~isempty(gg);
    V1=zeros(size(gg,1),3); V2=V1; V3=V1;
    % Get verticies of each cell whose centroid lies within the view
    V1=[pts(cells(gg,1),1) pts(cells(gg,1),2) pts(cells(gg,1),3)];
    V2=[pts(cells(gg,2),1) pts(cells(gg,2),2) pts(cells(gg,2),3)];
    V3=[pts(cells(gg,3),1) pts(cells(gg,3),2) pts(cells(gg,3),3)];
    cellarea=cellSA3d(V1,V2,V3);
    viewarea(i)=sum(cellarea);
    meancellarea(i)=mean(cellarea);
    stdcellarea(i)=std(cellarea);
  end
end

% count pixels in each direct view but not within the edges
for i=1:nviews
  gg=find(viewi(:,nviews+1)==i);
  if ~isempty(gg);
    [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
    Xigg=Xigg-pix(i,1)+1; Yigg=Yigg-pix(i,4)+1;
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
    directmask=zeros(nr,nc);
    directmask(sub2ind([nr,nc],round(Yigg),round(Xigg)))=1;  
    directmask=bwmorph(directmask,'close',Inf);  
    npixels(i)=length(find(directmask));
    averes(i)=sqrt(viewarea(i)/npixels(i));
  end
end

% count pixels per cell for each direct view
for i=1:nviews
  gg=find(viewi(:,nviews+1)==i);
  pixpercell=zeros(size(gg));
  for j=1:length(gg)
    pixusedpercell(j)=size(txtkernel{gg(j),i,1},1)*size(txtkernel{gg(j),i,1},2);
    pixpercell(j)=sum(sum(txtkernel{gg(j),i,2}));
  end
  meanpixelusedpercell(i)=mean(pixusedpercell);
  stdpixelusedpercell(i)=std(pixusedpercell);
  meanpixelpercell(i)=mean(pixpercell);
  stdpixelpercell(i)=std(pixpercell);  
end

%---------------------------------------------------------
% Write useful information to a readable textfile
infofid=fopen(infofname,'w');
fprintf(infofid,'%s, %s, %s\n',study,scanlabel,runlabel);
fprintf(infofid,'\n');
fprintf(infofid,'Surface mesh parameters:\n');
fprintf(infofid,'\n');
fprintf(infofid,'1) number of cells: %d\n',size(cells,1));
fprintf(infofid,'2) number of verticies: %d\n',size(pts,1));
fprintf(infofid,'3) total surface area: %4.3f cm\n',totalSA/100);
fprintf(infofid,'4) surface area mapped: %4.3f cm\n',sum(viewarea)/100);
fprintf(infofid,'5) surface area not mapped: %4.3f cm\n',totalSA/100-sum(viewarea)/100);

fprintf(infofid,'\n');
fprintf(infofid,'Mapped surface areas for each view:\n');
for i=1:nviews+2
  fprintf(infofid,'\n');
  if i<=nviews
    fprintf(infofid,'Direct view %d\n',i);
  elseif i==nviews+1
    fprintf(infofid,'Overlap view (number %d)\n',i);
  else 
    fprintf(infofid,'Edge view (number %d)\n',i);
  end
  fprintf(infofid,'1) total surface area mapped: %4.3f cm\n',viewarea(i)/100);
  fprintf(infofid,'2) mean surface area mapped by each cell: %4.3f mm\n',meancellarea(i));
  fprintf(infofid,'3) std surface area mapped by each cell: %4.3f mm\n',stdcellarea(i));
end

fprintf(infofid,'\n');
fprintf(infofid,'Additional parameters for each direct view:\n');
for i=1:nviews
  fprintf(infofid,'\n');
  fprintf(infofid,'Direct view %d\n',i);
  fprintf(infofid,'1) number of pixels in the view: %d\n',npixels(i));
  fprintf(infofid,'2) effective pixel resolution (computed as the sqrt of the total\n');
  fprintf(infofid,'   surface area mapped divided by number of pixels in the view: %4.3f mm\n',averes(i));
  fprintf(infofid,'3) mean number of pixels contributing to each cell: %4.3f\n',meanpixelusedpercell(i));
  fprintf(infofid,'4) std number of pixels contributing to each cell: %4.3f\n',stdpixelusedpercell(i));
  fprintf(infofid,'5) mean total number of pixels covered by each cell: %4.3f\n',meanpixelpercell(i));
  fprintf(infofid,'6) std total number of pixels covered by each cell: %4.3f\n',stdpixelpercell(i));
end  
fclose(infofid);

