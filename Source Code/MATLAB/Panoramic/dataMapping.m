clear all

%% - added fluorescence averaging within surface mesh cells. 6/4/04, MWKay
%% - removed signal normalization for subsequent signal quality measurement,
%%   filtering and normalization in pan_postprocess.m. 11/09/04, MWKay
%% - added setup_pan_texture.m. 11/15/04, MWKay
%% - added output of texturing information. 11/23/04, MWKay
%% - added pixel fine-tuning of mapping data positions on geometry. 12/03/04, MWKay
%% - corrected nonrobust neighbors manipulation. 2/23/04 MWKay

fprintf('\npan_textureF11.m \nVersion: 2/23/05 \n\n');

disp('Press Return to read setup_pan_texture.m');
disp('(located in current dir) and begin processing ...'); pause;

% Read user defined parameters and filenames
% setup_pan_textureMod_0527
setup_pan_textureMod

%-----------------------------------
% Load mapping camera calibration stuff
% cd([datadir studydir])
cd([datadir studydir])
nviews=size(mapcalfnames,1);

Rcmap=NaN.*ones(4,4,nviews);
Parmap=NaN.*ones(8,nviews);
Posmap=NaN.*ones(6,nviews);
xyzmap=NaN.*ones(nviews,4);
for i=1:nviews
    % Grab camera parameters for each CMOS camera
    pospar=load(mapcalfnames{i});
    % Compute variables need for rotation matrix
    wa=pospar(4)*pi/180;  % omega, Wx rotation (world basis, x axis)
    pa=pospar(5)*pi/180;  % psi, y rotation  (world basis, y axis)
    ra=pospar(6)*pi/180;  % kappa, z rotation (world basis, z axis)
    cw=cos(wa); sw=sin(wa);
    cp=cos(pa); sp=sin(pa);
    cr=cos(ra); sr=sin(ra);
    R=zeros(3,3);
    
    % Populate R with rotation matrix from the CMOS cameras
    R(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]';
    R(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
    R(:,3)=[-sp cp*sw cp*cw]';
    % Combine R with translation matrix to create tranformation matrix
    Rcmap(:,:,i)=[[R'; [0 0 0]] [pospar(1) pospar(2) pospar(3) 1]'];
    
    % Compute xyz positions of orbit from calibration matricies
    xyzmap(i,:)=(Rcmap(:,:,i)\[0 0 0 1]')';
    Parmap(:,i)=pospar(7:14);
    Posmap(:,i)=pospar(1:6);
end

% Visualize camera positions
scatter3(0,0,0,'r*')
text(0,0,0,'  Origin','FontSize',18)
hold on
view([0 90])
for n = 1:4
    if n == 5
        scatter3(xyzmap(n,1),xyzmap(n,2),xyzmap(n,3),'go')
        cmd =  sprintf('text(xyzmap(n,1),xyzmap(n,2),xyzmap(n,3),''  Geometry'',''Color'',''k'',''FontSize'',18)');
    else
        scatter3(xyzmap(n,1),xyzmap(n,2),xyzmap(n,3),'bo')
        cmd = sprintf('text(xyzmap(n,1),xyzmap(n,2),xyzmap(n,3),''  Mapping %s'',''Color'',''k'',''FontSize'',18)',scanlabel(n));
    end
    eval(cmd)
end
set(gca,'FontSize',14)
hold off

% Grab camera view directions
lookmap=ones(size(Rcmap,3),4);
looknmap=ones(size(Rcmap,3),3);
for i=1:size(Rcmap,3)
    lookmap(i,:)=(Rcmap(1:4,1:4,i)\[0 0 1 1]')';       % camera view direction
    lookmap(i,1:3)=lookmap(i,1:3)-xyzmap(i,1:3);            % camera view direction
    looknmap(i,:)=lookmap(i,1:3)./norm(lookmap(i,1:3));     % normalized view direction
end

%-------------------------------------------
% load geometry data
disp('Loading geometry data ...');
[centroids,~,~]=readdat(char(geomfnames(1)));
[norms,~,~]=readdat(char(geomfnames(4)));
[neighnum,neighs,~,~]=readneighs(char(geomfnames(5)));   % 2/23/05 MWKay
% neighs is a structure
[cells,~,~]=readdat(char(geomfnames(3)));
[pts,numparams,txtparams]=readdat(char(geomfnames(2)));
cells=cells+ones(size(cells));   % VTK starts with zero, matlab starts with one
for i=1:size(cells,1)            % VTK starts with zero, matlab starts with one
    for j=1:neighnum(i)
        neighs{i}(j)=neighs{i}(j)+1;
    end
end

%---------------------------------------
% load mapping data
disp('Loading the mapping data ...');
bgMaskCams = zeros(100,100,4);
figure
for n = 1:4
    % Load data from CMOS cameras
    load(datafnames{n})
    if n == 1
        bgMaskCams = zeros(size(bgMask,1),size(bgMask,2),4);
        cmosDataAll = cell(n,1);
        cmosDataRawAll = cmosDataAll;
        analogAll = zeros(n,size(analog,2));
    end
    % Place data in consolidated variables
    bgMaskCams(:,:,n) = bgMask;
    cmosDataAll{n} = cmosData;
%     cmosDataAll{n} = texture;
    cmosDataRawAll{n} = cmosRawData;
    analogAll(n,:) = analog;
    % Visualize data foreground binaries
    subplot(2,2,n)
    imagesc(bgMaskCams(:,:,n))
end

nr=size(bgMaskCams,1);
nc=size(bgMaskCams,2);

%--------------------------------------------
% Calculate geometry masks (should be the same resolution as mapping camera)
X = zeros(size(pts,1),4);
Xshift = cell(4,1);
Y = zeros(size(pts,1),4);
Yshift = cell(4,1);
shift = zeros(4,2);
geommasks = zeros(100,100,4);
camMap = figure;
geoMap = figure;
for n = 1:4
    % Projection of points onto mapping cameras
    [X(:,n),Y(:,n)] = pred(pts,Parmap(:,n),Posmap(:,n),camera);
    
    % Add user specified shift
    [Yshift{n},Xshift{n},shift(n,:)] = usershiftindMod(bgMaskCams(:,:,n),Y(:,n),X(:,n));
    
    figure(camMap)
    subplot(2,2,n)
    scatter(Xshift{n},Yshift{n},'bo')
    set(gca,'XLim',[1 100],'YLim',[1 100])
    
    % Round values and remove those outside the boundary
    x = round(Xshift{n}); y = round(Yshift{n});
    % Identify X values to remove
    rmX = (x > nc) + (x < 1);
    if sum(rmX) > 0
        rmX = rmX.*(1:length(rmX))';
        rmX = unique(rmX);
        rmX = rmX(2:end);
    else
        rmX = [];
    end
    % Identify Y values to remove
    rmY = (y > nr) + (y < 1);
    if sum(rmY) > 0
        rmY = rmY.*(1:length(rmY))';
        rmY = unique(rmY);
        rmY = rmY(2:end);
    else
        rmY = [];
    end
    % Combine points removed based on X and Y
    rm = [rmX;rmY];
    x(rm) = [];
    y(rm) = [];
    % Identify heart pixels to keep
    ind = sub2ind([100 100],x,y);
    ind = unique(ind);
    geommasks((n-1)*100^2+ind) = 1;
    geommasks(:,:,n) = imfill(geommasks(:,:,n))';
%     geommasks(:,:,n) = flip(geommasks(:,:,n),2);
    figure(geoMap)
    subplot(2,2,n)
    imagesc(geommasks(:,:,n))
    
% % %     [ri,ci] = find(geommasks(:,:,n));
% % %     mask(sub2ind([size(mask,1) size(mask,2)],newri,newci)+...
% % %         (size(mask,1)*size(mask,2)*(n-1))) = 1;
end

% % silhsData = load([datadir studydir '/silhs1.mat']);
% % geommask1 = silhsData.silhs(:,:,1);
% % geommask2 = silhsData.silhs(:,:,18);�
% % geommask3 = silhsData.silhs(:,:,36);
% % geommask4 = silhsData.silhs(:,:,54);



% Compare camera and geometry masks
figure
mask = bgMaskCams.*geommasks;
for n = 1:4
    subplot(2,2,n)
    imagesc(mask(:,:,n))
end


%---------------------------------------
% % % % Position and crop the masks according to user input
% % % mask = zeros(size(bgMaskCams,1),size(bgMaskCams,2),4);
% % % for n = 1:4
% % %     fprintf('Position mask for view %d\n',n)
% % %     [ri,ci] = find(geommasks(:,:,n));
% % %     [newri,newci] = usershiftindMod(bgMaskCams(:,:,n),ri,ci);
% % %     mask(sub2ind([size(mask,1) size(mask,2)],newri,newci)+...
% % %         (size(mask,1)*size(mask,2)*(n-1))) = 1;
% % % end

% % %---------------------------------------
% % % Detrend the data
% % disp('detrending dataset 1 ...');
% % [gr,gc]=find(mask1);
% % data1=op_detrend(data1,gr,gc);
% % disp('detrending dataset 2 ...');
% % [gr,gc]=find(mask2);
% % data2=op_detrend(data2,gr,gc);
% % disp('detrending dataset 3 ...');
% % [gr,gc]=find(mask3);
% % data3=op_detrend(data3,gr,gc);
% % disp('detrending dataset 4 ...');
% % [gr,gc]=find(mask4);
% % data4=op_detrend(data4,gr,gc);

%-------------------------------------------
% compute angles b/w normals and view directions
disp('Computing angles between normals and view directions .... ');
langle=zeros(size(norms,1),nviews+1);
fprintf('%d normals and %d views',size(langle,1),nviews);
two_comp=0;   % 1 for 2 components, 0 for three
if ~two_comp
    disp('Computing angles using x,y, and z components.');
    gg=norms*looknmap';
    gg=acos(gg);
    gg=real(gg);
    langle(:,1:nviews)=gg.*180/pi;  % x,y, and z components
elseif two_comp
    disp('Computing angles using x and y components.')
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
    fprintf('Computing weights for view %d ...',i);
    gg=find(langle(:,i)>90);  % find the front face
%     [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcmap(:,1:4,i),Parmap(:,i),camera);
    [Xigg,Yigg]=pred(centroids(gg,:),Parmap(:,i),Posmap(:,i),camera);
% %     Xigg=Xigg-pix(i,1)+1; Yigg=Yigg-pix(i,4)+1;
%     Xigg=Xigg-size(mask,2)+1; Yigg=Yigg-size(mask,2)+1; % adjusted for 100x100 cams
    Xigg = Xigg+shift(i,2);
    Yigg = Yigg+shift(i,1);
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
    com=sprintf('[conesilh,kernel]=conefilt(mask(:,:,%d),edgefiltval,1);',i);
    eval(com);
    ledge(gg,i)=interp2(conesilh,Xigg,Yigg);
    clear conesilh;
end

%-------------------------------------------
disp('Finding max edge weights ....');
for i=1:size(ledge,1)
    gg=find(ledge(i,1:nviews)==max(ledge(i,1:nviews)));
    if (~isempty(gg) && length(gg)==1)
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
    fprintf('Finding cells in DIRECT view %d ...',i);
    gg=find(langle(:,i)>=minangleDIRECT & ledge(:,i)>=minedgeDIRECT);
    if isempty(gg)
        fprintf('Found NO cells in view %d ...',i);
    else
        viewi(gg,i)=1;
        fprintf('Found %d cells in DIRECT view %d ...',length(gg),i);
        fprintf('Smoothing edges of DIRECT view %d...',i);
        [viewi(:,i)]=smoothregionbdr(viewi(:,i),neighnum,neighs);
    end
end
clear pvect;

fprintf('Assigning OVERLAP and EDGE views...');
for i=1:size(viewi,1)
    gg=find(viewi(i,1:nviews));
    if ~isempty(gg);
        if length(gg)==1
            viewi(i,nviews+1)=gg;          % DIRECT
        else
            viewi(i,nviews+1)=nviews+1;    % OVERLAP
        end
    else
        gg=find(langle(i,1:nviews)>=90 & ledge(i,1:nviews)>=minedgeEDGE);
        if ~isempty(gg)
            viewi(i,nviews+1)=nviews+2;    % EDGE
            viewi(i,gg)=1;  % average these views
        end
    end
end

disp('Finding indicies of cells with multiple views...');
mults=sum(viewi(:,1:nviews),2);
multsover1=find(mults>1);
multsii=zeros(length(find(viewi(multsover1,1:nviews))),1);
multsjj=multsii;
kount=0;
for i=1:size(viewi,1)
    gg=find(viewi(i,1:nviews));
    if ~isempty(gg) && length(gg)>1
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
    fprintf('Doing this for view %d ...',i);
    gg=find(viewi(:,i));
    if ~isempty(gg);
        [Xigg,Yigg]=pred(centroids(gg,:),Parmap(:,i),Posmap(:,i),camera);
%         Xigg=Xigg-pix(i,1)+1; Yigg=Yigg-pix(i,4)+1;
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
        V1=[pts(cells(gg,1),1) pts(cells(gg,1),2) pts(cells(gg,1),3)];
        V2=[pts(cells(gg,2),1) pts(cells(gg,2),2) pts(cells(gg,2),3)];
        V3=[pts(cells(gg,3),1) pts(cells(gg,3),2) pts(cells(gg,3),3)];
%         [XV1gg,YV1gg]=pred2([V1(:,1) V1(:,2) V1(:,3) ones(size(V1,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
        [XV1gg,YV1gg]=pred(V1,Parmap(:,i),Posmap(:,i),camera);
%         XV1gg=XV1gg-pix(i,1)+1;YV1gg=YV1gg-pix(i,4)+1;
%         [XV2gg,YV2gg]=pred2([V2(:,1) V2(:,2) V2(:,3) ones(size(V2,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
        [XV2gg,YV2gg]=pred(V2,Parmap(:,i),Posmap(:,i),camera);
%         XV2gg=XV2gg-pix(i,1)+1;YV2gg=YV2gg-pix(i,4)+1;
%         [XV3gg,YV3gg]=pred2([V3(:,1) V3(:,2) V3(:,3) ones(size(V3,1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
        [XV3gg,YV3gg]=pred(V3,Parmap(1:8,i),Posmap(:,i),camera);
%         XV3gg=XV3gg-pix(i,1)+1;YV3gg=YV3gg-pix(i,4)+1;
        for j=1:size(gg,1)
            % Define the triangle's bounding square and round up to pixel edges
            % To do this, add or subtract 0.5 to grid the integer-valued pixels
            VXmax=ceil(max([XV1gg(j) XV2gg(j) XV3gg(j)])-0.5)+0.5;
            VXmin=floor(min([XV1gg(j) XV2gg(j) XV3gg(j)])+0.5)-0.5;
            VYmax=ceil(max([YV1gg(j) YV2gg(j) YV3gg(j)])-0.5)+0.5;
            VYmin=floor(min([YV1gg(j) YV2gg(j) YV3gg(j)])+0.5)-0.5;
            if VXmax>nc+0.5
                VXmax=nc+0.5;
            end
            if VXmin<0.5
                VXmin=0.5;
            end
            if VYmax>nr+0.5
                VYmax=nr+0.5;
            end;
            if VYmin<0.5
                VYmin=0.5;
            end
            % Define the pixel grid
            [pvX,pvY]=meshgrid((VXmin:1:VXmax),(VYmin:1:VYmax));
            intarea=zeros(size(pvX,1)-1,size(pvX,2)-1);
            Q=[XV1gg(j),YV1gg(j);XV2gg(j),YV2gg(j);XV3gg(j),YV3gg(j)];
            Q=checkpoly(Q*100)/100;  % multiply and divide by 100 to avoid numerical imprecision with checkpoly.m
            txtyour(gg(j),i)=0;
            % plot(Q(:,1),Q(:,2),'k*');
            for pi=1:size(pvX,1)-1
                for pj=1:size(pvY,2)-1
                    % assign P for cw rotation (yes, cw! it works...)
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
            txtkernel{gg(j),i,1}=sub2ind([nr nc],pvY(1:end-1,...
                1:end-1)+0.5,pvX(1:end-1,1:end-1)+0.5);
            txtkernel{gg(j),i,2}=intarea;
        end
    end
end

% Open file for storing final texture values
tii = [];
if isempty(tii)
    txtfname=sprintf('%s.%s',char(txtyourfname(1)),char(txtyourfname(2)));
    tii=(1:size(cmosData,3));
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
fprintf(fid,'data_source=%s%s/%s%s\n',datadir,studydir,runnumber,runlabel);
dat_precision=sprintf('uint%d',8*4);

disp('Projection and mapping ...');
mapthresh=25;  % F value   THIS VALUE IS CURRENTLY NOT USED, MWK 5/31/04
top_thresh=25; % angle, degrees from z axis
h = waitbar(0,'Projection and mapping...');
dataProj = zeros(size(txtyour,1),length(tii));
% for k=1:length(tii)
for k = 1:length(tii)
    dat=zeros(nr,nc,4);
    dat(:,:,1)=cmosDataAll{1}(:,:,k);
    dat(:,:,2)=cmosDataAll{2}(:,:,k);
    dat(:,:,3)=cmosDataAll{3}(:,:,k);
    dat(:,:,4)=cmosDataAll{4}(:,:,k);
    txtyour=nan(size(norms,1),nviews+1);
    
    for i=1:nviews
        fprintf('Computing candidate DIRECT and EDGE textures for view %d at ti=%d\n',i,k);
        % find the faces visible from this view
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
    fprintf('Assigning textures for k=%d\n',k);
    for i=1:size(txtyour,1)
        if (viewi(i,nviews+1)~=0 && viewi(i,nviews+1)<=nviews)   % DIRECT
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
%     fwrite(fid,[0 txtyour(:,end)'].*1000,dat_precision);
%     fprintf('Saved texture in %s',txtfname);
    dataProj(:,k) = txtyour(:,end);
    
    waitbar(k/length(tii))
end
close(h)
fclose(fid);

%--------------------------------------------
% identify the edge cells
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
    end
    if neighnum(j)==2  % if cell is diff than either of its 2 neighs then on edge
        if (viewi(neighs{j}(1),end)-viewi(j,end))~=0
            edges(j)=1;
        elseif (viewi(j,end)-viewi(neighs{j}(2),end))~=0
            edges(j)=1;
        end
    end
end

%%%% THIS IS BROKEN - JUST SKIP IT FOR NOW 10/27/2006

%- % find points along the edges (find the indicies)
%- edgesi=find(edges);
%- edgepts2col=zeros(length(edgesi),2);
%- for i=1:length(edgesi)
%-   % lazy approach:
%-   % simply check all the edges of each previously identified cell
%-   % then 'unique' the list
%-   j=edgesi(i);
%-   if neighnum(j)==3
%-     if (viewi(j,end)-viewi(neighs{j}(1),end))~=0
%-       edgepts2col(i,:)=cells(j,[1 2]);
%-     elseif (viewi(j,end)-viewi(neighs{j}(2),end))~=0
%-       edgepts2col(i,:)=cells(j,[2 3]);
%-     elseif (viewi(j,end)-viewi(neighs{j}(3),end))~=0
%-       edgepts2col(i,:)=cells(j,[3 1]);
%-     end
%-   elseif neighnum(j)==2 % I do not know which edge these neighs correspond to
%-     if (viewi(j,end)-viewi(neighs{j}(1),end))~=0
%-       % then find the points in common
%-       edgepts2col(i,:)=intersect(cells(j,:),cells(neighs{j}(1)));
%-     elseif (viewi(j,end)-viewi(neighs{j}(2),end))~=0
%-       % then find the points in common
%-       edgepts2col(i,:)=intersect(cells(j,:),cells(neighs{j}(2)));
%-     end
%-   elseif neighnum(j)==1  % I do not know which edge this neigh corresponds to
%-      if (viewi(j,end)-viewi(neighs{j}(1),end))~=0
%-       % then find the points in common
%-       edgepts2col(i,:)=intersect(cells(j,:),cells(neighs{j}(1)));
%-     end
%-   end
%- end
%- edgepts=[edgepts2col(:,1); edgepts2col(:,2)];
%- edgepts=unique(edgepts);   % these are indicies into the pts array!

%---------------------------------------------
%- % save indicies of points on the edge. THESE ARE POINT INDICIES, NOT CELL INDICIES!
%- edgsfname=sprintf('%s.%s',char(edgesfname(1)),char(edgesfname(2)));
%- fid=fopen(edgsfname,'w');
%- fprintf(fid,'3\n');
%- fprintf(fid,'Edge point indicies, created by pan_textureF10.m\n');
%- fprintf(fid,'int\n');
%- fprintf(fid,'%s\n',datestr(now));
%- fprintf(fid,'%d\n',size(edgepts,1));
%- fprintf(fid,'%d\n',size(edgepts,2));
%- fwrite(fid,edgepts-1,'int');
%- fclose(fid);
%- disp(sprintf('Saved edge point indicies in %s',edgsfname));

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

% compute surface areas of each vpiew
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
%         [Xigg,Yigg]=pred2([centroids(gg,1) centroids(gg,2) centroids(gg,3) ones(size(centroids(gg,1),1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
        [Xigg,Yigg]=pred(centroids(gg,:),Parmap(:,i),Posmap(:,i),camera);
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

