function [status,fnametag] = cubeCarvingMod(geodir,silh,lims,dtheta,n_images,dofrontback,r,...
    rr,irsort,rsort,startlevel,maxlevel,msgCenter)

% Description: A modification of the occluding contours portion of whirl.m
% written by Dr. Matthew Kay. The silhsProcess_callback was getting
% excessively long so I decided to consolidate this code into a function.
% ATTENION: You will need to modify lines 129, 139, 149, 159, 169 to point
% to the Panoramic directory which should be downloaded with this set of
% functions. The recommended location for the Panoramic directory is in the
% Matlab root directory.
%
% Inputs:
% hdir = the home directory
% silh = silhouettes
% lims = limits of the bounding box for the silhouettes
% dtheta = degree of each step
% n_images = total number of images used
% dofrontback = (1) collapse redudant images, (0) do not collapse
% r = image angles
% rr = list of angles being used and their index
% irsort =
% inumsort =
% rsort =
% startlevel = 
% octMax
% msgCenter
%
% Outputs:
%
% Author: Christopher Gloschat
% Date: June 23, 2016
%
% Modification Log:
%
%
%% Code %%
% Load external camera calibration values
[filename,pathname] = uigetfile('*.mat','Select camera calibration file (e.g. cal_009.mat).');
if filename ~= 0
    status = 1;
    load([pathname '/' filename]);
% % %     disp('Camera parameters loaded.')
    
    
    %% Load octree basis %%
    [octb,octbi] = octree_basis;
    
    % Establish size and origin of initial carving cube
    % Keep it cubic!!!  This makes scaling volume and surface area much easier!
    % % xdim=6.0*25.4;  % mm
    xdim=6.0*20;  % mm
    ydim=xdim;    % mm
    zdim=xdim;    % mm
    
    X0=-0.5*xdim;
    Xn=0.5*xdim;
    Y0=-0.5*ydim;
    Yn=0.5*ydim;
    % % Z0=0.45*zdim;       GLO-EDIT
    % % Zn=-0.55*zdim;      GLO-EDIT
    Z0=-0.5*zdim;
    Zn=0.5*zdim;
    
    levs=[(1:10)' (Xn-X0)./(2.^(1:10)')];        % Level number, delt
%     disp([(1:10)' (Xn-X0)./(2.^(1:10)')])
    set(msgCenter,'String',sprintf('Octree levels:\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f\n%d\t%0.2f',...
        levs(1,1),levs(1,2),levs(2,1),levs(2,2),levs(3,1),levs(3,2),...
        levs(4,1),levs(4,2),levs(5,1),levs(5,2),levs(6,1),levs(6,2),...
        levs(7,1),levs(7,2),levs(8,1),levs(8,2),levs(9,1),levs(9,2),...
        levs(10,1),levs(10,2)),'FontSize',9);
%     startlevel=input('Start at level number: ');
%     maxlevel=input('Max level number: ');
    
    fnametag=inputdlg('Enter filename tag (ie: _povcyl_1mm_corrected):');
    % savedir=input('Save datafiles in this directory, leave off the "/" (ie: povcyl_05mm): ','s');
    if ~isempty(fnametag)
        savedir = geodir;
        
        % verticies of voxel h-1 of octree t is (remember, voxels are from 2:9):
        % vert(octree(t,find(octree(t,:,h)),1),:);
        
        level=startlevel;
        delt(1,1)=(Xn-X0)/(2^level);
        delt(1,2)=(Yn-Y0)/(2^level);
        delt(1,3)=(Zn-Z0)/(2^level);
        
        %% Record analysis parameters in textfile %%
        parafname=sprintf('%s/whirl%s.txt',savedir,fnametag{1});
        parafid=fopen(parafname,'w');
        fprintf(parafid,'dtheta, degree step: %4.3f \n',dtheta);
        fprintf(parafid,'n_images, total number of snapshots: %d \n',n_images);
        fprintf(parafid,'dofrontback, 0: use all silhouettes or 1: use only largest silhouettes: %d \n',dofrontback);
        fprintf(parafid,'xdim, size of clay cube in x dir (mm): %4.3f \n',xdim);
        fprintf(parafid,'(X0,Xn) mm: (%4.3f,%4.3f) \n',X0,Xn);
        fprintf(parafid,'ydim, size of clay cube in y dir (mm): %4.3f \n',ydim);
        fprintf(parafid,'(Y0,Yn) mm: (%4.3f,%4.3f) \n',Y0,Yn);
        fprintf(parafid,'zdim, size of clay cube in z dir (mm): %4.3f \n',zdim);
        fprintf(parafid,'(Z0,Zn) mm: (%4.3f,%4.3f) \n',Z0,Zn);
        fprintf(parafid,'Levels: \n');
        for i=1:size(levs,1)
            for j=1:size(levs,2)
                fprintf(parafid,'%4.6f   ',levs(i,j));
            end
            fprintf(parafid,'\n');
        end
        fprintf(parafid,'startlevel, start at this spacing level: %d \n',startlevel);
        fprintf(parafid,'maxlevel, stop after completion of this spacing level: %d \n',maxlevel);
        fprintf(parafid,'fnametag: %s \n',fnametag{1});
        fprintf(parafid,'savedir: %s \n',savedir);
        % % fprintf(parafid,'Camera calibration, Position file: %s \n',Rcfname);
        fprintf(parafid,'Camera calibration Position and Parameters file: %s \n',[geodir '/Calibration/' filename]);
        fprintf(parafid,'Angle array (r):\n');
        for i=1:length(r)
            fprintf(parafid,'%4.6f \n',r(i));
        end
        fprintf(parafid,'Hemisphere array (rr):\n');
        for i=1:size(rr,1)
            for j=1:size(rr,2)
                fprintf(parafid,'%4.6f   ',rr(i,j));
            end
            fprintf(parafid,'\n');
        end
        fprintf(parafid,'Sorted snapshot array (irsort):\n');
        for i=1:size(irsort,1)
            for j=1:size(irsort,2)
                fprintf(parafid,'%4.6f   ',irsort(i,j));
            end
            fprintf(parafid,'\n');
        end
        
        %%
        goodthresh=0.95;
        % Determine rotation matrix for cube based on directionality of rotation
        if strcmp(cwtxt,'ccw')
            Rz = [cosd(-dtheta) -sind(-dtheta) 0; sind(-dtheta) cosd(-dtheta) 0; 0 0 1];
        else
            Rz = [cosd(dtheta) -sind(dtheta) 0; sind(dtheta) cosd(dtheta) 0; 0 0 1];
        end
        
        vert_out = [];
        vox_rm = [];
        
        count = 0;
        h = waitbar(count/(maxlevel-startlevel+1),sprintf('Creating octree for level %d...',level));
        while level<=maxlevel
            %--------------------------------------------
            % Create the initial verticies matrix (verts) and octree matrix (octree)
            % vert(:,1)=x
            % vert(:,2)=y
            % vert(:,3)=z
            % vert(:,4)=v (0 for outside, 1 for inside, NaN if not checked
            % vert(:,5)=level, resolution level when this vertex was added
            % octree is 3 dimensional: (Num of octrees) x (27 indicies denoting rows in vert)
            % x (9, where octree(:,:,1) contains vert indicies and octree(:,:,2:9) denotes voxel verticies
            % octree(octree_num,:,1)=vert indicies
            % octree(octree_num,:,2:9)=voxel verticies
            % At each level the octree matrix is rebuilt and the vert matrix is expanded.
            
            % Chris' Code %
            % The bounding box itself will not have a vertex inside the   %
            % heart volume, so the first step is to crate and subdivide   %
            % the bounding box.                                           %
% % %             
% % %             % Create the first subdivision from the the bounding box
% % %             if level == startlevel
% % %                 vert = octb*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
% % %                 vert = vert + repmat([X0 Y0 Z0],[size(vert,1) 1]);
% % %                 vert = [vert nan(size(vert,1),1) ones(size(vert,1),1)];
% % %                 % Create the octree variable referencing all voxels
% % %                 octree = octbi ~= 0;
% % %                 octree = octree.*reshape((1:size(octree,1)*size(octree,2)),[size(octree,1) size(octree,2)]);
% % %                 octree = find(octree);
% % %                 octree = reshape(octree,[8 8]);
% % %                 tmp = repmat((1:27)',[1 8]);
% % %                 octree = tmp(octree);
            % Matt's Code %
            if (level==1 && level==startlevel)
                %             fprintf('Creating octree for level 1...\n')
                vert=octb*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(27,3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert=vert+sshift;
                clear sshift
                vert(:,4)=NaN;
                vert(:,5)=1;
                octree=zeros(1,27);
                octree=1:27;
            elseif (level==2 && level==startlevel)
                %             waitbar(1/maxlevel,'Creating octree for level 1...');
                %             fprintf('Loading level 2 octree from datafile...\n');
                load /Users/Chris/Documents/MATLAB/Panoramic/level2
                vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(size(vert,1),3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert(:,1:3)=vert(:,1:3)+sshift;
                clear sshift
            elseif (level==3 && level==startlevel)
                %             fprintf('Loading level 3 octree from datafile...\n');
                load /Users/Chris/Documents/MATLAB/Panoramic/level3
                vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(size(vert,1),3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert(:,1:3)=vert(:,1:3)+sshift;
                clear sshift
            elseif (level==4 && level==startlevel)
                %             fprintf('Loading level 4 octree from datafile...\n');
                load /Users/Chris/Documents/MATLAB/Panoramic/level4
                vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(size(vert,1),3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert(:,1:3)=vert(:,1:3)+sshift;
                clear sshift
            elseif (level==5 && level==startlevel)
                %             fprintf('Loading level 5 octree from datafile...\n');
                load /Users/Chris/Documents/MATLAB/Panoramic/level5
                vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(size(vert,1),3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert(:,1:3)=vert(:,1:3)+sshift;
                clear sshift
            elseif (level==6 && level==startlevel)
                %             fprintf('Loading level 6 octree from datafile...\n');
                load /Users/Chris/Documents/MATLAB/Panoramic/level6
                vert(:,1:3)=vert(:,1:3)*[(Xn-X0) 0 0; 0 (Yn-Y0) 0; 0 0 (Zn-Z0)];
                sshift=zeros(size(vert,1),3);
                sshift(:,1)=X0;
                sshift(:,2)=Y0;
                sshift(:,3)=Z0;
                vert(:,1:3)=vert(:,1:3)+sshift;
                clear sshift
                %--------------------------------------------
            else
                % which voxels were inside AND outside the volume?
                % find them, then subdivide and repeat for each octree
                % to which voxels were those that were inside belong
                % 'good' tells me which verticies were inside the volume
                % and which verticies were outside the volume of the verticies
                % which have not been checked.
                
                %             sprintf('Creating octree for level %d.... ',level)
                count = count+1;
                waitbar(count/(maxlevel-startlevel),h,sprintf('Creating octree for level %d...',level));
                
                good=find(vert(:,4)>=goodthresh);  % Voxels to be subdivided have at least
                % one (but not all) vertex greater than 0.99
                % Use 0.95 to provide good surface definition.
                
                % CHRIS CODE %
                
% % %                 % Find the voxels which contain the vertices in 'good'
% % %                 if length(good) > 1
% % %                     % Reshape to use parallel computing to find voxels
% % %                     goodMat = reshape(good,[1 1 length(good)]);
% % %                 else
% % %                     goodMat = good;
% % %                 end
% % %                 % Replicate to match octree dimensions for logical compare
% % %                 goodMat = repmat(goodMat,[size(octree,1) size(octree,2) 1]);
% % %                 % Create temporary variable of replicated octree
% % %                 tmp = repmat(octree,[1 1 size(goodMat,3)]);
% % %                 % Perform logical comparison
% % %                 tmp = goodMat == tmp;
% % %                 % Sum across all 'good' vertices (3rd dimension)
% % %                 tmp = sum(tmp,3);
                % Sum across all vertices in voxels (columns) to find 
                % voxels of interest, voxels that are both split and
                % completely inside the volume

                
                
                % Voxels that are completely inside the heart volume will
                % have a summed value of 0 while voxels completely inside
                % the volume will have a summed value of 8. Neither of
                % these eligible to be subdivied into a new octree and
                % should be removed from consideration.
% % %                 vox_oi = sum(voxCheck);
% % % % % %                 vox_oi = vox_oi ~= 0;
% % %                 % Check for voxels that are completely outside the volume
% % %                 out = vox_oi == 0;
% % %                 vox_rm = [vox_rm octree(:,out)];
% % %                 octree(:,out) = [];
% % %                 vox_oi(out) = [];
% % %                 % Check for voxels that are completely inside the volume
% % %                 in = vox_oi == 8;
% % %                 vox_rm = [vox_rm octree(:,in)];
% % %                 octree(:,in) = [];
                
% % %                 if sum(in) ~= 0
% % %                     % Keep split voxels
% % %                     vox_oi = find((vox_oi+in)==1);
% % %                 end
                
% % %                 % Find column indices of all other voxels
% % %                 vox_rm = (1:size(octree,2));
% % %                 vox_rm(vox_oi) = [];
% % %                 % Move all other voxels to octree out variable
% % %                 % FIX CONCANTENATION %
% % %                 if ~isempty(vox_rm)
% % %                     oct_out = [oct_out octree(:,vox_rm)];
% % %                     % Remove from octree variable
% % %                     octree(:,vox_rm) = [];
% % %                 end
% % %                 
% % %                 % Subdivide voxels to create new octree variable
% % %                 numOldVoxels = size(octree,2);
% % %                 for n = 1:numOldVoxels
% % %                     % Get vertices associated with voxel
% % %                     newVert = vert(octree(:,n),1:3);
% % %                     % Get dimensions of voxel by subtracting vertices from
% % %                     % opposite corners of voxel
% % %                     tmp = newVert(8,:)-newVert(1,:);
% % %                     tmp = octb*[tmp(1) 0 0; 0 tmp(2) 0; 0 0 tmp(3)];
% % %                     newVert = tmp + repmat(newVert(1,:),[size(tmp,1) 1]);
% % %                     % Check vert variable for duplicates keep track of
% % %                     % index values for newVert old and new
% % %                     tmp = reshape(newVert',[1 size(newVert,2) size(newVert,1)]);
% % %                     tmp = repmat(tmp,[size(vert,1) 1]);
% % %                     newVertInd = repmat(vert(:,1:3),[1 1 size(tmp,3)]);
% % %                     newVertInd = tmp == newVertInd;
% % %                     newVertInd = sum(newVertInd,2) == 3;
% % %                     newVertInd = reshape(newVertInd,[size(newVertInd,1) size(newVertInd,3)]);
% % %                     newVertInd = newVertInd.*repmat((1:size(newVertInd,1))',[1 size(newVertInd,2)]);
% % %                     newVertInd = sum(newVertInd)';
% % %                     % New vertices will have a newVertInd value of 0, add %
% % %                     % these vertices to the vert variable and place their %
% % %                     % new vert index values in the newVertInd variable    %
% % %                     % New vertices will have a 0 value
% % %                     tmp = newVertInd == 0;
% % %                     % A sum not equal to zero means new vertices created
% % %                     if sum(tmp) ~= 0
% % %                         % Add their future index in vert to newVertInd
% % %                         newVertInd(tmp) = (1:sum(tmp))+size(vert,1);
% % %                         % Grab the new vertex values
% % %                         tmp = reshape(newVert(repmat(tmp,[1 3])),[sum(tmp) 3]);
% % %                         % Add the check flag (col 4) and the level the 
% % %                         % vertex was added during (col 5)
% % %                         tmp = [tmp nan(size(tmp,1),1) repmat(level,[size(tmp,1) 1])];
% % %                         % Add new vertices to vert variable
% % %                         vert = [vert;tmp];
% % %                     end
% % %                     % Use the octree template with the new indices to 
% % %                     % create new voxels and add them to octree                  
% % %                     newVertInd = repmat(newVertInd,[1 8]);
% % %                     newVertInd = reshape(newVertInd(octbi ~= 0),[8 8]);
% % %                     octree = [octree newVertInd];
% % %                 end
% % %                 % Remove subdivided voxels from octree
% % %                 vox_rm = [vox_rm octree(:,1:numOldVoxels)];
% % %                 octree(:,1:numOldVoxels) = [];
% % %                 % Reset flags of inside vs outside (col. 4) to 0
% % %                 vert(:,4) = 0;

                
                % Matt Code %
                % Grab the vertex index from good.
                % Find which octrees contain that index
                clear ovi;
                for j=1:size(octree,1)
                    for i=1:length(good)
                        [ovi1,ovi2]=find(octree(j,:)==good(i));
                        if ~isempty(ovi1)
                            if ~exist('ovi')
                                ovi=[j ovi2];
                            else
                                ovi=[ovi; j ovi2];
                            end
                        end
                    end
                end
                
                % ovi- octree and vertex indicies
                % ovi(:,1) contains octree numbers for verticies inside the volume
                % ovi(:,2) contains octree vertex numbers for the verticies
                
                % with which octrees are we dealing?
                octs=unique(ovi(:,1));
                
                % now divide and conquer
                for i=1:length(octs)   % loop through each octree
                    % get vertex numbers for this octree that are inside the volume
                    vs=ovi(ovi(:,1)==octs(i),2);
                    % which voxels have these vertex numbers?
                    % first built a 'mini' octbi in os
                    os=octbi(vs,:);
                    % next find the voxel numbers
                    voxelsi=find(max(os,[],1)>0);  % these are actual voxel numbers for the respective octree!
                    % is a voxel totally inside the volume?
                    vflag=0;
                    for j=1:length(voxelsi)
                        %if length(find(ismember([1 2 3 4 5 6 7 8],os(:,voxelsi(j)))))==8     % CAN I SPEED THIS UP????
                        if length(unique(os(:,voxelsi(j))))==9
                            voxelsi(j)=NaN;
                            vflag=1;   % We have dropped one or more!
                            %sprintf('Dropped one!')
                        end
                    end
                    % drop voxels totally inside the volume
                    if vflag
                        voxelsi=voxelsi(find(~isnan(voxelsi)));
                    end
                    % we now have the voxel number of each voxel of the octree
                    % (stored in voxelsi) that is both inside and outside the volume
                    % now subdivide each voxel, rebuild the 'octree' matrix
                    % and append to the 'vert' matrix
                    
                    for j=1:length(voxelsi)              % loop through each voxel of interest in the current octree, all voxels are the same size
                        vni=find(octbi(:,voxelsi(j))>0);   % find vertex number indicies of the voxel of interest
                        vvi=octree(octs(i),vni);           % actual indicies into vert for each vertex  ***REVISED 1/29/02***
                        vxyz=vert(vvi,1:3);                % vertex locations are in order!
                        
                        % now create an octree inside the voxel
                        % first fill in the verticies that already exist!
                        if i==1 && j==1
                            okount=1;
                            new_octree=zeros(1,27).*NaN;
                            new_octree([1 4 9 11 17 19 24 27])=vvi;
                            delt(level,1)=(vxyz(2,1)-vxyz(1,1))/2;
                            delt(level,2)=(vxyz(5,2)-vxyz(1,2))/2;
                            delt(level,3)=(vxyz(3,3)-vxyz(1,3))/2;
                        else
                            okount=okount+1;
                            new_octree(okount,[1 4 9 11 17 19 24 27])=vvi;  % zeros elsewhere
                        end
                        
                        % New octree limits
                        xn=vxyz(2,1);
                        x0=vxyz(1,1);
                        yn=vxyz(5,2);
                        y0=vxyz(1,2);
                        zn=vxyz(3,3);
                        z0=vxyz(1,3);
                        
                        new_vert=zeros(19,5).*NaN;
                        new_vert(:,5)=level;
                        new_vert(:,1:3)=octb([2 3 5 6 7 8 10 12 13 14 15 16 18 20 21 22 23 25 26],:)*[(xn-x0) 0 0; 0 (yn-y0) 0; 0 0 (zn-z0)];
                        sshift=zeros(19,3);
                        sshift(:,1)=x0;
                        sshift(:,2)=y0;
                        sshift(:,3)=z0;
                        new_vert(:,1:3)=new_vert(:,1:3)+sshift;
                        
                        new_vert_ii=[1:19]+size(vert,1);
                        
                        % do any of these verticies already exist in vert?
                        % if so then delete from new_vert and update new_vert_ii
                        [C,rnew_vert,rvert]=intersect(new_vert(:,1:3),vert(:,1:3),'rows');  % rows of new_vert that are already in vert
                        clear C;
                        
                        if ~isempty(rnew_vert)   % if redundant verticies exist
                            new_vert_ii(rnew_vert)=rvert;  % This makes sense
                            new_new_vert=new_vert;
                            new_new_vert(rnew_vert,:)=[];  % drop redundant rows in new_vert
                            [C,r1,r2]=intersect(new_new_vert(:,1:3),new_vert(:,1:3),'rows');
                            clear C;
                            new_vert_ii(r2)=r1+size(vert,1);
                            new_vert=new_new_vert;
                        end
                        
                        new_octree(okount,[2 3 5 6 7 8 10 12 13 14 15 16 18 20 21 22 23 25 26])=new_vert_ii;
                        vert=[vert; new_vert];
                        
                        if length(unique(new_octree(okount,:)))~=27
                            sprintf('Junk!')
                            return
                        end
                    end  % End voxel loop (j)
                end % End octree loop (i), okount is combination of i and j
                octree=new_octree;
            end
%             --------------------------------------------
            
            % vert(:,4) of NaN: vertex not checked yet
            % vert(:,4) of 0: vertex is definately not in the volume
            % Matt's Code %
            sites2check=find(isnan(vert(:,4))); % indicies into vert
            check_num=length(sites2check);
            % Variable for tracking if vertex is inside the volume
% %             check=ones(check_num,2);       
            check = ones(check_num,73);
            check(:,1)=sites2check;
            % Loop for identifying which vertices are inside the volume
            rot_vert = vert(:,1:3);
            for i=1:length(rsort)
%             for i = 1       
                %for i=1:1
%                 good=find((check(:,2))>0);   % must shave anything greater than zero that is not inside b/c projection
                good=find((check(:,i+1))>0);   % must shave anything greater than zero that is not inside b/c projection
                if isempty(good)             % places verticies that are outside inside, depending upon view!
                    sprintf('No points inside volume. Aborting Snapshots. Restart and increase initial level.')
                    return
                else
                    % Project vertices onto camera model
%                     [Xi,Yi]=pred([vert(check(good,1),1) vert(check(good,1),2) vert(check(good,1),3)],par,pos,camera);
                    [Xi,Yi]=pred([rot_vert(check(good,1),1) rot_vert(check(good,1),2) rot_vert(check(good,1),3)],par,pos,camera);
                    % Which points of Xi and Yi are inside the silhouette?
                    vicinity=ones(length(good),1);
                    vicinity(Xi<lims(i,1) | Xi>lims(i,2) | Yi<lims(i,3) | Yi>lims(i,4))=0;
%                     check(good(~vicinity),2)=0;
                    check(good(~vicinity),i+1)=0;
                    vici=find(vicinity);
                    % Interpolation provides weight to points on the edge
                    Zi=interp2(silh(:,:,i),Xi(vici),Yi(vici));
%                     check(good(vici),2)=Zi;
                    check(good(vici),i+1)=Zi;
                    % Rotate cube
                    if i == 28
                        disp(i)
                    end
                    rot_vert = (Rz*rot_vert')';
% % %                     % Rotate cube
% % %                     tmp = Rz*vert(:,1:3)';
% % %                     vert(:,1:3) = tmp';
                end
            end
            vert(check(:,1),4)=check(:,2);
            %             %

            % Chris' Code %
% % %             % Grab vertices of voxels to be interrogated
% % %             sites2check = unique(octree);
% % %             voxCheck = zeros(size(octree,1),size(octree,2),length(inumsort));
% % %             vertCheck = zeros(length(sites2check),2,length(inumsort));
% % %             rot_vert = vert(:,1:3);
% % %             for n = 1:length(inumsort)
% % %                 % Project vertices of voxels onto camera model
% % %                 [Xi,Yi]=pred([rot_vert(sites2check,1) rot_vert(sites2check,2) rot_vert(sites2check,3)],par,pos,camera);
% % %                 % Check to see if points are inside the heart silhouettte
% % %                 vicinity=ones(length(sites2check),1);
% % %                 vicinity(Xi<lims(inumsort(i),1) | Xi>lims(inumsort(i),2) | Yi<lims(inumsort(i),3) | Yi>lims(inumsort(i),4))=0;
% % %                 % Grab index value in sites2check variable
% % %                 tmp = repmat(reshape(sites2check,[1 1 length(sites2check)]),[size(octree,1) size(octree,2)]);
% % %                 tmp = tmp == repmat(octree,[1 1 size(tmp,3)]);
% % %                 tmp = tmp.*repmat(reshape(1:length(vicinity),[1 1 length(vicinity)]),[size(tmp,1) size(tmp,2) 1]);
% % %                 tmp = sum(tmp,3);
% % %                 % Save to voxel check variable
% % %                 voxCheck(:,:,n) = vicinity(tmp);
% % %                 % Grab inside vertices and mark as inside (1)
% % %                 vici = find(vicinity);
% % %                 vertCheck(vici,1,n) = 1;
% % %                 % Interpolation provides weight to points on the edge
% % %                 Zi=interp2(silh(:,:,inumsort(i)),Xi(vici),Yi(vici));
% % %                 vertCheck(vici,2,n) = Zi;
% % %                 % Rotate cube
% % %                 rot_vert = (Rz*rot_vert')';
% % %             end
% % %             % Grab vertices of voxels to be interrogated
% % %             sites2check = unique(octree);
% % %             % Find vertices that are inside the heart volume
% % %             [inOut] = volumeCheck(silh,vert,sites2check,par,pos,camera,Rz);
% % %             % Update check column (4)
% % %             vert(sites2check,4) = inOut;
% % %             voxCheck = vert(octree,4);
% % %             voxCheck = reshape(voxCheck,[8 length(voxCheck)/8]);
            
            % keep up with the points that have been tested
            % vert(:,4) of ~NaN means that the vertex has been tested
            % vert(:,4) of NaN means that the vertex has not been tested
            
            % Update which level of subdivision algorithm is on
            level=level+1;
        end % End level loop
        close(h)
        
        %% Close hole in top
        vertIn = vert(:,4) == 1;
        vertIn = vertIn.*(1:length(vertIn))';
        vertIn = unique(vertIn);
        vertIn = vertIn(2:end);
        vertOut = vert(:,4) == 0;
        vertOut = vertOut.*(1:length(vertOut))';
        vertOut = unique(vertOut);
        vertOut = vertOut(2:end);
        % % %     figure,scatter3(vert(vertIn),vert(vertIn+size(vert,1)),vert(vertIn+size(vert,1)*2),'g')
        % hold on
        % scatter3(vert(vertOut),vert(vertOut+size(vert,1)),vert(vertOut+size(vert,1)*2),'r')
        
        % capTop = vert(:,3) == 15;          %Find value in plot and set accordingly
        % capTop = capTop.*(1:size(vert,1))';
        % capTop = unique(capTop);
        % capTop = capTop(2:end);
        % vert(capTop,4) = 0;
        
        
        %% Normalization %%
        % No need to normalize but code stays for posterity
        
        % Normalize and save data
        vertn=zeros(size(vert,1),size(vert,2)-1);
        vertn(:,1)=vert(:,1)-min(vert(:,1));
        vertn(:,2)=vert(:,2)-min(vert(:,2));
        vertn(:,3)=vert(:,3)-min(vert(:,3));
        vertn(:,4)=vert(:,4);
        
        scale_x=max(vertn(:,1));
        scale_y=max(vertn(:,2));
        scale_z=max(vertn(:,3));
        
        vertn(:,1)=vertn(:,1)./scale_x;
        vertn(:,2)=vertn(:,2)./scale_y;
        vertn(:,3)=vertn(:,3)./scale_z;
        
        savecom=sprintf('save %s/octvert%s octree vert delt level xdim ydim zdim scale_x scale_y scale_z rr inumsort rsort',savedir,fnametag{1});
        eval(savecom);
        
        scalesfname=sprintf('%s/scales%s.dat',savedir,fnametag{1});
        fid=fopen(scalesfname,'w');
        fprintf(fid,'%4.3f\n',scale_x);
        fprintf(fid,'%4.3f\n',scale_y);
        fprintf(fid,'%4.3f\n',scale_z);
        fclose(fid);
        
        xyzvfname=sprintf('%s/xyzv%s.dat',savedir,fnametag{1});
        fid=fopen(xyzvfname,'w');
        fwrite(fid,vert(:,1:4)','float');   % Edited 3/12/02 from fwrite(fid,vertn','float');
        fclose(fid);
        
        fprintf(parafid,'Saved %d points in %s \n',size(vert,1),xyzvfname);
        fclose(parafid);
        
        return
        
        %% Rendering %%
        disp('Rendering ....')
        
        close all
        clear all
        
        load octvert
        
        xsign=sign(delt(end,1));
        ysign=sign(delt(end,2));
        zsign=sign(delt(end,3));
        
        inside=find(vert(:,4)>=0.80);  % Get the size of volume for meshgrid
        minx=min(xsign.*vert(inside,1)).*xsign;
        maxx=max(xsign.*vert(inside,1)).*xsign;
        miny=min(ysign.*vert(inside,2)).*ysign;
        maxy=max(ysign.*vert(inside,2)).*ysign;
        minz=min(zsign.*vert(inside,3)).*zsign;
        maxz=max(zsign.*vert(inside,3)).*zsign;
        
        [x,y,z]=meshgrid(minx:delt(end,1)/2:maxx,miny:delt(end,2)/2:maxy,minz:delt(end,3)/2:maxz);
        X=[x(:), y(:), z(:)];
        v=griddatan(vert(:,1:3),vert(:,4),X);
        v=reshape(v,size(x));
        vs=smooth3(v);
        
        gg=find(vert(:,4)>=0.90);
        figure
        plot3(vert(gg,1),vert(gg,2),vert(gg,3),'b.')
        hold on
        
        p=patch(isosurface(x,y,z,vs,0.90));
        isonormals(x,y,z,vs,p);
        %set(p,'FaceColor','red','EdgeColor','none');
        set(p,'FaceColor','red');
        
        pcaps=patch(isocaps(x,y,z,vs,0.90));
        set(pcaps,'FaceColor','red');
        
        view(3);
        camlight;
        lighting phong
    else
        status = 0;
    end
else
    status = 0;
end

end