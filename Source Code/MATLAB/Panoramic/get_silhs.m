function [silh1,lims1,centroids1,orientations1,areas1,thresharr,analyzed]=get_silhs(bdir,bfilename,sfilename,ndigits,sdigit,n_images)
%
% Revised for increased user-directed analysis with consistent
% saving of analyzed silhouettes as well as user-selected image numbers
% for analysis. 12/08/2004, MWKay
%
% analyzed set to 1 if the corresponding image has been analyzed, zero otherwise.
%

fid=fopen('thresharr.dat');
if fid~=-1
    disp('Found thresharr.dat!');
    fclose(fid);
    loadthresh='';
else
    disp('Could not find thresharr.dat!');
    loadthresh='N';
end

go=1;
while go
    ab=input('Is the object ABOVE the threshold (bright object, dark background) \nor BELOW the threshold (dark object, bright background) (A/B) [A]: ','s');
    if (isempty(ab) || ab=='a' || ab=='A')
        go=0;
        aabb=1;
    elseif (ab=='b' || ab=='B')
        go=0;
        aabb=0;
    end
end

% Write the image filenames to a textfile
imnamesfid=fopen(sprintf('%s_fnames.txt',bfilename),'w');
for i=1:n_images
    fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,i+(sdigit-1),sfilename);',ndigits);
    eval(fnamecom);
    fprintf(imnamesfid,'%s\n',fname);
end
fclose(imnamesfid);

% Determine the default threshold OR load thresholds from file thresharr.dat
if isempty(loadthresh)
    goflag=1;
    while (goflag)
        loadthresh=input('Load thresholds from thresharr.dat? (Y/N) [Y]: ','s');
        if isempty(loadthresh); loadthresh='Y'; end;
        if (loadthresh=='Y' || loadthresh=='y' || loadthresh=='N' || loadthresh=='n')
            goflag=0;
        end
    end
end

if (loadthresh=='Y' || loadthresh=='y')
    load thresharr.dat
else % determine default threshold
    figh=figure;
    j=1;
    fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,j+(sdigit-1),sfilename);',ndigits);
    eval(fnamecom);
    %   a=readtiff(fname);
    a = double(imread(fname));
    a = a(:,:,1)/max(a(:));
    subplot(1,2,1); subimage(a); title(sprintf('%s image %d',bfilename,j+(sdigit-1)));
    subplot(1,2,2); imhist(a);
    set(figh,'Position',[23 585 1195 365]);
    
    goflag=1;
    def_thresh=0.35;
    while goflag
        diffflag=1;
        while diffflag
            new_def_thresh=input(sprintf('Enter default threshold [%4.3f]: ',def_thresh));
            if isempty(new_def_thresh)
                diffflag=0;
            elseif (new_def_thresh>=0 && new_def_thresh<=1.0)
                diffflag=0;
            end
        end
        if ~isempty(new_def_thresh)
            def_thresh=new_def_thresh;
        end
        thresharr=zeros(1,n_images)+def_thresh;
        [defsilh,deflims,defcentroid,deforientation,defarea]=getsilh(a,def_thresh,aabb);
        deffig=figure;
        subplot(1,3,1), imagesc(a); title(sprintf('%s image %d',bfilename,j+(sdigit-1)));
        subplot(1,3,2), imagesc(defsilh); title(sprintf('Threshold level: %2.3f',def_thresh));
        subplot(1,3,3), imhist(a);
        set(gcf,'Position',[11 665 1235 267]);
        
        olinefig=figure;
        set(olinefig,'Position',[387 60 873 624]);
        outline=bwperim(defsilh,8);
        [or,oc]=find(outline);
        imagesc(a); colormap('gray'); hold on;
        plot(oc,or,'y.');
        
        answflag=1;
        while answflag
            answer=input(sprintf('Keep default threshold at %4.3f? (Y/N): ',def_thresh),'s');
            if ~isempty(answer)
                if (answer=='Y' || answer=='y' || answer=='N' || answer=='n')
                    answflag=0;
                end
            end
        end
        if (answer=='Y' || answer=='y')
            goflag=0;
        end
        close(deffig);
        close(olinefig); clear olinefig;
    end
    disp(sprintf('Default threshold set to %4.3f!',def_thresh));
    close(figh);
end

fid=fopen('silhs1.mat');
if fid~=-1
    fclose(fid);
    disp(sprintf('Loading previously analyzed data from silhs1.mat ....'));
    load silhs1
    if (loadthresh~='Y' & loadthresh~='y')
        disp(sprintf('Resetting default thresholds for images not analyzed to %4.3f.',def_thresh));
        thresharr(find(~analyzed))=def_thresh;
    end
else
    disp(sprintf('Previously analyzed data not found. Initializing data arrays...'));
    % read the first image file to get the image dimensions.
    % do that so that the silh1 matrix can be initialized.
    j=1;
    fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,j+(sdigit-1),sfilename);',ndigits);
    eval(fnamecom);
%     a=readpgm(fname);
    a = double(imread(fname));
    a = a(:,:,1)/max(a(:));
    
    silh1=zeros(size(a,1),size(a,2),n_images);
    lims1=zeros(n_images,4);
    centroids1=zeros(n_images,2);
    orientations1=zeros(1,n_images);
    areas1=zeros(1,n_images);
    analyzed=zeros(1,n_images);  % set this to 1 if the image has been analyzed, zero otherwise.
end

notdone=1;
j=1;            % j ranges from 1 to n_images
imnum=sdigit-1;   % image numbers range from sdigit to n_images+(sdigit-1)
increment=1;
%---------------------------------------------
%for j=1:n_images
while notdone
    if ~exist('barwin')
        barwin=figure;
    else
        figure(barwin);
    end
    bar([sdigit:n_images+sdigit-1],analyzed); hold on;
    plot([sdigit:n_images+sdigit-1],analyzed,'r.'); axis([sdigit-1 n_images+sdigit -0.1 max(analyzed)+0.1]);
    disp('------------------------------');
    disp(sprintf('%d images left to analyze!',length(find(~analyzed))));
    
    if increment
        imnum=imnum+1;
        if (imnum>n_images+(sdigit-1)); imnum=n_images+(sdigit-1) ; end;
    end
    newimnum=-3;
    while (newimnum<-2 | newimnum>n_images+(sdigit-1))
        disp(sprintf('Enter image number to analyze (%d to %d)',sdigit,n_images+(sdigit-1)));
        disp('   or -1 to list images not analyzed');
        newimnum=input(sprintf('   or -2 to quit [%d]: ',imnum));
        if isempty(newimnum); newimnum=imnum; end
    end
    
    goflag=1;
    if newimnum==-1
        numstring='';
        disp(' ');
        disp(sprintf('Images yet to be analyzed (n=%d):',length(find(~analyzed))));
        for i=1:n_images
            if ~analyzed(i)
                numstring=[numstring,sprintf('%d, ',i+(sdigit-1))];
            end
        end
        disp(sprintf('%s',numstring));
        disp(' ');
        goflag=0;
        increment=0;
    elseif newimnum==-2
        disp('Exiting ....');
        goflag=0;
        notdone=0;  % Exit!
    else
        imnum=newimnum;
        j=imnum+1-sdigit;
        increment=1;
    end
    
    fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,j+(sdigit-1),sfilename);',ndigits);
    eval(fnamecom);
%     a=readpgm(fname);
    a = double(imread(fname));
    a = a(:,:,1)/max(a(:));
    while goflag
        [silh1(:,:,j),lims1(j,:),centroids1(j,:),orientations1(j),areas1(j)]=getsilh(a,thresharr(j),aabb);
        if ~exist('anwin')
            anwin=figure;
        else
            figure(anwin);
        end
        subplot(1,3,1), subimage(a);
        titlecom=sprintf('figtitle=sprintf(''%%s image %%0%dd'',bfilename,imnum);',ndigits);
        eval(titlecom); title(figtitle);
        subplot(1,3,2), subimage(silh1(:,:,j)); title(sprintf('Threshold level: %2.3f',thresharr(j)));
        subplot(1,3,3), imhist(a);
        %set(gcf,'Position',[11 665 1235 267]);
        
        if ~exist('olinefig')
            olinefig=figure;
        else
            figure(olinefig);
        end
        % set(olinefig,'Position',[387 60 873 624]);
        outline=bwperim(silh1(:,:,j),8);
        [or,oc]=find(outline);
        imagesc(a); colormap('gray'); hold on;
        plot(oc,or,'y.');
        titlecom=sprintf('figtitle=sprintf(''image %%0%dd'',imnum);',ndigits);
        eval(titlecom);
        title(figtitle);
        
        disp(sprintf('Current image: %s',fname));
        if j==1
            diffflag=1;
            while diffflag
                difflevel=input(sprintf('Enter a different threshold level or press enter to continue [%4.3f]: ',thresharr(j)));
                if isempty(difflevel)
                    diffflag=0;
                elseif (difflevel>=0 & difflevel<=1.0)
                    diffflag=0;
                end
            end
        else
            disp(sprintf('Previous image threshold (#%d) was %4.3f.',imnum-1,thresharr(j-1)));
            diffflag=1;
            while diffflag
                difflevel=input(sprintf('Enter a different threshold level or press enter to continue [%4.3f]: ',thresharr(j)));
                if isempty(difflevel)
                    diffflag=0;
                elseif (difflevel>=0 & difflevel<=1.0)
                    diffflag=0;
                end
            end
        end
        if isempty(difflevel)
            goflag=0;
            analyzed(j)=analyzed(j)+1;  % marked as analyzed b/c getsilh was called and the user liked the threshold
        else
            thresharr(j)=difflevel;
        end
    end
end
%---------------------------------------------

disp('Saving silhs1.mat....');
save silhs1 silh1 lims1 centroids1 orientations1 areas1 thresharr analyzed

disp('Saving thresharr.dat....');
fid=fopen('thresharr.dat','w');
for i=1:length(thresharr)
    fprintf(fid,'%d\n',thresharr(i));
end
fclose(fid);
if exist('anwin'); close(anwin); end
if exist('barwin'); close(barwin); end
if exist('olinefig'); close(olinefig); end



