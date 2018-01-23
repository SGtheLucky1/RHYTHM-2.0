function [newri,newci,shift]=usershiftindMod(img,ri,ci)
% User-directed pixel shift of indicies in ri,ci for
% optimal alignment with region of interest in image img.
%
% [newri,newci,shift]=usershiftind(img,ri,ci);
%
% Inputs:
% img = mapping mask
% ri = row indices of geometry 
% ci = column indices of geometry
%
% Outputs:
% newri = new row index value
% newci = new column index value
%
% Author: Matthew W. Kay
% Date: 2004
%
% Maintained by: Christopher Gloschat
%
% Revisions Log:
% August 10, 2016 - Added an output (shift) that is populated by the user 
% specified shift for each image.
%
%% Code %%

doneflag=0;
newri=ri;
newci=ci;

while ~doneflag
    h1=figure;
    imagesc(img); hold on; plot(ci,ri,'r.');
    defrshift=0; defcshift=0;
    rshift=0; cshift=0;
    disp('Enter column (x) and row (y) shift in pixels');
    disp('No changes indicate acceptance of bracketed (current) values.');
    while (~isempty(rshift) || ~isempty(cshift))
        rshift=input(sprintf('Enter row (y) shift [%d]: ',defrshift));
        cshift=input(sprintf('Enter column (x) shift [%d]: ',defcshift));
        if ~isempty(rshift)
            defrshift=defrshift+rshift;
            newri=ri+defrshift; 
        end
        if ~isempty(cshift)
            defcshift=defcshift+cshift;
            newci=ci+defcshift;
        end
        figure(h1); hold off;
        imagesc(img);axis square; hold on; plot(newci,newri,'r.');
    end
    goflag=1;
    while goflag
        answer=input(sprintf('Use final values of row shift=%d and column shift=%d? (Y/N): ',defrshift,defcshift),'s');
        if ~isempty(answer)
            if (answer=='Y' || answer=='y'); doneflag=1; goflag=0; end;
            if (answer=='N' || answer=='n'); close(h1); goflag=0; end;
        end
    end
end

% Save the output of the shift
shift = [defrshift defcshift];

disp(' ');
disp('Checking for any new pixel positions that are outside the image ...');
gout=find(newri<1 | newri>size(img,1));
h2=[];
if ~isempty(gout)
    h2=figure; plot(newci,newri,'y.'); set(gca,'YDir','Reverse'); hold on;
    disp('Found indicies outside range of rows (red circles).');
    fprintf('Removing %d sets of indicies\n',length(gout));
    plot(newci(gout),newri(gout),'ro');
    newri(gout)=[];
    newci(gout)=[];
else
    disp('Found no indicies outside range of rows.');
end

gout=find(newci<1 | newci>size(img,2));
if ~isempty(gout)
    if isempty(h2);
        h2=figure; plot(newci,newri,'r.'); set(gca,'YDir','Reverse'); hold on;
    else
        figure(h2);
    end
    disp('Found indicies outside range of columns (black circles).');
    fprintf('Removing %d sets of indicies\n',length(gout));
    plot(newci(gout),newri(gout),'ko')
    newri(gout)=[];
    newci(gout)=[];
else
    disp('Found no indicies outside range of columns.');
end

% set figure axis to match that of the original image
set(gca,'XLim',[1 size(img,1)],'YLim',[1 size(img,2)])

disp(' ');
disp('Press return to close window and exit...');
pause;
close(h1);
if ~isempty(h2)
    close(h2)
end

end

