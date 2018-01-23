function [data,sn,numparams,txtparams]=readpandat(fname)
%
% ie: [data,sn,numparams,txtparams]=readpandat('pan02a/bin/map02a010A_2-601.bin');
%
% Read panoramic mapping data files. Headers are expected.
% 
% MWKay, 2002
%

fid=fopen(fname,'r','n');
l=lower(fgets(fid));
hlinesi=findstr('header lines:',l)+length('header lines:');
hlines=str2num(l(hlinesi:end));

for i=2:hlines
  l=lower(fgets(fid));
  if findstr('xdim',l) 
    ii=findstr('xdim:',l)+length('xdim:')+1;
    xdim=str2num(l(ii:end-1));
  elseif findstr('ydim',l) 
    ii=findstr('ydim:',l)+length('ydim:')+1;
    ydim=str2num(l(ii:end-1));
  elseif findstr('tdim',l) 
    ii=findstr('tdim:',l)+length('tdim:')+1;
    tdim=str2num(l(ii:end-1));
  elseif findstr('datatype',l) 
    ii=findstr('datatype:',l)+length('datatype:')+1;
    datatype=l(ii:end-1);
  elseif findstr('source',l) 
    ii=findstr('source file:',l)+length('source file:')+1;
    sourcefilename=l(ii:end-1);
  elseif findstr('created',l) 
    ii=findstr('created:',l)+length('created:')+1;
    created=l(ii:end-1);
  end
end

dat=fread(fid,[xdim,ydim*tdim],datatype);
sn=fread(fid,[xdim,ydim],datatype);
fclose(fid);

data=reshape(dat,xdim,ydim,tdim);
numparams=[xdim ydim tdim];
txtparams={sourcefilename created datatype};


