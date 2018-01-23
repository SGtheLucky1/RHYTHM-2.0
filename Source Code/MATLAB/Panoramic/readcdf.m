function [data,numparams,txtparams]=readcdf(fname)
%
% [data,numparams,txtparams]=readcdf(fname);
%
% MWKay, 2003

fid=fopen(fname,'r','b');
l=fgets(fid);
hlinesi=findstr('header_lines=',l)+length('header_lines=');
hlines=str2num(l(hlinesi:end));

for i=2:hlines
  l=lower(fgets(fid));
  if findstr('nchan',l) 
    ii=findstr('nchannels=',l)+length('nchannels=');
    nchannels=str2num(l(ii:end));
  elseif findstr('delta',l) 
    ii=findstr('delta-t=',l)+length('delta-t=');
    deltat=str2num(l(ii:end));
    Fs=1.0/deltat;
  elseif findstr('n_samp',l) 
    ii=findstr('n_samples=',l)+length('n_samples=');
    nsamples=str2num(l(ii:end));
  elseif findstr('word',l) 
    ii=findstr('word_size=',l)+length('word_size=');
    wsize=str2num(l(ii:end));
  elseif findstr('event',l) 
    ii=findstr('event=',l)+length('event=');
    event=str2num(l(ii:end));
  elseif findstr('date',l) 
    ii=findstr('date=',l)+length('date=');
    ddate=l(ii:end-1);          % end-1: don't want newline char
  elseif findstr('data',l) 
    ii=findstr('data_source=',l)+length('data_source=');
    data_source=l(ii:end-1);    % end-1: don't want newline char
  end
end
dat_precision=sprintf('int%d',8*wsize);
dat=fread(fid,(nchannels+1)*nsamples,dat_precision);
fclose(fid);
data=reshape(dat',(nchannels+1),nsamples);
numparams=[hlines,nchannels,deltat,nsamples,wsize,event];
txtparams={ddate data_source};

