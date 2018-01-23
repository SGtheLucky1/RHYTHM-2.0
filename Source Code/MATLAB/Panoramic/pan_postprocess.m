clear all
close all

% read filenames and parameters from setup file.
% setup file usually created by setup_pan_texture.m
% setup_pan_texture.m should have been copied into local dir.

setup_pan_postprocess

%---------------------------------------

disp(sprintf('Reading %s ...',fname1));
[data,numparamsorig,txtparams]=readcdf(fname1);
data=data./1000;  % divide by 1000 b/c of cdf file "convention".
fs=1/numparamsorig(3);

disp('Filtering data ...');
[pdata]=lpfilter(data,fs,W1,W2);
pdata(:,1:5)=[];
pdata(:,end-5:end)=[];

disp('Normalizing ...');
pdata=op_norm(pdata,1,100);

disp(sprintf('Saving %s ...',fname2));
numparams3=numparamsorig;
numparams3(3)=1/fs;
numparams3(4)=size(pdata,2);
txtparams(1)={sprintf('%s,%s',txtparams{1},date)};
txtparams(2)={sprintf('%s,lpbutterfilt(%d,%d),normalize',txtparams{2},W1,W2)};
writecdf(fname2,pdata.*1000,numparams3,txtparams);
