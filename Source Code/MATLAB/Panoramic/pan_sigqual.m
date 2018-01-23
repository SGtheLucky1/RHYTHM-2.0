clear all
close all

% read filenames and parameters from setup file.
% setup file usually created by setup_pan_texture.m
% setup_pan_texture.m should be copied into local dir.

setup_pan_sigqual

%---------------------------------------

disp(sprintf('Reading %s ...',fname1));
[data,numparamsorig,txtparams]=readcdf(fname1);
data=data./1000;  % divide by 1000 b/c of cdf file "convention".

disp(sprintf('Performing signal quality analysis ...',fname1));
fs=1/numparamsorig(3);
pfreq1=zeros(size(data,1),1);
pfreq2=zeros(size(data,1),1);
for i=1:size(data,1)
  [pf1,pf2,freq,p]=pwelch_op(data(i,:),fs,winsize,nfft,minfreq,maxfreq);
  pfreq1(i)=pf1(1);
  pfreq2(i)=pf2;
end

sigstd=zeros(size(data,1),1);
sigmean=zeros(size(data,1),1);
sigp2p=zeros(size(data,1),1);
for i=1:size(data,1)
   sigstd(i)=std(data(i,:));   
   sigmean(i)=mean(data(i,:)); 
   sigp2p(i)=2*sqrt(2)*sigstd(i);
end

gg=find(~pfreq2);
npfreq2=-abs(median(pfreq2)-pfreq2);
npfreq2=op_norm(npfreq2,0,1);
npfreq2(gg)=0;
nsigstd=op_norm(sigstd,0,1);
%qscale=npfreq2.^2 + nsigstd.^2;
qscale=sigp2p;

sigqual=[pfreq1,pfreq2,npfreq2,sigmean,sigstd,sigp2p,qscale];
numparams2=numparamsorig;
numparams2(3)=1;
numparams2(4)=size(sigqual,2);
txtparams(2)={sprintf('%s, %s',study,runlabel)};
writecdf(fname2,sigqual.*1000,numparams2,txtparams)

numparams3=numparamsorig;
numparams3(3)=1;
numparams3(4)=1;
txtparams(2)={sprintf('%s, %s',study,runlabel)};
writecdf(fname3,op_norm(qscale,0,255).*1000,numparams3,txtparams)

