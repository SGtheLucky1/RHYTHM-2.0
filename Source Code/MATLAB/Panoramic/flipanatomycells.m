function flipanatomycells(cdffile,cellsfile,anatcode,cdffile2);
% flips anatomy cells according to anatcode
% 
%
% anatcode: 1 for RV/LV
%           2 for anterior/posterior
%           3 for apex/base
% anatcode corresponds to the sample sets in cdffile.
%
% flipanatomycells(cdffile,cellsfile,anatcode,cdffile2);
%
% flipanatomycells('scan3_coloranatomy.cdf','vcells.dat',1,'scan3_coloranatomyX.cdf');

[data,numparams,txtparams]=readcdf(cdffile);
writecdf(sprintf('%s.OLD',cdffile),data,numparams,txtparams);
disp(sprintf('copied %s into %s.OLD',cdffile,cdffile));

eval(sprintf('load %s',cellsfile));
[path,dataname,ext]=fileparts(cellsfile);
eval(sprintf('thecells=%s(:,1)+1;',dataname));  % increase by 1 b/c vtk indicies start with zero
thecells=unique(thecells);
thecells=thecells+1;   % increase by 1 again b/c first row of data from cdf files is zero

data=data./1000;
% Reassign LV/RV cells according to cellsfile and anatcode
data(thecells,anatcode)=~data(thecells,anatcode);

% Reassign anatomy codes
for i=2:size(data,1)
  binnum=[data(i,1) data(i,2) data(i,3)];
  data(i,4)=dot(binnum,[1 2 4]);
end

writecdf(cdffile,data.*1000,numparams,txtparams);
disp(sprintf('wrote %s',cdffile));
numparams(4)=1;  % num of channels
writecdf(cdffile2,data(:,4).*1000,numparams,txtparams);
disp(sprintf('wrote %s',cdffile2));

