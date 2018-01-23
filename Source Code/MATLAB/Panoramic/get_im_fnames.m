function [fnames,nogo]=get_im_fnames(expname,scanname);
%
% Construct image filenames from user specifications; with smart defaults.
%
% nogo=0 then all files found
% nogo=1 then some or all files NOT found
%
% 2004, MWKay

nogo=0;
defaultbdir=sprintf('/disk2/raw_data/%s/Geometry/%s/',expname,scanname);
bdir=input(sprintf('images are stored in this directory [%s]: ',defaultbdir),'s');
if isempty(bdir); bdir=defaultbdir; end
defaultbfilename=sprintf('%s_',scanname);
bfilename=input(sprintf('basis of image filenames(ie: scan_ is basis for scan_000.pgm )[%s]: ',defaultbfilename),'s');
if isempty(bfilename); bfilename=defaultbfilename; end
defaultsfilename='pgm';
sfilename=input(sprintf('suffix of image filenames(ie: pgm is suffix for scan_000.pgm )[%s]: ',defaultsfilename),'s');
if isempty(sfilename); sfilename=defaultsfilename; end
defaultndigits=3;
ndigits=input(sprintf('number of digits in image filenames(ie: 3 digits in scan_001.pgm)[%d]: ',defaultndigits));
if isempty(ndigits); ndigits=defaultndigits; end
defaultsdigit=0;
sdigit=input(sprintf('starting digit (ie: 0 or 1)[%d]: ',defaultsdigit));
if isempty(sdigit); sdigit=defaultsdigit; end
defaultnimgs=72;
nimgs=input(sprintf('number of images [%d]: ',defaultnimgs));
if isempty(nimgs); nimgs=defaultnimgs; end

for j=1:nimgs
  fnamecom=sprintf('fname=sprintf(''%%s%%s%%0%dd.%%s'',bdir,bfilename,j+(sdigit-1),sfilename);',ndigits);
  eval(fnamecom);
  fnames{j}=fname;
  fid=fopen(fname,'r');
  if fid==-1
    disp(sprintf('unable to open %s',fname));
    nogo=1;
  else
    fclose(fid);
  end
end
