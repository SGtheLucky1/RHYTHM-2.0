% extract calibration images from calibration movie for Andor cameras
% 10/8/04
%
clear all
close all

for i=1:4
  if i==1
    ddir=('/cdrom2/SpeedPace');
    fname=('calSpeed.sif');
    BMPfname=('calspeed.bmp');
    PIXfname=('calspeed.pix');
    rr=2;
  elseif i==2
    ddir=('/cdrom2/SpritlePace');
    fname=('calSpritle.sif');
    BMPfname=('calspritle.bmp');
    PIXfname=('calspritle.pix');
    rr=3; 
  elseif i==3
    ddir=('/cdrom2/ChimChimPace');
    fname=('calChimChim.sif');
    BMPfname=('calchimchim.bmp');
    PIXfname=('calchimchim.pix');
    rr=2;
  elseif i==4
    ddir=('/cdrom2/TrixiePace');
    fname=('calTrixie.sif');
    BMPfname=('caltrixie.bmp');
    PIXfname=('caltrixie.pix');
    rr=3;
  end
  
  %--------------------------
  filename=sprintf('%s/%s',ddir,fname);
  fid=fopen(filename);
  if fid==-1
    disp(sprintf('Unable to open %s, file does not exist!',filename));
    return
  end
  fclose(fid);

  [datastruct]=sifopen2(filename);
  data=datastruct.Signal(:,:,1);
  p=datastruct.Signal_lims;
  a=rot90(data(:,:,1));

  vals=a(find(a>min(min(a)) & a<max(max(a))));

  llim=median(vals)-rr*std(vals);
  ulim=median(vals)+rr*std(vals);
  a(find(a<llim))=llim;
  a(find(a>ulim))=ulim;
  a=im_norm(a,0,1);

  imwrite(a,BMPfname,'bmp');
  disp(sprintf('wrote %s',BMPfname));

  fid=fopen(PIXfname,'w');
  fprintf(fid,'%d  %d  %d  %d\n',p);
  fclose(fid);
  disp(sprintf('wrote %s',PIXfname));
  %--------------------------
  
end



