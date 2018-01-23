function scalecdf(filename,scale)

% function a=scalecdf(filename,scale_factor)
% multiply all data in input cdf file by scale_factor.
% If input filename is 'blat.cdf', output is written to
% 'blat_scale.cdf'


fid = fopen(filename,'r','b');
line = fgetl(fid)
nh=sscanf(line,'header_lines=%d');
nh=nh-1;
for n=1:nh,
  line = fgetl(fid);
  if (findstr(line,'nchannels'))
    nchannels = sscanf(line,'nchannels=%d')
  elseif (findstr(line,'n_samples'))
    nsamples = sscanf(line,'n_samples=%d')
  elseif (findstr(line,'word_size'))
    word_size = sscanf(line,'word_size=%d')
	if word_size==2,
      precision='int16';
	else
      precision='int32';
    end
  end
end
a=fread(fid,[nchannels+1,nsamples],precision);
fclose(fid);

a=a.*scale;

[path,file,ext,ver] = fileparts(filename);
newfile = fullfile(path,[file '_scale' '.cdf']);
fid=fopen(newfile,'w','b');
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',nchannels);
fprintf(fid,'delta-t=0.001\n');
fprintf(fid,'n_samples=%d\n',nsamples);
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=whocares\n');
fprintf(fid,'data_source=whocares\n');
fwrite(fid,a,'int32');
fclose(fid);
