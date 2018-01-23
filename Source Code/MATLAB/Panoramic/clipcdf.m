function clipcdf(cdffile,clipfile)

% function clipcdf(cdffile,clipfile)
% remove subset of channels from a cdf file.
% clipfile (e.g., the output of clipRegion) has
% two columns. The second column lists the channels that
% are to be *retained* from the input file. The first colum
% is the number of each retained channel in the output file.
% If the input file is named 'blat.cdf', output goes to 
% 'blat_clipped.cdf'


fid = fopen(cdffile,'r','b');
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

map = load(clipfile);

nnew = size(map,1);
out(nnew+1,nsamples)=0;
for n=1:nnew,
	nold = map(n,2)+2;
	out(n+1,:)=a(nold,:);
end

[path,file,ext,ver] = fileparts(cdffile);
newfile = fullfile(path,[file '_clipped' '.cdf']);
fid=fopen(newfile,'w','b');
fprintf(fid,'header_lines=8\n');
fprintf(fid,'nchannels=%d\n',nnew);
fprintf(fid,'delta-t=0.001\n');
fprintf(fid,'n_samples=%d\n',nsamples);
fprintf(fid,'word_size=4\n');
fprintf(fid,'event=0\n');
fprintf(fid,'date=whocares\n');
fprintf(fid,'data_source=whocares\n');
fwrite(fid,out,'int32');
fclose(fid);
