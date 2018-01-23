function writepandat(data,sn,save_fname,orig_datafilename)
% Saves processed panoramic data
%   processed data: data
%   raw signal to noise ratios: sn
%
%  writepandat(data,sn,save_fname);
%  
% 2002, MWKay

xdim=size(data,1);
ydim=size(data,2);
tdim=size(data,3);

fid=fopen(save_fname,'w');
fprintf(fid,'header lines: 7\n');
fprintf(fid,'xdim: %d\n',xdim);
fprintf(fid,'ydim: %d\n',ydim);
fprintf(fid,'tdim: %d\n',tdim);
fprintf(fid,'datatype: float\n');
fprintf(fid,'source file: %s\n',orig_datafilename);
fprintf(fid,'created: %s\n',datestr(now));
fwrite(fid,data,'float');
fwrite(fid,sn,'float');
fclose(fid);

