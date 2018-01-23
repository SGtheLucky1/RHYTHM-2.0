%-----------------------------------
% user defined parameters for pan_textureF.m v11 and above, 11/15/04

% Will usually only need to edit the next 3 lines!

% study='map24';
% runtype='VF';
% runlabel='vf7';
runnumber = '04';
runlabel = 'cam';
% scanlabel='scan4';
scanlabel = 'ABCD';
anum = ['045';'315';'225';'135'];

datadir = '/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/';
studydir = 'Rabbit_2016_0527';
geoName = 'xyzv_try24_povcyl_1mm_corrected_';

%tii=[1:500];       % snapshot number for mapping
%tii=[1:1500];
%tii=[];  % use empty tii to process all frames 

% Everything below here can probably be left as is. 
%---------------------------------------------------

% camera='andor_ixon_2.2mm';
% fps=750;
camera='brainvision_ultimaL';
fps=1000;

minangleDIRECT=95;  % front face angles must be greater than this, use 180 to force max angle assignment
minedgeDIRECT=0.9;
minedgeEDGE=0.55;     % 0.25 is too small!! use at least 0.5!! but 0.25 is same as 0.5; 0.25 is too small to make a difference
edgefiltval=16;      % was 20. reduced for comparison to 100 for spinny images
nuketopdata=0;       % if 1 then no texture mapping on the top of the model. Should be 0 for decimated models

% INPUT FILES
	    
% % datafnames={sprintf('%sChimChim%s/%sChimChim.sif',datadir,runtype,runlabel);  
% %             sprintf('%sSpeed%s/%sSpeed.sif',datadir,runtype,runlabel);
% %             sprintf('%sSpritle%s/%sSpritle.sif',datadir,runtype,runlabel);
% %             sprintf('%sTrixie%s/%sTrixie.sif',datadir,runtype,runlabel);};
datafnames = cell(4,1);
for n = 1:4
    datafnames{n} = sprintf('%s%s/Analysis/%s_%s%s.mat',datadir,studydir,runnumber,runlabel,scanlabel(n));
end

% % mapcalfnames={'../cal/calchimchim.bmp.pospar';   
% %               '../cal/calspeed.bmp.pospar';     
% %               '../cal/calspritle.bmp.pospar';     
% %               '../cal/caltrixie.bmp.pospar';};   
mapcalfnames = cell(4,1);
for n = 1:4
   mapcalfnames{n} = sprintf('%s%s/cal%s_%s.pospar',datadir,studydir,scanlabel(n),anum(n,:));
end

% % maskfnames={sprintf('../%s/chimchim_%smask.bmp',scanlabel,scanlabel);
% %             sprintf('../%s/speed_%smask.bmp',scanlabel,scanlabel);
% % 	    sprintf('../%s/spritle_%smask.bmp',scanlabel,scanlabel);
% % 	    sprintf('../%s/trixie_%smask.bmp',scanlabel,scanlabel);};

% % geomfnames={sprintf('../%s/%s_d_centroids.dat',scanlabel,scanlabel);
% %     sprintf('../%s/%s_d_pts.dat',scanlabel,scanlabel);
% %     sprintf('../%s/%s_d_cells.dat',scanlabel,scanlabel);
% %     sprintf('../%s/%s_d_smoothnorms.dat',scanlabel,scanlabel);
% %     sprintf('../%s/%s_d_neighs.dat',scanlabel,scanlabel);};
geomfnames = {sprintf('%s%s/Geometry/%scentroids.dat',datadir,studydir,geoName);
    sprintf('%s%s/Geometry/%spts.dat',datadir,studydir,geoName);
    sprintf('%s%s/Geometry/%scells.dat',datadir,studydir,geoName);
    sprintf('%s%s/Geometry/%snormals.dat',datadir,studydir,geoName);
    sprintf('%s%s/Geometry/%sneighs.dat',datadir,studydir,geoName);};

	    	        
% OUTPUT FILES
	    
txtyourfname={sprintf('%s_d',runnumber),'cdf'};
edgesfname={sprintf('%sedges_d',runnumber),'dat'};
viewsfname={sprintf('%sviews_d',runnumber),'dat'};
infofname=sprintf('%s_textureinfo.txt',runnumber);

% Create a setup_pan_sigqual.m on the chance that pan_sigqual.m is to be run
fid=fopen('setup_pan_sigqual.m','w');
fprintf(fid,'study=''%s'';\n',studydir);
fprintf(fid,'runlabel=''%s'';\n',runnumber);
fprintf(fid,'fname1=''%s_d.cdf'';\n',runnumber);
fprintf(fid,'fname2=''%s_sigqual.cdf'';\n',runnumber);
fprintf(fid,'fname3=''%s_d.qscale'';\n',runnumber);
fprintf(fid,'winsize=512;  %% 512\n');
fprintf(fid,'nfft=1024;     %% 1024\n');
fprintf(fid,'minfreq=5;   %% Hz, should be higher than the pacing rate\n');
fprintf(fid,'maxfreq=60;  %% Hz\n');
fclose(fid);

% Create a setup_pan_postprocess.m 
fid=fopen('setup_pan_postprocess.m','w');
fprintf(fid,'study=''%s'';\n',studydir);
fprintf(fid,'runlabel=''%s'';\n',runnumber);
fprintf(fid,'fname1=''%s_d.cdf'';\n',runnumber);
fprintf(fid,'fname2=''%s_d+f.cdf'';\n',runnumber);
fprintf(fid,'W1=55;  %% Hz\n');
fprintf(fid,'W2=65;  %% Hz\n');
fclose(fid);

