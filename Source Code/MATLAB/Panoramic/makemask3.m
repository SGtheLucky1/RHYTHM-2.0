%
% Creates mapping masks using geometry.
%
% revised 9/22/03
% revised for Andor cameras 10/15/04
%
% MWKay
%

clear all
close all

camera='andor_ixon_2.2mm';

geomfnames={'scan4_d_centroids.dat'};
mapcalfnames={'../cal/calchimchim.bmp.pospar';   
              '../cal/calspeed.bmp.pospar';     
              '../cal/calspritle.bmp.pospar';     
              '../cal/caltrixie.bmp.pospar';};   
savefnames={'chimchim_scan4mask.bmp';
            'speed_scan4mask.bmp';
	    'spritle_scan4mask.bmp';
	    'trixie_scan4mask.bmp';};

nr=128;
nc=128;


%---------------------------------
% Load geometry data
[centroids,numparams,txtparams]=readdat(char(geomfnames(1)));
	    
%-------------------------------
% Build mapping camera model	      
Rcmap=NaN.*ones(4,4,length(mapcalfnames));
Parmap=NaN.*ones(8,length(mapcalfnames));
xyzmap=NaN.*ones(length(mapcalfnames),4);
pixtot=0;
for i=1:length(mapcalfnames)
  clear Xi;
  clear Yi;
  
  pospar=load(char(mapcalfnames(i)));
  wa=pospar(4)*pi/180;  % omega, Wx rotation (world basis, x axis) 
  pa=pospar(5)*pi/180;  % psi, y rotation  (world basis, y axis)  
  ra=pospar(6)*pi/180;  % kappa, z rotation (world basis, z axis) 
  cw=cos(wa); sw=sin(wa);
  cp=cos(pa); sp=sin(pa);
  cr=cos(ra); sr=sin(ra);
  R=zeros(3,3);

  R(:,1)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw]';             % build the 4x4 R matrix
  R(:,2)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw]';
  R(:,3)=[-sp cp*sw cp*cw]';
  Rcmap(:,:,i)=[[R'; [0 0 0]] [pospar(1) pospar(2) pospar(3) 1]'];  % store each R in a 3d Rc matrix for easy referencing later

  xyzmap(i,:)=(inv(Rcmap(:,:,i))*[0 0 0 1]')';   % compute xyz positions of orbit from calibration matricies
  Parmap(:,i)=pospar(7:14);
  
  mask=zeros(nr,nc);
  [Xi,Yi]=pred2([centroids(:,1) centroids(:,2) centroids(:,3) ones(size(centroids(:,1),1),1)],Rcmap(:,1:4,i),Parmap(1:8,i),camera);
  out=find(Xi>nc|Xi<1|Yi>nr|Yi<1);
  if ~isempty(out)
    Xi(out)=[]; Yi(out)=[];
  end
  mask(sub2ind([nr,nc],round(Yi),round(Xi)))=1;  
  mask=bwmorph(mask,'close',Inf);
  imwrite(double(mask),char(savefnames(i)),'bmp');
  pixtot=pixtot+length(find(mask));
  disp(sprintf('Wrote %s, %d pixels in mask',char(savefnames(i)),length(find(mask))));
end

disp(sprintf('%d total pixels',pixtot));
