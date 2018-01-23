%-----------------------------------
% user defined parameters for pan_coloranatomy.m, 02/21/05
% MWKay

% Will usually only need to edit theses lines...

study='map18';
scanlabel='scan4';

datadir='/disk2/raw_data/map18/Geometry/';

% These filenames have colored (blue) regions of the LV
imnames={'scan4anatomy_000.bmp'
'scan4anatomy_007.bmp'
'scan4anatomy_023.bmp'
'scan4anatomy_033.bmp'
'scan4anatomy_039.bmp'};

imnums=[0 7 23 33 39];

anteriorimnum=33;
posteriorimnum=68;

% Reassign cells to the proper anterior/posterior view:
% use the next 2 lines to read a textfile of cell numbers for reassignment.
% those cells will be reassigned to the anterior (1) or posterior (0) view 
%   indicated by the integer in acellscode.
acellsflag=1;   % 0 if no cells to reassign
acellsfname='acells.dat';
acellscode=0;

% Everything below here can probably be left as is. 
%---------------------------------------------------

camera='watec_with_f8.5';

minangle=95;   % front face angles for a view must be greater than this


% INPUT FILES

spincalfnames={'../cal/Rc.dat';       % Rc first
               '../cal/Par.dat';};    % Par second
	    
geomfnames={sprintf('../%s/%s_d_centroids.dat',scanlabel,scanlabel);
            sprintf('../%s/%s_d_pts.dat',scanlabel,scanlabel);
            sprintf('../%s/%s_d_cells.dat',scanlabel,scanlabel);
	    sprintf('../%s/%s_d_smoothnorms.dat',scanlabel,scanlabel);
	    sprintf('../%s/%s_d_neighs.dat',scanlabel,scanlabel);};
	    	        
% OUTPUT FILES
	    
txtyourfname={sprintf('%s_coloranatomy',scanlabel),'cdf'};
txtyourfname2={sprintf('%s_coloranatomyX',scanlabel),'cdf'};
edgesfname={sprintf('%sedges_coloranatomy',scanlabel),'dat'};
viewsfname={sprintf('%sviews_coloranatomy',scanlabel),'dat'};

