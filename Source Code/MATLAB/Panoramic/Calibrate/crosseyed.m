function data=crosseyed(a,calpts,ifilename,method,zrange,zthresh,visinterp)

% a: calibration image
% calpts: unit normals (cols 1-3) and 3d positions (cols 4-6) of control points
% ifilename: calibration image filename, only for display purposes
% method: 1 for click center, 2 for line intersect
% zrange: zoom this many pixels around initial control point selection for final control point selection
% zthresh: threshold local control point areas by this value in an attempt to hightlight line intersections

if nargin==6
  visinterp=1;
end

[any,anx]=size(a);   % yes, rows=any and columns=anx!
fhpts=figure;
imshow(a);
set(fhpts,'Position',[600 400 600 500]);
hold on;
title('Calibration points','FontSize',14);
fhqc=figure;
imshow(a);
imagesc(a); colormap('gray');
set(fhqc,'Position',[600 400 600 500]);
title('Quality control','FontSize',14);
hold on;

[calM,calN]=size(calpts);
qquit=0;
i=0;
jj=0;
while jj<calM
  jj=jj+1;
  if ~qquit
    skip=0;
    done=0;
    mskip=0;
    while ~done
      disp(sprintf('Image %s\nClick on calibration point #%d:\nUnit normal [%d,%d,%d]\nPosition[%2.2f,%2.2f,%2.2f]',ifilename,jj,calpts(jj,:)));
      fhimpixel=figure;
      prexi1=[]; preyi1=[];
      while isempty(prexi1)
        [prexi1,preyi1,P]=impixel(a);
      end
      xi1(jj)=prexi1; yi1(jj)=preyi1;
      set(fhimpixel,'Position',[600 400 600 500]);
      close(fhimpixel);
      figure(fhpts);
      text(xi1(jj),yi1(jj),sprintf('%d',jj),'Color','b');
      ss=input('OK? [y] or Single Skip? [s] or Multiple Skip? [ms] or Quit? [q]: ','s');
      if (isempty(ss) || ss=='Y' || ss=='y') 
	done=1;
      elseif (ss=='s' || ss=='S') 
	done=1; 
	skip=1; 
	text(xi1(jj),yi1(jj),sprintf('%d',jj),'Color','r'); 
      elseif (strcmp(ss,'ms') || strcmp(ss,'MS'))
        done=1;
	mskip=1;
	skip=1;	
      elseif (strcmp(ss,'q') || strcmp(ss,'Q'))
	done=1; 
	skip=1;
	qquit=1;
      end
    end
    
    if mskip
      skipnum=0;
      mskipend=0;
      while ~mskipend && ~qquit
        defaultprompt='n';
	prompt='g';
	while prompt~='n' && prompt~='N' && prompt~='Y' && prompt~='y'
          prompt=input(sprintf('Skip point #%d: Unit normal [%d,%d,%d] Position[%2.2f,%2.2f,%2.2f]? [n]: ',jj,calpts(jj,:)),'s');
	  if isempty(prompt)
	    prompt=defaultprompt;
	  end;
	end
	if prompt=='Y' || prompt=='y'
	  jj=jj+1;
	  skipnum=skipnum+1;
	else
	  mskipend=1;
	  jj=jj-1;
	end
	if skipnum==1
	  text(xi1(jj-1),yi1(jj-1),sprintf('%d',jj-1),'Color','r');
	end
	if jj>calM
	  qquit=1;
	end
      end
    end

    if ~skip
      i=i+1;
      halfrange=ceil(zrange/2);
      if xi1(jj)-halfrange<1
          xzmin=1; 
      else 
          xzmin=xi1(jj)-halfrange; 
      end
      if xi1(jj)+halfrange>anx
          xzmax=anx; 
      else
          xzmax=xi1(jj)+halfrange; 
      end
      if yi1(jj)-halfrange<1
          yzmin=1; 
      else
          yzmin=yi1(jj)-halfrange; 
      end
      if yi1(jj)+halfrange>any
          yzmax=any; 
      else
          yzmax=yi1(jj)+halfrange; 
      end
      azoom=a(yzmin:yzmax,xzmin:xzmax);
      if visinterp
        interpfactor=4;
        [azoomNx,azoomNy]=size(azoom);
        azoom2=imresize(azoom,[azoomNx,azoomNy]*interpfactor,'bilinear');
      else
        interpfactor=1;
        azoom2=azoom;
      end
      azoom3=zeros(size(azoom2));
      azoom3(find(histeq(azoom2)>=zthresh))=1;

      done=0;
      fhzoom=figure;
      set(fhzoom,'Position',[30 250 1100 700]);
      while ~done
          subplot(1,2,1);
          hold off; imshow(azoom2); imagesc(azoom2); colormap('gray'); axis('square'); title('Click here','FontSize',14); hold on;
          subplot(1,2,2);
          hold off; imshow(azoom2); imagesc(azoom3); colormap('gray'); axis('square'); title('Or Click here','FontSize',14); hold on;
          
          if method==1
              fprintf('Click again on calibration point #%d:',jj);
              [xi(i),yi(i)]=ginput(1);
              subplot(1,2,1); plot(xi(i),yi(i),'ro'); subplot(1,2,2); plot(xi(i),yi(i),'ro');
          elseif method==2
              sprintf('Define calibration point #%d\nby clicking 4 points that define crosshairs.\nConsecutive points will be taken as colinear.',jj)
              precx=[]; while isempty(precx); [precx,precy]=ginput(1); end; cx(1)=precx; cy(1)=precy;
              subplot(1,2,1); plot(cx(1),cy(1),'ro'); subplot(1,2,2); plot(cx(1),cy(1),'ro');
              precx=[]; while isempty(precx); [precx,precy]=ginput(1); end; cx(2)=precx; cy(2)=precy;
              subplot(1,2,1); plot(cx(2),cy(2),'ro'); subplot(1,2,2); plot(cx(2),cy(2),'ro');
              precx=[]; while isempty(precx); [precx,precy]=ginput(1); end; cx(3)=precx; cy(3)=precy;
              subplot(1,2,1); plot(cx(3),cy(3),'ro'); subplot(1,2,2); plot(cx(3),cy(3),'ro');
              precx=[]; while isempty(precx); [precx,precy]=ginput(1); end; cx(4)=precx; cy(4)=precy;
              subplot(1,2,1); plot(cx(4),cy(4),'ro'); subplot(1,2,2); plot(cx(4),cy(4),'ro');
              absdiff=abs([cx(2)-cx(1),cy(2)-cy(1)]);   % Perform schnagagins to avoid infinite fits
              maxdiffi1=find(absdiff==max(absdiff));
              if maxdiffi1==1                           % then fit y=my1(x)+by1 for line1
                  [Py1,Sy1]=polyfit(cx(1:2),cy(1:2),1);
                  my1=Py1(1);
                  by1=Py1(2);
                  ppx=cx(1)-(cx(2)-cx(1))/10:(cx(2)-cx(1))/10:cx(2)+(cx(2)-cx(1))/10;
                  ppy=polyval(Py1,ppx);
                  subplot(1,2,1); plot(ppx,ppy,'r-'); subplot(1,2,2); plot(ppx,ppy,'r-');
              else %then maxdiffi1=2                    % then fit x=mx1(y)+bx1 for line1
                  [Px1,Sx1]=polyfit(cy(1:2),cx(1:2),1);
                  mx1=Px1(1);
                  bx1=Px1(2);
                  ppy=cy(1)-(cy(2)-cy(1))/10:(cy(2)-cy(1))/10:cy(2)+(cy(2)-cy(1))/10;
                  ppx=polyval(Px1,ppy);
                  subplot(1,2,1); plot(ppx,ppy,'r-'); subplot(1,2,2); plot(ppx,ppy,'r-');
              end
              
              absdiff=abs([cx(4)-cx(3),cy(4)-cy(3)]);   % Perform schnagagins to avoid infinite fits
              maxdiffi2=find(absdiff==max(absdiff));
              if maxdiffi2==1                           % then fit y=my2(x)+by2 for line2
                  [Py2,Sy2]=polyfit(cx(3:4),cy(3:4),1);
                  my2=Py2(1);
                  by2=Py2(2);
                  ppx=cx(3)-(cx(4)-cx(3))/10:(cx(4)-cx(3))/10:cx(4)+(cx(4)-cx(3))/10;
                  ppy=polyval(Py2,ppx);
                  subplot(1,2,1); plot(ppx,ppy,'r-'); subplot(1,2,2); plot(ppx,ppy,'r-');
              else %then maxdiffi2=2                    % then fit x=mx2(y)+bx2 for line2
                  [Px2,Sx2]=polyfit(cy(3:4),cx(3:4),1);
                  mx2=Px2(1);
                  bx2=Px2(2);
                  ppy=cy(3)-(cy(4)-cy(3))/10:(cy(4)-cy(3))/10:cy(4)+(cy(4)-cy(3))/10;
                  ppx=polyval(Px2,ppy);
                  subplot(1,2,1); plot(ppx,ppy,'r-'); subplot(1,2,2); plot(ppx,ppy,'r-');
              end
              
              if maxdiffi1==1 && maxdiffi2==1       % then y=my1(x)+by1 and y=my2(x)+by2
                  xi(i)=(by2-by1)/(my1-my2);
                  yi(i)=my2*xi(i)+by2;
              elseif maxdiffi1==1 && maxdiffi2==2   % then y=my1(x)+by1 and x=mx2(y)+bx2
                  yi(i)=(my1*bx2+by1)/(1-my1*mx2);
                  xi(i)=mx2*yi(i)+bx2;
              elseif maxdiffi1==2 && maxdiffi2==1   % then x=mx1(y)+bx1 and y=my2(x)+by2
                  yi(i)=(my2*bx1+by2)/(1-my2*mx1);
                  xi(i)=mx1*yi(i)+bx1;
              elseif maxdiffi1==2 && maxdiffi2==2   % then x=mx1(y)+bx1 and x=mx2(y)+bx2
                  yi(i)=(bx2-bx1)/(mx1-mx2);
                  xi(i)=mx2*yi(i)+bx2;
              end
              
              subplot(1,2,1); plot(xi(i),yi(i),'go'); subplot(1,2,2); plot(xi(i),yi(i),'go');
          end  % end method 2
          
          xi(i)=(xi(i)-0.5)/interpfactor;      % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
          yi(i)=(yi(i)-0.5)/interpfactor;      % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
          ss=input('OK? [y] ','s');
          if (isempty(ss) || ss=='Y' || ss=='y')
              done=1;
          end
      end
      close(fhzoom);
      xi(i)=xi(i)+xzmin-1;           % xi range is 0:anx-1
      yi(i)=yi(i)+yzmin-1;           % yi range is 0:any-1
      
      figure(fhqc)
      plot(xi(i)+0.5,yi(i)+0.5,'ro');  % image command plots images in range [0.5 N+0.5 0.5 M+0.5]
      
      if (i==1)
          data=[calpts(jj,4:6), xi(i), yi(i), calpts(jj,1:3)];
      else
          data=[data;calpts(jj,4:6), xi(i), yi(i), calpts(jj,1:3)];
      end;
      
    end % end skip condition
  end   % end quit condition
end     % end control point for loop

sprintf('Press return to exit crosseyed...')
pause
close(fhqc);
close(fhpts);

