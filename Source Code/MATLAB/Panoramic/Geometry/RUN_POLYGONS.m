disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%   DRIVER PROGAMME FOR CONVEX POLYGON INTERSECTION   %%')
disp('    %%                                                     %%')
disp('    %%   (c) Martin Kaeser, TU Muenchen, 07.02.2002        %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

clear, close all;

disp(sprintf('\n\t DRAW POLYGONS P and Q '));
[p,q,n,m] = drawpoly(1);

disp(sprintf('\n\t CHECKING CONVEXITY ... '));
po = checkpoly(p);
qo = checkpoly(q);
if ( length(po)<1 | length(qo)<1 )
    disp(sprintf('\t *** POLYGONS NOT CONVEX !!! ***\n\n'));
    break
end
disp(sprintf('\t OKAY !'));

disp(sprintf('\n\t COMPUTE INTERSECTION ... '));
[pint,k,A] = convex_intersect(po,qo);
disp(sprintf('\t OKAY !'));

figure
fill(p(:,1),p(:,2),'g'),hold on,axis([-1 1 -1 1])
fill(q(:,1),q(:,2),'b')
if k > 2
    fill(pint(:,1),pint(:,2),'r')
end
axis square, grid on
title('POLYGON INTERSECTION')

disp(sprintf('\n\t AREA OF P:      %g',polyarea(p(:,1),p(:,2)) ));
disp(sprintf(  '\t AREA OF Q:      %g',polyarea(q(:,1),q(:,2)) ));
disp(sprintf(  '\t AREA OF PINT:   %g\n\n',A ));
    
    
        
            
        
        


