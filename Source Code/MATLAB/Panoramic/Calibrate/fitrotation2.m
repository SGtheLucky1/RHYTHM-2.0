function mse=fitrotation(theta)

global Rc3p par Rc2 rw r0 Xia Yia n calpts

[result,RR_pos]=rotate3([1 1 1 1],theta,[0 0 0 1],n');     
Rc=(RR_pos*(rw-[r0' 1]')+[r0' 1]');    
Rc3p(:,4)=-Rc3p*Rc;                       
[Xi,Yi]=pred2([calpts(:,4:6).*25.4 ones(size(calpts(:,4:6),1),1)],Rc3p,par);

mse=sqrt(sum((Xi-Xia).^2+(Yi-Yia).^2));
[theta mse]

