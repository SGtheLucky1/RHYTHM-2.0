function mse=fitrotation(rotorigin)

global C0
global n

r0=rotorigin(1:3);
r=n;

p=C0(:,find(C0(1,:)==0));
p(1)=[];
q=C0;
q(:,find(C0(1,:)==0))=[];

for i=1:size(q,2)

  theta=q(1,i);
  ct=cos(theta*pi/180);
  st=sin(theta*pi/180);
  r=r./sqrt(r(1)^2+r(2)^2+r(3)^2);  % normalize r

  R=[ct+(1-ct)*r(1)^2  (1-ct)*r(1)*r(2)+r(3)*st  (1-ct)*r(1)*r(3)-r(2)*st
     (1-ct)*r(1)*r(2)-r(3)*st  ct+(1-ct)*r(2)^2  (1-ct)*r(2)*r(3)+r(1)*st
     (1-ct)*r(1)*r(3)+r(2)*st  (1-ct)*r(2)*r(3)-r(1)*st  ct+(1-ct)*r(3)^2];

  q_pred(:,i)=R*(p-r0)+r0;
  eradius(i)=sqrt( (q(2,i)-q_pred(1,i))^2 + (q(3,i)-q_pred(2,i))^2 + (q(4,i)-q_pred(3,i))^2 );
end

%mse=sum((eradius).^2);
%mse=mean((eradius).^2);
mse=sum((eradius));
  

