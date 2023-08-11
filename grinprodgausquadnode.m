function int=grinprodgausquadnode(g,sigma);

% The function int=grinprodgausquadnode(g,sigma);
% calculates the gradient part in the linear 
% FEM in EIT.

% P. Ronkanen and M. Vauhkonen 10.5. 1996


w=[1/6*ones(3,1)];
ip=[1/2 0;1/2 1/2;0 1/2];


int=0;
 for ii=1:3
  S=[1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
  L=[4*(ip(ii,1)+ip(ii,2))-3, -8*ip(ii,1)-4*ip(ii,2)+4, ...
     4*ip(ii,1)-1, 4*ip(ii,2), 0, -4*ip(ii,2); ...
     4*(ip(ii,1)+ip(ii,2))-3, -4*ip(ii,1), ...
     0, 4*ip(ii,1), 4*ip(ii,2)-1, -8*ip(ii,2)-4*ip(ii,1)+4];
  Jt=L*g;
  iJt=inv(Jt);
  dJt=abs(det(Jt));
  G=iJt*L;
  int=int+w(ii)*(S'*sigma)*G'*G*dJt;
 end








