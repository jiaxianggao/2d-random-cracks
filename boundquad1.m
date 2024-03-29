function int=boundquad1(g);

%boundquad1 calculates the boundary integral
% of one quadratic basis function over the curve defined by the coordinates in g.
% Function int=boundquad1(g,ind) calculates the boundary integral
% of one quadratic basis function over the curve defined by the coordinates in g.
%
% INPUT
%
% g = integration points
%
% OUTPUT
%
% int = value of the integral

% 10.5. 1996 P. Ronkanen and M. Vauhkonen
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


w=[1/2,1/2];
ip=[1/2-1/6*sqrt(3),1/2+1/6*sqrt(3)];
int=0;
 for ii=1:2
  S=[2*ip(ii)^2-3*ip(ii)+1; ...
     -4*ip(ii)^2+4*ip(ii); ...
     2*ip(ii)^2-ip(ii)];
  dJt=sqrt((g(1,1)*(4*ip(ii)-3)+g(2,1)*(4-8*ip(ii))+g(3,1)*(4*ip(ii)-1))^2+ ...
            (g(1,2)*(4*ip(ii)-3)+g(2,2)*(4-8*ip(ii))+g(3,2)*(4*ip(ii)-1))^2);
  int=int+w(ii)*S*dJt;
 end





