function [int]=electrlen(g);

% The function int=bound2_2(g) calculates the length of the electrode
% from g(1,:) to g(2,:).

% 10.5. 1996 P. Ronkanen and M. Vauhkonen

%w=[1/2,1/2];
%ip=[1/2-1/6*sqrt(3),1/2+1/6*sqrt(3)];
dJt=sqrt((g(2,1)-g(1,1))^2+(g(2,2)-g(1,2))^2); 
%int=0;
% for ii=1:2
%  int=int+w(ii);
% end
%int=int*dJt;
int=dJt;