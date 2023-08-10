function [Gamma_pr] = GaussianDistanceSmoothnessPrior(g,sigvar,corrlength)

sN = length(g);
c = .001*sigvar;
a = sigvar - c;
b = sqrt(-corrlength^2/(2*log(.01)));

%Gamma_pr = sparse(sN,sN);
Gamma_pr = zeros(sN,sN);

for ii = 1:sN
  for jj = ii:sN

  dist_ij = norm(g(ii,:)-g(jj,:));
  gamma_ij = a*exp(-dist_ij^2/(2*b^2));
  if ii == jj
    gamma_ij = gamma_ij + c;
  end

  Gamma_pr(ii,jj) = gamma_ij; 
  Gamma_pr(jj,ii) = gamma_ij; 

  end
end







