function [ind] = findnearestpoint( g,pts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ind = [];
for ii = 1:size(pts,1)
   
  dist_ii = (g(:,1)-pts(ii,1)).^2 + (g(:,2)-pts(ii,2)).^2;
  [mindist_ii,ind_ii] = min(dist_ii);
  ind_ii = ind_ii(1);
  ind = [ind; ind_ii];
    
end

end

