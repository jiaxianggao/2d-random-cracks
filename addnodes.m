function [g,HH,Indb1]=addnodes(g,H,bg);

% Function [g,HH]=addnodes(g,H);
% adds the nodes in the middle of each edge of the 2D triangular
% element. g includes the coordiates and H the connectivity.
% HH is the new connectivity matrix. New nodes are added in the end
% of the matrix g.

% M. Vauhkonen 21.10.1999, University of Kuopio, Dept. of Applied
% Physics, Finland

%r=[1 0.7 0.4 0.2 0];
%N=[32 16 8 4 1];
%ntas=2;
%el=[1,1];
%[g,gp,H,E]=cirgrid_eit(r,N,el);

%keyboard

g2=[];
for ii=1:length(H)
  g1=g(H(ii,:),:);
  g21=(g1(1,:)+g1(2,:))/2;
  g22=(g1(2,:)+g1(3,:))/2;
  g23=(g1(1,:)+g1(3,:))/2;
  g2=[g2;g21;g22;g23];
end

g2=intersect(g2,g2,'rows');
g=[g;g2];
inds=zeros(size(g));

HH=zeros(size(H,1),6);
 for ii=1:size(H,1);
   g1=g(H(ii,:),:);
   g21=(g1(1,:)+g1(2,:))/2;
   g22=(g1(2,:)+g1(3,:))/2;
   g23=(g1(1,:)+g1(3,:))/2;
   H1=find(g(:,1)==g21(1) & g(:,2)==g21(2));
    if find(bg==H(ii,1)) & find(bg==H(ii,2))
      inds(H1)=1;
    end
   H2=find(g(:,1)==g22(1) & g(:,2)==g22(2));
    if find(bg==H(ii,2)) & find(bg==H(ii,3))
      inds(H2)=1;
    end
   H3=find(g(:,1)==g23(1) & g(:,2)==g23(2));
    if find(bg==H(ii,3)) & find(bg==H(ii,1))
      inds(H3)=1;
    end 
   HH(ii,:)=[H(ii,1),H1,H(ii,2),H2,H(ii,3),H3];
 end
 
Indb1=[bg;find(inds)];






























