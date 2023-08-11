function [g1,H1,Indb1,E1,Element1,Nodelist1,Node1] = ModifyMesh2nd(g1,H1,Indb1,E1)

% The mesh is modified for 2nd order basic functions

[g1,H1,Indb1]=addnodes(g1,H1,Indb1);
%if some of the boundary nodes are not in boundary
  %figure,trimesh(H1,g1(:,1),g1(:,2),zeros(size(g1,1),1))    
  %title('click the marks that are not on the boundary');
  %view(2),axis equal,hold on                  
  %ttt1=g1(Indb1,:);       
  %plot(ttt1(:,1),ttt1(:,2),'r+'),drawnow;            
  %[xi,yi] = ginput;
  %rm=[];
  %for k=1:length(xi)
  %   remov(k) = find(((g1(:,1)-xi(k)).^2+(g1(:,2)-yi(k)).^2)==min((g1(:,1)-xi(k)).^2+(g1(:,2)-yi(k)).^2));
  %   rm=[rm find(Indb1==remov(k))];
  %end;

  %Indb1(rm)=[];

[Element1,Nodelist1]=MakeElement2nd(H1,Indb1,E1);
[Node1]=MakeNode(Element1,Nodelist1,g1);




