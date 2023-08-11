function [Node]=MakeNode(Element,Nodelist,g);

% Function [Node]=MakeNode(Element,g);
% computes the Node data for MeshData.
% Node is a structure including all the nodal coordinates and
% for each node there is information to which nodes (NodeConnection) 
% and elements (ElementConnection) the node is
% connected.  

% M. Vauhkonen, University of Kuopio, Finland, 11.8.1999 

[rg,cg]=size(g);
msE=max(size(Element));
for ii=1:rg
 ElementConnection=[];
 Node(ii).Coordinate=[g(ii,:)];
  for jj=1:msE
   if find(Element(jj).Topology==ii)
    ElementConnection=[ElementConnection,jj];
   end
  end
  Node(ii).ElementConnection=ElementConnection;
 Nc=[];
 [I,J]=find(Nodelist==ii);
  for kk=1:size(I,1)
    if J(kk)==2, i2=1; else i2=2;end
     if kk==1
       Nc=[Nc,Nodelist(I(kk),i2)];
     else
       if find(Nc==Nodelist(I(kk),i2));
       else
        Nc=[Nc,Nodelist(I(kk),i2)];
       end  
     end  
  end
  Node(ii).NodeConnection=Nc;
end
  