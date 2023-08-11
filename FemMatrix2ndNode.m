function [A]=FemMatrix2nd(Node,Element,sig,z,C);

% function [AA]=FemMatrix2nd(Node,Element,sig,z,C);
% supplies matrices
% needed for the FEM solution in the EIT forward problem.
%
% Node=node data
% Element=element data
% sig=complex conductivity
% z=complex contact impedance
% C=voltage reference matrix

% M. Vauhkonen 11.5.1994, modified from the version of J. Kaipio
% 25.4.1994. Modified 5.9.1994 by M. Vauhkonen for EIT.
% Modified for EIDORS by M. Vauhkonen 20.10.1999
% University of Kuopio, Dept. of Applied Physics, Finland

Nel=max(size(z)); %The number of electrodes.
gN=max(size(Node));
HN=max(size(Element));
M=full(sparse(zeros(gN,Nel)));
K=full(sparse(zeros(gN,gN)));
s=zeros(Nel,1);
g=reshape([Node.Coordinate],2,gN)'; %Nodes

for ii=1:HN
    % Go through all triangles
    ind=(Element(ii).Topology); % The indices to g of the ii'th triangle.
    gg=reshape([Node(ind).Coordinate],2,6)';% A 6x2 matrix of triangle nodes in (x,y) coord.
    indsig=ind(1:2:6);
    grint=grinprodgausquadnode(gg,sig(indsig));
    
    if any([Element(ii).Face{:,3}]),        %Checks if the triangle ii is the triangle that is
        % under the electrode.
        Ind=find([Element(ii).Face{:,3}]);
        a=g(Element(ii).Face{Ind,1}(1),:);
        b=g(Element(ii).Face{Ind,1}(2),:);
        c=g(Element(ii).Face{Ind,1}(3),:);
        InE=Element(ii).Face{Ind,3};         %Electrode index.
        s(InE)=s(InE)+1/z(InE)*electrlen([a;c]);% Assumes straight electrodes.
        bb1=boundquad1([a;b;c]);
        bb2=boundquad2([a;b;c]);
        
        for il=1:6
            eind=find(Element(ii).Topology(il)==Element(ii).Face{Ind,1});
            if eind
                M(ind(il),InE)=M(ind(il),InE)-1/z(InE)*bb1(eind);
            end,
            for im=1:6
                eind1=find(Element(ii).Topology(il)==Element(ii).Face{Ind,1});
                eind2=find(Element(ii).Topology(im)==Element(ii).Face{Ind,1});
                if eind1 & eind2
                    K(ind(il),ind(im))=K(ind(il),ind(im))+...
                        grint(il,im)+1/z(InE)*bb2(eind1,eind2);
                else
                    K(ind(il),ind(im))=K(ind(il),ind(im))+ ...
                        grint(il,im);
                end
            end
        end
    else %The triangle isn't under the electrode.
        for il=1:6,
            for im=1:6,
                K(ind(il),ind(im)) = K(ind(il),ind(im)) + ...
                    grint(il,im);
            end
        end
    end
end
%Calculate next the matrix S

S=full(sparse(diag(s)));
%%
C=C(:,1:Nel-1);
S=C'*S*C;
M=M*C;
A=[K,M;M',S];






