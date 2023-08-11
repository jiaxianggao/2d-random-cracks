function J=jacob_nodesigma(Nodesig,Elementsig,Node,Element,Agrad,U,U0,style);

% Function J=jacob_nodesigma(Node,Element,Agrad,U,U0,style);
% computes the Jacobian of the measured voltages with respect to the
% nodal values of the conductivity. If only THe Node nad Element inormation is
% given then the integral term is computed, which need to be computed
% only once for each mesh.
%
% Node=node data
% Nodesig=node data for sigma
% Element=element data
% Agrad=\int_{Element(ii) \nabla\phi_i\cdot\nabla\phi_j
% U=voltages of the injected currents
% U0=voltages of the measurement field


% M. Vauhkonen 22.10.1999, Univeristy of Kuopio, Finland


gNs=max(size(Nodesig));
gN=max(size(Node));
Hn=max(size(Elementsig));
Hn=max(size(Element));


if nargin<5
    
    Agrad=sparse(gN^2,gNs); % Gradients of the basis functions integrated over
    % each element.
    
    for jj=1:gNs
        Aa=sparse(gN,gN);
        El=Nodesig(jj).ElementConnection;
        for ii=1:max(size(El))
            indsig=Elementsig(El(ii)).Topology; % Indices of the element
            ind=Element(El(ii)).Topology; % Indices of the element
            gg=reshape([Node(ind).Coordinate],2,6)'; % A 3x2 matrix of triangle vertices in (x,y) coord.
            I=find(jj==indsig);
            if ~isempty(I)
                anis=grinprodgaus2ndnode(gg,I);
                for il=1:6
                    for im=1:6
                        Aa(ind(il),ind(im))=Aa(ind(il),ind(im))+anis(il,im);
                    end
                end
            end
        end
        Agrad(:,jj)=Aa(:);
    end
    J=Agrad;
else
    if style=='real'
        J=zeros(size(U,2)*size(U0,2),size(Agrad,2));
        for ii=1:size(Agrad,2)
            JJ=-U0.'*reshape(Agrad(:,ii),gN,gN)*U;
            JJ=JJ(:);
            J(:,ii)=JJ;
        end
    end
end








