
function [Element,Nodelist,H,Indb,IndbEl,E,Node,g,g1,H1,IndbEl1,E1,Element1,Nodelist1,Node1,Agrad] = ...
         MakeEITMeshSquareRealEl_Cop(squarewidth,squareheight,L_x,L_y,el_width,el_depth,eldist_x,eldist_y,L,Nelem_el,DRAW)
     
%keyboard

%%%%%%%%%%%%%%%%%%%%%
%% Mesh parameters %%
%%%%%%%%%%%%%%%%%%%%%

%     squarewidth  = 18;
%     squareheight = 4;
%     L_x          = 12;
%     L_y          = 4;
%     el_width     = .75;             % Electrode width
%     el_depth     = .5;              % Electrode depth
%     eldist_x     = 1.4;
%     eldist_y     = 1;
%     L = 2*(L_x + L_y);    % Number of electrodes
%     Nelem_el     = 2;     % Number of elements under electrode


%el = [Nelem_el,Nelem_bel];
%minelemsize = el_width/Nelem_el;

[g,H,E,Indb,IndbEl] = squaregrid_eit_realel_Cop(squarewidth,squareheight,L_x,L_y,el_width,el_depth,eldist_x,eldist_y,L,Nelem_el);
%[g,H,E,Indb,IndbEl] = squaregrid_eit_realel_RAND(squarewidth,squareheight,L_x,L_y,el_width,el_depth,eldist_x,eldist_y,L,Nelem_el);

if DRAW
    figure(66),clf,trimesh(H,g(:,1),g(:,2),zeros(length(g),1),'FaceColor',[1 1 1],'EdgeColor',[1 0 0]);
    view(2),axis image;
    drawnow
end
%keyboard

% check
eltest = zeros(size(H,1),1);
for jj=1:L
  eltest(E(jj,:)) = jj;
end

%keyboard 
% figure(5),clf, markgrid(eltest,g,H), axis image

[Element,Nodelist]=MakeElement(H,IndbEl,E);
[Node]=MakeNode(Element,Nodelist,g);

[g1,H1,IndbEl1,E1,Element1,Nodelist1,Node1] = ModifyMesh2nd(g,H,IndbEl,E);

Agrad=jacob_nodesigma2nd(Node,Element,Node1,Element1);

