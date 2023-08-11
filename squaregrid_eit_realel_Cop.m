function [gnew,Hnew,E,Indb,IndbEl] = sqaregrid_eit_realel_Cop(squarewidth,squareheight,L_x,L_y,el_width,el_depth,eldist_x,eldist_y,L,Nelem_el,bel_width_x,bel_width_y,side_dist_x,side_dist_y,Nelem_bel_x,Nelem_bel_y,Nelem_side_x,Nelem_side_y);

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

%keyboard

bel_width_x = eldist_x - el_width;
bel_width_y = eldist_y - el_width;
side_dist_x = (squarewidth - sum(eldist_x) - el_width)/2;
side_dist_y = (squareheight - sum(eldist_y) - el_width)/2;
side_width1_x = el_depth;
side_width2_x = side_dist_x - side_width1_x;
side_width1_y = el_depth;
side_width2_y = side_dist_y - side_width1_y;


% bel_width_x = eldist_x - el_width;
% bel_width_y = eldist_y - el_width;
% side_dist_x = (squarewidth - L_x*el_width - (L_x-1)*bel_width_x)/2;
% side_dist_y = (squareheight - L_y*el_width - (L_y-1)*bel_width_y)/2;
% side_width1_x = el_depth;
% side_width2_x = side_dist_x - side_width1_x;
% side_width1_y = el_depth;
% side_width2_y = side_dist_y - side_width1_y;

Nelem_bel_x = round(bel_width_x/(el_width/Nelem_el));
Nelem_bel_y = round(bel_width_y/(el_width/Nelem_el));
Nelem_side1_x = round(side_width1_x/(el_width/Nelem_el));
Nelem_side2_x = round((side_width2_x)/(el_width/Nelem_el));
Nelem_side1_y = round(side_width1_y/(el_width/Nelem_el));
Nelem_side2_y = round((side_width2_y)/(el_width/Nelem_el));

%%%%%%%%%%%%%

% element sizes
el_elwidth = el_width/Nelem_el;
el_belwidth_x = bel_width_x./Nelem_bel_x;
el_belwidth_y = bel_width_y./Nelem_bel_y;
%el_sidewidth_x = side_dist_x/Nelem_side_x;
%el_sidewidth_y = side_dist_y/Nelem_side_y;
el_sidewidth1_x = side_width1_x/Nelem_side1_x;
el_sidewidth2_x = side_width2_x/Nelem_side2_x;
el_sidewidth1_y = side_width1_y/Nelem_side1_y;
el_sidewidth2_y = side_width2_y/Nelem_side2_y;

%keyboard

x = linspace(0,side_width1_x,Nelem_side1_x+1);
xnew = linspace(max(x)+el_sidewidth2_x,side_width1_x+side_width2_x,Nelem_side2_x);
x = [x, xnew];
for iel = 1:(L_x-1)
    xnew = linspace(max(x)+el_elwidth,max(x)+el_width,Nelem_el);
    x = [x, xnew];
    xnew = linspace(max(x)+el_belwidth_x(iel),max(x)+bel_width_x(iel),Nelem_bel_x(iel));
    x = [x, xnew];
end
xnew = linspace(max(x)+el_elwidth,max(x)+el_width,Nelem_el);
x = [x, xnew];
xnew = linspace(max(x)+el_sidewidth2_x,max(x)+side_width2_x,Nelem_side2_x);
x = [x, xnew];
xnew = linspace(max(x)+el_sidewidth1_x,squarewidth,Nelem_side1_x);
x = [x, xnew];

y = linspace(0,side_width1_y,Nelem_side1_y+1);
ynew = linspace(max(y)+el_sidewidth2_y,side_width1_y+side_width2_y,Nelem_side2_y);
y = [y, ynew];
for iel = 1:(L_y-1)
    ynew = linspace(max(y)+el_elwidth,max(y)+el_width,Nelem_el);
    y = [y, ynew];
    ynew = linspace(max(y)+el_belwidth_y(iel),max(y)+bel_width_y(iel),Nelem_bel_y(iel));
    y = [y, ynew];
end
ynew = linspace(max(y)+el_elwidth,max(y)+el_width,Nelem_el);
y = [y, ynew];
ynew = linspace(max(y)+el_sidewidth2_y,max(y)+side_width2_y,Nelem_side2_y);
y = [y, ynew];
ynew = linspace(max(y)+el_sidewidth1_y,squareheight,Nelem_side1_y);
y = [y, ynew];

% figure,plot(x, ones(length(x),1),'*')
% figure,plot(y, ones(length(y),1),'*')
[X,Y] = meshgrid(x,y);
g = [X(:),Y(:)];
%figure,plot(g(:,1),g(:,2),'+'),axis image
H = delaunay(g(:,1),g(:,2));
% figure(7),clf,trimesh(H,g(:,1),g(:,2),zeros(size(g,1),1)),
% view(2), axis image

HlineDelete = [];
gx = g(:,1);
gy = g(:,2);
elementcenters = [mean(gx(H),2),mean(gy(H),2)];

%keyboard

y1 = -el_depth;
y2 = el_depth - .01*el_elwidth;
y1B = squareheight - el_depth + .01*el_elwidth;
y2B = squareheight + .01*el_elwidth;
for iel = 1:L_x
    
   %x1 = side_width1_x + side_width2_x + (iel-1)*(el_width+bel_width_x) + .1*el_elwidth;
   %x2 = side_width1_x + side_width2_x + (iel-1)*(el_width+bel_width_x) + el_width - .1*el_elwidth;
   x1 = side_width1_x + side_width2_x + (iel-1)*(el_width)+sum(bel_width_x(1:iel-1))  + .1*el_elwidth;
   x2 = side_width1_x + side_width2_x + (iel-1)*(el_width)+sum(bel_width_x(1:iel-1)) + el_width - .1*el_elwidth;
   und_el_ind_x{2*L_x+1-iel} = find( g(:,1) > x1 & g(:,1) < x2 & g(:,2) > y1 & g(:,2) < y2);
   und_el_ind_x{iel} = find( g(:,1) > x1 & g(:,1) < x2 & g(:,2) > y1B & g(:,2) < y2B);
   und_el_Hind1 = find( elementcenters(:,1) > x1 & elementcenters(:,1) < x2 & elementcenters(:,2) > y1 & elementcenters(:,2) < y2);
   und_el_Hind2 = find( elementcenters(:,1) > x1 & elementcenters(:,1) < x2 & elementcenters(:,2) > y1B & elementcenters(:,2) < y2B);
   HlineDelete = [HlineDelete;und_el_Hind1;und_el_Hind2];
   elpts_x1 = H(und_el_Hind1,:);
   elpts_x1 = elpts_x1(:);
   elpts_x1 = setdiff(elpts_x1,und_el_ind_x{2*L_x+1-iel});
   elpts_x{2*L_x+1-iel} = elpts_x1;
   elpts_x2 = H(und_el_Hind2,:);
   elpts_x2 = elpts_x2(:);
   elpts_x2 = setdiff(elpts_x2,und_el_ind_x{iel});
   elpts_x{iel} = elpts_x2;
end

% figure(100), hold on,
% for iel = 1:L_x*2
%     plot(g(und_el_ind_x{iel},1),g(und_el_ind_x{iel},2),'bo')
%     pause
% end

x1 = -el_depth;
x2 = el_depth - .01*el_elwidth;
x1B = squarewidth - el_depth + .01*el_elwidth;
x2B = squarewidth + .01*el_elwidth;
for iel = 1:L_y    
   y1 = side_width1_y + side_width2_y + (iel-1)*(el_width)+sum(bel_width_y(1:iel-1)) + .1*el_elwidth;
   y2 = side_width1_y + side_width2_y + (iel-1)*(el_width)+sum(bel_width_y(1:iel-1)) + el_width - .1*el_elwidth;
   und_el_ind_y{L_y+iel} = find( g(:,2) > y1 & g(:,2) < y2 & g(:,1) > x1 & g(:,1) < x2);
   und_el_ind_y{L_y+1-iel} = find( g(:,2) > y1 & g(:,2) < y2 & g(:,1) > x1B & g(:,1) < x2B);
   und_el_Hind1 = find( elementcenters(:,2) > y1 & elementcenters(:,2) < y2 & elementcenters(:,1) > x1 & elementcenters(:,1) < x2);
   und_el_Hind2 = find( elementcenters(:,2) > y1 & elementcenters(:,2) < y2 & elementcenters(:,1) > x1B & elementcenters(:,1) < x2B);
   HlineDelete = [HlineDelete;und_el_Hind1;und_el_Hind2];
   elpts_y1 = H(und_el_Hind1,:);
   elpts_y1 = elpts_y1(:);
   elpts_y1 = setdiff(elpts_y1,und_el_ind_y{L_y+iel});
   elpts_y{L_y+iel} = elpts_y1;
   elpts_y2 = H(und_el_Hind2,:);
   elpts_y2 = elpts_y2(:);
   elpts_y2 = setdiff(elpts_y2,und_el_ind_y{L_y+1-iel});
   elpts_y{L_y+1-iel} = elpts_y2;
end

% figure(100), hold on,
% for iel = 1:L_y*2
%     plot(g(und_el_ind_y{iel},1),g(und_el_ind_y{iel},2),'ro')
%     pause
% end

for iel = 1:L_x 
    und_el_ind{iel} = und_el_ind_x{iel};
    und_el_ind{L_x+L_y+iel} = und_el_ind_x{L_x+iel};
    elpts{iel} = elpts_x{iel};
    elpts{L_x+L_y+iel} = elpts_x{L_x+iel};
end
for iel = 1:L_y
    und_el_ind{L_x+iel} = und_el_ind_y{iel};
    und_el_ind{2*L_x+L_y+iel} = und_el_ind_y{L_y+iel};
    elpts{L_x+iel} = elpts_y{iel};
    elpts{2*L_x+L_y+iel} = elpts_y{L_y+iel};
end

% figure(100),clf,trimesh(H,g(:,1),g(:,2),zeros(size(g,1),1)),
% view(2), axis image, hold on
% for iel = 1:L
%     plot(g(und_el_ind{iel},1),g(und_el_ind{iel},2),'ko')
%     pause
% end

% figure(100),clf,trimesh(H,g(:,1),g(:,2),zeros(size(g,1),1)), 
% view(2), axis image, hold on
% for iel = 1:L
%     plot(g(elpts{iel},1),g(elpts{iel},2),'ko')
%     pause
% end

glineDelete = [];
IndbEl = [];
for iel = 1:L
    glineDelete = [glineDelete;und_el_ind{iel}];
end

gnew = g;
Hnew = H;
gnew(glineDelete,:) = [];
Hnew(HlineDelete,:) = [];

%keyboard

[newind] = findnearestpoint(gnew,g);
Hnew = newind(Hnew);
%  figure(100),clf,trimesh(Hnew,gnew(:,1),gnew(:,2),zeros(size(gnew,1),1)),
%  view(2),axis image

% E = [];
% for iel = 1:L
%     newelpts{iel} = findnearestpoint(gnew,g(elpts{iel},:));
%     ee = [];
%     for iH = 1:size(Hnew,1)
%         fnd = intersect(Hnew(iH,:),newelpts{iel});
%         if length(fnd) == 2
%             ee = [ee, iH];
%         end
%     end
%     E(iel,:) = ee;
% end

E = [];
for iel = 1:L
    newelpts{iel} = findnearestpoint(gnew,g(elpts{iel},:));
    pt = newelpts{iel};
    Htmp = zeros(size(Hnew,1),1);
    for ipts = 1:length(pt)
        [ii,jj] = find(Hnew == pt(ipts));
        Htmp(ii) = Htmp(ii) + ones(length(ii),1);
    end
    E(iel,:) = find(Htmp == 2);
end
   
% figure(100),clf,trimesh(Hnew,gnew(:,1),gnew(:,2),zeros(size(gnew,1),1)),
% view(2),axis image, hold on
% for iel = 1:L
%     plot(gnew(newelpts{iel},1),gnew(newelpts{iel},2),'ko')
%     pause
% end

IndbEl = [];
for iel = 1:L
    IndbEl = [IndbEl;newelpts{iel}];
end

Indbout = find(gnew(:,1) < .01*el_elwidth | gnew(:,1) > squarewidth - .01*el_elwidth | ...
               gnew(:,2) < .01*el_elwidth | gnew(:,2) > squareheight - .01*el_elwidth);

Indb = unique([Indbout; IndbEl]);

% figure(100),clf,trimesh(Hnew,gnew(:,1),gnew(:,2),zeros(size(gnew,1),1)),
% view(2),axis image, hold on
% plot(gnew(Indb,1),gnew(Indb,2),'ko')

% eltest = zeros(size(Hnew,1),1);
% for jj=1:L
%   eltest(E(jj,:)) = jj;
% end
% figure, markgrid(eltest,gnew,Hnew), %axis image






