% 2nd version, modified, difference imaging 
clear
% SetPaths
% run D:\matlabwork\eidors-v3.11-ng\eidors\startup.m

%%Define the Geomotry
squarewidth  = 18; % Width (cm)
squareheight = 4.3; % Height (cm)
L_x          = 12; % Number of horizontal eletrodes (each side)
L_y          = 2; % Number of vertical eletrodes (each side)
el_width     = .25;             % Electrode width
el_depth     = .15;              % Electrode depth
eldist_x     = [1.5*ones(1,5), 2, 1.5*ones(1,5)]; % horizontal spacing
eldist_y     = [2.3]; % vertical spacing
L = 2*(L_x + L_y);    % total number of electrodes
Nelem_el     = 1; % Number of elements under electrode

% [Element,Nodelist,H,Indb_all,Indb,E,Node,g,g1,H1,Indb1,E1,Element1,Nodelist1,Node1,Agrad] = ...
%     MakeEITMeshSquareRealEl_Cop(squarewidth,squareheight,L_x,L_y,el_width,el_depth,eldist_x,eldist_y,L,Nelem_el,1);
% 
% sN = size(g,1);

% save('D:\matlabwork\Probabilistic crack\mesh\finemesh_f_c.mat','Element','Nodelist','H','Indb_all','Indb','E','Node','g','g1','H1',...
%         'Indb1','E1','Element1','Nodelist1','Node1','Agrad',...
%         'el_width','Nelem_el','sN')  
% 
% load('D:\matlabwork\Probabilistic crack\mesh\finemesh_f_c.mat')

% save('D:\matlabwork\Probabilistic crack\mesh\coarsemesh_f_c.mat','Element','Nodelist','H','Indb_all','Indb','E','Node','g','g1','H1',...
%         'Indb1','E1','Element1','Nodelist1','Node1','Agrad',...
%         'el_width','Nelem_el','sN')  

load('D:\matlabwork\crack-liang\mesh\coarsemesh_f_c.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Target distribution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prior Gaussian background conductivity information
% to form a continuous distribution within the domain
sig0 = 10;
sigexp = sig0*ones(sN,1); % Prior expectation
sigmamax_coef = 2; % Define max conductivity (in 99,7 pros probability sense...)
sigmamax = sigmamax_coef*sig0;
sigmamin = 1; %define the minimum conductivity in the distribution, can be player around
sigma_range = [sigmamin sigmamax];
sig_var = (diff(sigma_range)/6)^2;

% isotropic smoothness with a correlating length of 4 cm 
% to incorporate spatial inhomogeneity
corrlength_x = 4; 
corrlength_y = 4;

if corrlength_x == corrlength_y % isotropic smoothness
    Gamma_pr = GaussianDistanceSmoothnessPrior(g,sig_var,corrlength_x);
else % anisotropic smoothness
    Gamma_pr = GaussianDistanceAnisotropicSmoothnessPrior(g,sig_var,corrlength_x,corrlength_y);
end

Igammasigma = inv(Gamma_pr);  %  is the inverse of the covariance matrix
L_pr_s =  chol(Igammasigma);    

% sigmaT = L_pr_s\abs(randn(sN,1))+ sigmamin*ones(sN,1);

% figure(1),clf, PlotSolScaled(g,H,sigmaT,min(sigmaT),max(sigmaT)),
% colormap(hot),
% axis image,colorbar,drawnow

% fprintf('The minimum conductivity value is: %f\n', min(sigmaT));
% fprintf('The maximum conductivity value is: %f\n', max(sigmaT));
% fprintf('The average conductivity value is: %f\n', mean(sigmaT));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate EIT measurements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = 0.171689268264545*ones(L,1);  % EXPECTED contact impedances (INVERSION)
% I = toeplitz([1; zeros(L/2-1,1); -1; zeros(L/2-1,1)]',...
%      [1 zeros(1,L/2-1)]);   % Opposite current injections
  
I = []; % All injections against electrode 1, 6, 12, 21
Cinj = [1 6 12 21];  %injection on four different eletrodes
% 4*27=108 rounds
I1 = [ones(1,L-1);-eye(L-1)];
for ipat = 1:length(Cinj)
    cel_ipat = Cinj(ipat);
    Iipat{ipat} = [I1(2:cel_ipat,:);I1(1,:);I1(cel_ipat+1:end,:)];
    I = [I,Iipat{ipat}];
end

% MeasPatt = [ones(1,L-1);-eye(L-1)];% DrawCurrents.m visualization works only with this choice
MeasPatt = -toeplitz([1;-1;zeros(L-2,1)],[1,zeros(1,L-2),-1]); % Adjacent pattern
% 108*28=3024 measurements
 
% A sample from a smoothness prior
%  blobtype = 'sharp';
 blobtype = 'smooth';
 CONDUCTIVE_BLOB = 0;
 min_sigma =  1;
 max_sigma =  10; 
 cp_x = [8];
 cp_y = [2.5];
 width_x =4;
 width_y =2.5;
 band = .2;

if strcmp(blobtype,'sharp');
    sigmaT = GenerateEllipsesInGrid(g,H,min_sigma,max_sigma,width_x,width_y,...
    cp_x,cp_y,band,CONDUCTIVE_BLOB);
elseif strcmp(blobtype,'smooth');
    sigmaT = GenerateAnistropicBubblesInGrid(g,H,min_sigma,max_sigma,width_x,width_y,...
    cp_x,cp_y,1,CONDUCTIVE_BLOB);
end

%%%%% Generate training data %%%%
N = 10; % number of training samples
sigTdata = [];
U_T = [];
Uel_T = [];
sHom = 5;

crack_loc1 = [];
crack_loc2 = [];

for j = 1:N
    %%% A sample from a smoothness prior
    sigmaT(:,j) = L_pr_s\abs(randn(sN,1))+ sigmamin*ones(sN,1);
    NoCrack = 2; % number of cracks
    crack_loc1 = [];
    crack_loc2 = [];

    if NoCrack == 1
        %%% Crack Parameter
        c_length =randi([50,300],1); %The Total Drunkstep the crack takes
        c1 =(0.14-0.11) + 0.01; %drunkstep length the crack takes in y-direction
        c2 =(0.14-0.11) + 0.01; %drunkstep length the crack takes in x-direction
        %%% crack 1
        a = 0;
        b = 3.625;
        c_startx = (b-a).*rand + a; % random number between a and b
        csx = dsearchn(g(:,1),c_startx); %start node idx for crack, dsearch returns index
        g_start1 = [g(csx,1),0]; % coordinate/location of start node

        for i = 1:c_length
            if i == 1
                crack_loc1(i,:) = g_start1;
            else
                drunkstepy1 = c1*rand;
                drunkstepx1 = c2*rand;
                drunkstep1 = [drunkstepx1 drunkstepy1];
                drunkstep1 = crack_loc1(i-1,:) + drunkstep1;
                crack_loc1(i,:) = drunkstep1;
                if crack_loc1(i,:) > squarewidth | crack_loc1(i,:) > squareheight
                    break
                end
            end
        end

        crackidx1 = dsearchn(g,crack_loc1);
        crackx1 = g(crackidx1,1);
        cracky1 = g(crackidx1,2);
        crack1 = [crackx1,cracky1];

        sigmaT_crack(:,j) = sigmaT(:,j);
        sigmaT_crack(crackidx1,j)= 1e-6;

    else if NoCrack == 2
            c_length =randi([50,300],1); %The Total Drunkstep the crack takes

            c1 = 0.04; %drunkstep length the crack takes in y-direction
            c2 = 0.04; %drunkstep length the crack takes in x-direction
            %%% crack 1
            a = 3.625;
            b = 14.375;
            c_startx = (b-a).*rand + a;
            csx = dsearchn(g(:,1),c_startx); %start node idx for crack
            g_start1 = [g(csx,1),0];

            for i = 1:c_length
                if i == 1
                    crack_loc1(i,:) = g_start1;
                else
                    drunkstepy1 = c1*rand;
                    drunkstepx1 = c2*randn;
                    drunkstep1 = [drunkstepx1 drunkstepy1];
                    drunkstep1 = crack_loc1(i-1,:) + drunkstep1;
                    crack_loc1(i,:) = drunkstep1;
                    if crack_loc1(i,:) > squarewidth | crack_loc1(i,:) > squareheight
                        break
                    end
                end
            end

            crackidx1 = dsearchn(g,crack_loc1);
            crackx1 = g(crackidx1,1);
            cracky1 = g(crackidx1,2);
            crack1 = [crackx1,cracky1];

            sigmaT_crack(:,j) = sigmaT(:,j);
            sigmaT_crack(crackidx1,j)= 1e-6;
            %%%crack2
            c_length =randi([50,300],1); %The Total Drunkstep the crack takes
            c1 = 0.04; %drunkstep length the crack takes in y-direction
            c2 = 0.04; %drunkstep length the crack takes in x-direction
            a = 3.625;
            b = 14.375;
            c_startx = (b-a).*rand + a;
            csx = dsearchn(g(:,1),c_startx);
            g_start2 = [g(csx,1),0];

            for i = 1:c_length
                if i == 1
                    crack_loc2(i,:) = g_start2;
                else
                    drunkstepy2 = c1*rand;
                    drunkstepx2 = c2*randn;
                    drunkstep2 = [drunkstepx2 drunkstepy2];
                    drunkstep2 = crack_loc2(i-1,:) + drunkstep2;
                    crack_loc2(i,:) = drunkstep2;
                    if crack_loc2(i,:) > squarewidth | crack_loc2(i,:) > squareheight
                        break
                    end
                end
            end

            crackidx2 = dsearchn(g,crack_loc2);
            crackx2 = g(crackidx2,1);
            cracky2 = g(crackidx2,2);
            crack2 = [crackx2,cracky2];

            sigmaT_crack(crackidx2,j)= 1e-6;

    end
    end
end
% figure(2),clf, PlotSolScaled(g,H,sigmaT_crack,min(sigmaT_crack),max(sigmaT_crack)),
% colormap(hot),
% axis image,colorbar,drawnow
% title('true conductivity')
Uel_Y=[];



parfor j = 1:N     
    
   %%%  Add random noise onto the traning data
    aR = 1e-3; %noise pararmeters of the traning samples
    bR = 1e-4;
    aR2 = 1e-2;
    bR2 = 1e-3;
    
    meas_noise_coef = (bR-aR).*rand + aR;  % (meas_noise_coef*(max(Uel)-min(Uel)))^2
    meas_noise_coef2 = (bR2-aR2).*rand + aR2; %  var(l) = (meas_noise_coef2*Uel(l))^2
    
    UT = ForwardSolution2ndNode(Node1,Element1,I,sigmaT_crack(:,j),z,MeasPatt,'real');
    Uel_d = UT.Electrode(:);
    Uel_Y(:,j)=Uel_d;
    Uel_nonoise = Uel_Y(:,j);

    %     %%% Add solid noise onto the training data
    %     meas_noise_coef = 1e-3;  % (meas_noise_coef*(max(Uel)-min(Uel)))^2
    %     meas_noise_coef2 = 1e-2; %  var(l) = (meas_noise_coef2*Uel(l))^2
    % add noise to measurements
    Uel1 = Uel_nonoise + meas_noise_coef2.*abs(Uel_nonoise).*randn(size(Uel_nonoise));
    Uel = Uel1 + (meas_noise_coef.*(max(Uel1)-min(Uel1))).*randn(size(Uel1))
    
    sigTdata(:,j) = sigmaT_crack(:,j); %%% conductivity data
    U_T(:,j) = Uel; %%% voltage training data
end

  
%%%%Generating Background Conductivity Distribution%%

parfor j = 1:N   
       %%%  Add random noise onto the traning data
    aR = 1e-3; %noise pararmeters of the traning samples
    bR = 1e-4;
    aR2 = 1e-2;
    bR2 = 1e-3;
    
    meas_noise_coef = (bR-aR).*rand + aR;  % (meas_noise_coef*(max(Uel)-min(Uel)))^2
    meas_noise_coef2 = (bR2-aR2).*rand + aR2; %  var(l) = (meas_noise_coef2*Uel(l))^2
    
    UT_I = ForwardSolution2ndNode(Node1,Element1,I,sigmaT(:,j),z,MeasPatt,'real');
    Uel_T_I = UT_I.Electrode(:);
    Uel_nonoise_I = Uel_T_I;
    
    Uel1_I = Uel_nonoise_I + meas_noise_coef2*abs(Uel_nonoise_I).*randn(size(Uel_nonoise_I));
    Uel_I = Uel1_I + (meas_noise_coef*(max(Uel1_I)-min(Uel1_I)))*randn(size(Uel1_I));
    
    sigTdata_I(:,j) = sigmaT(:,j);
    U_T_I(:,j) = Uel_I;
    
end

RRRr_Vol2 = U_T - U_T_I;
RRRr_Con2 = sigTdata - sigTdata_I;
% RRRr_Con2 = (RRRr_Con2 ~= 0);
Split = 0.8;

XTrain = RRRr_Vol2(:,1:Split*N);
YTrain = RRRr_Con2(:,1:Split*N);
XValidation = RRRr_Vol2(:,Split*N+1:N);
YValidation = RRRr_Con2(:,Split*N+1:N);   


save('Diff_Flextural_Crack_Train.mat','XTrain','YTrain')
save('Diff_Flextural_Crack_Test.mat','XValidation','YValidation')
