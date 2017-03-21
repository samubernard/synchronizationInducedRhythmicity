% MAIN FILE TO RUN SIMULATIONS OF SCN NEURONS
% S. Bernard et al. (2007) PLoS Comput. Biol.
% samubernard@gmail.com

clc;
clear;

%%%% Parameters

geom='SCN';   % Geometry of the system
                %   'SCN'= 3D SCN, roughly two intersecting ellipsoids
                %   'SCNslice'= 3D slice of the SCN along the VL-DM axis, 3 layer thick
                %   '2Dslice'= 2D slice (1 layer thick)
D=10;        % dimension of single oscillator DO NOT CHANGE
S=13;        % System size (this is not the number of cells)
dmax=3.5;        % max distance of neighbors
trans=72;   % Transient
tend=72;    % length of simulation
K=0.9;    % coupling factor strength
Ktype=1;    % Type of coupling, see KOUPLING3 for details 
                % 1=nearest neighbours; 
                % 2=modulated nearest neighbours; 
                % 3=uniformly random.
                % 4=SCN-like
                % 5=Diagonal coupling (self-coupling)
                % 0=global
hext=1.0;   % heterogeneity, set HEXT and HINT = 1 for no heterogeneity
hint=1.0;   % heterogeneity 
c0=0.1; % nominal connectivity of random coupling
klocal=1; % local coupling strength factor (only for Ktype=4)

% Determine the period of a single neuron
[sper, scper, Rper]=scn3dBWnt(K,1,1,trans,tend,0,24,0,[],1);
resol=0.1;
xint=[0:resol:tend];
yint=deval(sper,xint);
Per=periode(xint,yint(1,:)); % a simple code to calculate the circadian period
disp(['Period of a single oscillator: ' num2str(Per)]);

%return;
%%% Light/Dark Cycle Entrainment
L0=0.0; % L0=0 no light (DD), L0=0.22 default value for light entrainment
light=12; % duration of the light phase
dark=12;  % duration of the dark phase, LIGHT+DARK does not necessarily have to add up to 24h

%%%% Task:

fprintf('Define the geometry...\n');
[Koup,Klight,N,G,H,Left,Right,VLleft,VLright,DMleft,DMright]=...
    koupling3(geom,S,dmax,[Ktype hext hint c0 klocal]); % This constructs the geometry of the system
connectivity = full((Koup>0)*ones(N,1)/N);
c=repmat(connectivity,1,N);
fprintf('Connectivity= %g \n',mean(connectivity));
KK=Koup/N./c;
KK(find(isnan(KK)))=0;  % KK is a normalized coupling matrix 
L=L0*Klight;

fprintf('Run the simulation...\n');
[sol, sc, R]=scn3dBWnt(K*KK,N,H,trans,tend,L,light,dark);

% if R<0.2,
%     fprintf('Warning, low synchronization...\n');
% %     return
% end

fprintf('And this is it! Now printing stats...\n');
stats_scn  % prints and plots various statistics and figure. The code is not robust so there are errors with some simulations

