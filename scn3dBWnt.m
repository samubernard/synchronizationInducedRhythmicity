function [sol, sc, R]=scn3dBWnt(K,N,H,trans,tend,L0,light,dark,IC,scold)
% SCN3DBWNT runs the extended coupled Becker Weimann model with activation cascade for the
% neuropeptide
% To use with GOBWNT.M

%%%%%% Initial conditions

spread=0.05; % This is the relative (%) standard deviation periods of the single oscillators

if nargin<6,
    ic=repmat([0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35 ]',N,1); % initial conditions 
    sc = 1./(1+spread*randn(N,1))./H;    % period scaling
    L0=0;
    light=12;
    dark=12;
elseif nargin==8,
    ic=2*rand(10*N,1).*repmat([0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35 ]',N,1); % initial conditions 
    sc = 1./(1+spread*randn(N,1))./H;    % period scaling
else
    if isempty(IC),
        ic=repmat([0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35 ]',N,1); % initial conditions 
    else
        ic=IC;
    end
    sc=scold;
end

%%%%% Simulation
figure(2); clf;
outx=10*floor(N*rand(5,1))+1;
options=odeset('OutputS',outx,'OutputF','odeplot');
% options=[];
if (trans>0)
    [t,y] = ode45(@BWnt,[0 trans/2 trans],ic,[],sc,K,L0,light,dark);
    ic = y(end,:);
end
sol = ode45(@BWnt,[0 tend],ic,options,sc,K,L0,light,dark);        % ODE



%%%%%% Results with high resolution
resol=0.1;
xint=[0:resol:tend];
yint=deval(sol,xint);

%%%%%% Order parameter
b=sol.y(1:10:end,:);     % to compute the order param on the whole time serie
M=mean(b,1);
num=mean(M.^2)-mean(M).^2;
B=mean(b.^2,2)-mean(b,2).^2;
den=mean(B);
R=num/den;
fprintf('Order parameter: R=%g \n',R);