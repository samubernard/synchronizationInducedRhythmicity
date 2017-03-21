function y=BWntcf(t,x,sc,c);
% USE ONLY WITH CONSTANTFORCING.M

% number of variables
NV=10;
N=length(x)/NV;

%%% Parameters of the single cell model:

v1b=9; 
k1b=1;
k1i=0.56;
% c=.8; %0.05
p=3; 
pc=2; 
k1d=0.18; %0.18; 
k2b=0.3;
q=2;
k2d=0.1; %0.07; 
k2t=0.36; %0.24;
k3t=0.02;
k3d=0.18; 
v4b=1; %0.3;
k4b=2.16;
r=3;
k4d=1.1; 
k5b=0.24;
k5d=0.09;
k5t=0.45;
k6t=0.06;
k6d=0.18;
k6a=0.09;
k7a=0.003;
k7d=0.13;
k8=3;       % in light: k8=5 or k8=10
k8d=4;
k9=.25;
k9d=10;
k10=1;
k10d=4;
x1t=15;
x2t=15;
%K=1.5;

%%% Variables:

y1=x(1:NV:end);
y2=x(2:NV:end);
y3=x(3:NV:end);
y4=x(4:NV:end);
y5=x(5:NV:end);
y6=x(6:NV:end);
y7=x(7:NV:end);
x1=x(8:NV:end);
x2=x(9:NV:end);
v=x(10:NV:end);

%%% Meanfield and coupling:

%mf=mean(v);
% coupling=K*x2;

%%% Coupling
%coupling=K*v;
% coupling=K*v./(1+K*v);

%%% Light

%%% Equations (FRQ):

fpercry=v1b*(y7+c+x2.^pc)./(k1b*(1+(y3/k1i).^p)+(y7+c+x2.^pc));
fbmal=v4b*y3.^r./(k4b^r+y3.^r);

y(1:NV:NV*N) = sc.*(fpercry-k1d*y1);                    % per/cry mRNA
y(2:NV:NV*N) = sc.*(k2b*y1.^q-(k2d+k2t)*y2+k3t*y3);       % PER/CRY cyto
y(3:NV:NV*N) = sc.*(k2t*y2-k3t*y3-k3d*y3);                % PER/CRY nucl
y(4:NV:NV*N) = sc.*(fbmal-k4d*y4) ;                       % bmal1 mRNA
y(5:NV:NV*N) = sc.*(k5b*y4-(k5d+k5t)*y5+k6t*y6);          % BMAL1 cyto
y(6:NV:NV*N) = sc.*(k5t*y5-(k6t+k6d)*y6+k7a*y7-k6a*y6);   % Bmal1 nucl
y(7:NV:NV*N) = sc.*(k6a*y6-(k7a+k7d)*y7);                 % Bmal1 actif
y(8:NV:NV*N) = sc.*(- k8d*x1);                            % PKA
y(9:NV:NV*N) = sc.*(k9*x1.*(x2t-x2) - k9d*x2);            % CREB
y(10:NV:NV*N) = sc.*(k10*y2 - k10d*v);                    % neurotransmitter

%%% Output:

y=y';
