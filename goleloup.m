function [t,x]=goleloup(K,L,ic);

N=40;
if nargin==2,
    ic=rand(16,N);
%     ic=repmat([0.92 0.668 5.74 0.58 1.82 0.10 0.58 0.39 0.09 0.16 0.06 1.12 0.76 0.68 0.27 0.01]',N,1);
end

sc = 1./(1+0.05*randn(N,1));
% K=0.1;
tend=960;

options=odeset('Events','eventLeloup');
% options=[];
[t,x]=ode23(@leloupcoupled,[0 tend],ic,options,N,sc,K,L);

xm=x(:,1:16:end);
mf=mean(xm');

figure(1);
subplot(2,1,1);
plot(t,x(:,1:2*16:end));
subplot(2,1,2);
plot(t,mf);




% ==============================================================

function y=leloupcoupled(t,x,N,sc,K,L);

NV=16;   % number of variables
% N = number of oscillators
% sc = period rescaling
% K = coupling strength

%%% Model parameters:

vsp=1.5;
kap=0.7;
nn=1.9; % 4;
vmp=1.1;
kmp=0.31;
kdn=0.01;
vsc=1.1;
kac=0.6;
vmc=1.0;
kmc=0.4;
vsb=1;
kib=2.2;
mm=1.3; % 2;
vmb=0.8;
kmb=0.4;
ksp=0.6;
v1p=1.0;
kp=0.1;
v2p=0.3;
kdp=0.1;
k4=0.2;
k3=0.4;
ksc=1.6;
v1c=0.6;
v2c=0.1;
vdpc=0.7;
vdcc=0.7;
kd=0.3;
v1pc=1.0;
v2pc=0.1;
k2=0.2;
k1=0.4;
v3pc=1.0;
v4pc=0.1;
k7=0.5;
k8=0.1;
vdpcc=0.7;
vdpcn=0.7;
ksb=0.12;
v1b=0.5;
v2b=0.1;
k5=0.4;
k6=0.2;
vdbc=0.5;
v3b=0.5;
v4b=0.2;
vdbn=0.6;
vdin=0.8;
kdnc=0.12;
vphos=0.4;
vstot=1.0;

%%% Coupling parameters

vc=1;
Kc=1;


%%% Variables:

Mp=x(1:NV:end);
Mc=x(2:NV:end);
Mb=x(3:NV:end);
Pc=x(4:NV:end);
Cc=x(5:NV:end);
Pcp=x(6:NV:end);
Ccp=x(7:NV:end);
Pcc=x(8:NV:end);
Pcn=x(9:NV:end);
Pccp=x(10:NV:end);
Pcnp=x(11:NV:end);
Bc=x(12:NV:end);
Bcp=x(13:NV:end);
Bn=x(14:NV:end);
Bnp=x(15:NV:end);
In=x(16:NV:end);

%%% Coupling:

mf=mean(Mp);
% coupling=vc*K*mf/(Kc+K*mf);
% coupling=K*mf;
coupling=K*mf*((t<=192 | t>=336)*(t<=720))+...
    .6*K*mf*(t>=192 & t<=336); %+...
%     *K*mf*(t>=288 & t<=336);
L=1.25*K*(t>=720)+10*K*(t>=336 & t<=340); 

%%% Equations:

w1=vstot.*vsp.*(Bn.^nn)./(kap.^nn+Bn.^nn);
w2=vmp.*Mp./(kmp+Mp)+kdn.*Mp;
w3=vstot.*vsc.*(Bn.^nn)./(kac.^nn+Bn.^nn);
w4=vmc.*Mc./(kmc+Mc)+kdn.*Mc;
w5=vstot.*vsb.*(kib.^mm)./(kib.^mm+Bn.^mm);
w6=vmb.*Mb./(kmb+Mb)+kdn.*Mb;
w7=ksp.*Mp;
w8=v1p.*vphos.*Pc./(kp+Pc);
w9=v2p.*Pcp./(kdp+Pcp);
w10=k4.*Pcc;
w11=k3.*Pc.*Cc;
w12=kdn.*Pc;
w13=ksc.*Mc;
w14=v1c.*Cc./(kp+Cc);
w15=v2c.*Ccp./(kdp+Ccp);
w16=kdnc.*Cc;
w17=vdpc.*Pcp./(kd+Pcp)+kdn.*Pcp;
w18=vdcc.*Ccp./(kd+Ccp)+kdn.*Ccp;
w19=v1pc.*vphos.*Pcc./(kp+Pcc);
w20=v2pc.*Pccp./(kdp+Pccp);
w21=k2.*Pcn;
w22=k1.*Pcc;
w23=kdn.*Pcc;
w24=v3pc.*vphos.*Pcn./(kp+Pcn);
w25=v4pc.*Pcnp./(kdp+Pcnp);
w26=k7.*Bn.*Pcn;
w27=k8.*In;
w28=kdn.*Pcn;
w29=vdpcc.*Pccp./(kd+Pccp)+kdn.*Pccp;
w30=vdpcn.*Pcnp./(kd+Pcnp)+kdn.*Pcnp;
w31=ksb.*Mb;
w32=v1b.*Bc./(kp+Bc);
w33=v2b.*Bcp./(kdp+Bcp);
w34=k5.*Bc;
w35=k6.*Bn;
w36=kdn.*Bc;
w37=vdbc.*Bcp./(kd+Bcp)+kdn.*Bcp;
w38=v3b.*Bn./(kp+Bn);
w39=v4b.*Bnp./(kdp+Bnp);
w40=kdn.*Bn;
w41=vdbn.*Bnp./(kd+Bnp)+kdn.*Bnp;
w42=vdin.*In./(kd+In)+kdn.*In;

y(1:NV:NV*N) = sc.*(w1-w2)+coupling+L;                  % dMp/dt
y(2:NV:NV*N) = sc.*(w3-w4);                           % dMc/dt
y(3:NV:NV*N) = sc.*(w5-w6);                           % dMb/dt
y(4:NV:NV*N) = sc.*(w7-w8+w9+w10-w11-w12);            % dPc/dt
y(5:NV:NV*N) = sc.*(w13-w14+w15+w10-w11-w16);         % dCc/dt
y(6:NV:NV*N) = sc.*(w8-w9-w17);                       % dPcp/dt
y(7:NV:NV*N) = sc.*(w14-w15-w18);                     % dCcp/dt
y(8:NV:NV*N) = sc.*(-w19+w20-w10+w11+w21-w22-w23);    % dPcc/dt
y(9:NV:NV*N) = sc.*(-w24+w25-w21+w22-w26+w27-w28);    % dPcn/dt
y(10:NV:NV*N) = sc.*(w19-w20-w29);                     % dPccp/dt
y(11:NV:NV*N) = sc.*(w24-w25-w30);                     % dPcnp/dt
y(12:NV:NV*N) = sc.*(w31-w32+w33-w34+w35-w36);         % dBc/dt
y(13:NV:NV*N) = sc.*(w32-w33-w37);                     % dBcp/dt
y(14:NV:NV*N) = sc.*(-w38+w39+w34-w35-w26+w27-w40);    % dBn/dt 
y(15:NV:NV*N) = sc.*(w38-w39-w41);                     % dBnp/dt
y(16:NV:NV*N) = sc.*(-w27+w26-w42);                    % dIn/dt

%%% Output:

y=y';


