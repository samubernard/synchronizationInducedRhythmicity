function [K,Klight,N,G,H,Left,Right,VLleft,VLright,DMleft,DMright]=koupling3(R,S,dmax,options)
% KOUPLING generates a coupling matrix K
% R=geometry
%   'SCN'= 3D SCN, roughly two intersecting ellipsoids
%   'SCNslice'= 3D slice of the SCN along the VL-DM axis, 3 layer thick
%   '2Dslice'= 2D slice (1 layer thick)
% S=size
% dmax=max number of neighbors
% DF=Type of coupling K 
%    1= Nearest neighbours, equally coupled
%    2= Nearest neighbours, decreasing coupling
%    3= Uniformly random coupling, fixed connectivity
%    4= VL->DM non-symetric coupling
%    5= Diagonal coupling (self-couping) to generate a pop of dispersed
%    neurons
%    0= global coupling

globalk=0;

if nargin<4,
    DF=1;
    hext=1.01;
    hint=0.99;
    c0=0.05;
    klocal=1;
else
    DF=options(1);
    hext=options(2);
    hint=options(3);
    c0=options(4);
    klocal=options(5);
end

[G,Left,Right]=numgrid3(R,S);
neighbour=dmax;
k=find(G);
N=length(k);
S1=size(G,1);
S2=size(G,2);
S3=size(G,3);

[xi,yi,zi]=ind2sub(size(G),k);
xi=xi;
yi=yi;
zi=zi;
XI=repmat(xi,1,length(xi));
YI=repmat(yi,1,length(yi));
ZI=repmat(zi,1,length(zi));
XJ=XI';
YJ=YI';
ZJ=ZI';
pack
d=sqrt((XI-XJ).^2+(YI-YJ).^2+(ZI-ZJ).^2);
%% VL and DM parts
ctrL=[ceil(.22*S1),ceil(S2/2),ceil(0.68*S3)];
ctrR=[ceil(.78*S1),ceil(S2/2),ceil(0.68*S3)];
disttocenterL=sqrt(2*(xi-ctrL(1)).^2+(yi-ctrL(2)).^2+(zi-ctrL(3)).^2);
disttocenterR=sqrt(2*(xi-ctrR(1)).^2+(yi-ctrR(2)).^2+(zi-ctrR(3)).^2);
PhL=zeros(N,1);
PhR=zeros(N,1);
PhL=(disttocenterL>=0.4*S3/2);
PhR=(disttocenterR>=0.4*S3/2);
VLleft=zeros(size(G));
VLleft(k)=~PhL & Left(k);
VLright=zeros(size(G));
VLright(k)=~PhR & Right(k);
DMleft=zeros(size(G));
DMleft(k)=PhL & Left(k);
DMright=zeros(size(G));
DMright(k)=PhR & Right(k);

vl_left=G(find(VLleft));
vl_right=G(find(VLright));
dm_left=G(find(DMleft));
dm_right=G(find(DMright));
Klight=zeros(N,1);
Klight(vl_left)=1;
Klight(vl_right)=1;
switch DF, % Type of coupling
    case 1, %% Nearest neighbours, equally coupled
        K=d<=neighbour;
    case 2, %% Nearest neighbours, decreasing coupling
        K=(d<=neighbour)./(d+1)+2*(d==0);
    case 3, %% Uniformly random coupling, fixed connectivity
        K=rand(N)<=c0;
    case 4, %% VL->DM non-symetric coupling
        K=klocal*(d<=neighbour); % local coupling
        K(vl_left,:)=0; % no coupling for VL
        K(vl_right,:)=0;
        vl_left_proj=randsubset(vl_left,1); % projections
        vl_right_proj=randsubset(vl_right,1);
        dm_left_proj=randsubset(dm_left,length(vl_left_proj));
        dm_right_proj=randsubset(dm_right,length(vl_right_proj));
%         K=K+sparse(dm_left_proj,vl_left_proj,1,N,N);
%         K=K+sparse(dm_right_proj,vl_right_proj,1,N,N);
        K(dm_left_proj,vl_left_proj)=1;
        K(dm_right_proj,vl_right_proj)=1;
    case 5, %% Self-coupling
        K=diag(ones(1,N));
    otherwise
        disp('Using global coupling K=1\n')
        K=ones(1,N);
        globalk=1;
end

K=sparse(K);

%% Definition of heterogeneity in the SCN

mdl=[ceil(S1/2),ceil(S2/2),ceil(S3/3)];
% circular gradient
disttomiddle=sqrt((xi-mdl(1)).^2+(yi-mdl(2)).^2+(zi-mdl(3)).^2); 
% disttomiddle=zi; % linear gradient
distmax=max(disttomiddle);

H=zeros(N,1);

h1=(hext-hint)/distmax;
H=h1*disttomiddle+hint;
H2=zeros(size(G));
H2(k)=H;

figure(1);
clf;
% set(gcf,'WindowStyle','docked')
if hint~=hext,
    subplot(131)
else
    subplot(121)
end
if DF==4,
    w=sparse(repmat(vl_left,1,length(vl_left)),repmat(vl_left,1,length(vl_left))',1,N,N,length(vl_left)^2);
    spy(K.*w,[0.1 0.1 0.7])
    hold on
    w=sparse(repmat(vl_right,1,length(vl_right)),repmat(vl_right,1,length(vl_right))',1,N,N,length(vl_right)^2);
    spy(K.*w,[0.7 0.1 0.1])
    w=sparse(repmat(dm_left,1,length(dm_left)),repmat(dm_left,1,length(dm_left))',1,N,N,length(dm_left)^2);
    spy(K.*w,[0.5 0.5 0.9])
    w=sparse(repmat(dm_right,1,length(dm_right)),repmat(dm_right,1,length(dm_right))',1,N,N,length(dm_right)^2);
    spy(K.*w,[0.9 0.5 0.5])
    w=sparse(repmat(dm_left,1,length(vl_left)),repmat(vl_left,1,length(dm_left))',...
        1,N,N,length(dm_left)*length(vl_left));
    spy(K.*w,[0.3 0.3 0.3])
    w=sparse(repmat(dm_right,1,length(vl_right)),repmat(vl_right,1,length(dm_right))',...
        1,N,N,length(dm_right)*length(vl_right));
    spy(K.*w,[0.3 0.3 0.3])
else
    spy(K);  % show the coupling matrix
end
set(gca,'XTick',[])
set(gca,'YTick',[])   
title('Coupling matrix')
ylabel(['N=' num2str(N)]) 
if hint~=hext,
    subplot(132)
else
    subplot(122)
end
spy(squeeze(VLleft(:,ceil(S2/2),:))',[0.1 0.1 0.7])
hold on
spy(squeeze(VLright(:,ceil(S2/2),:))',[0.7 0.1 0.1])
spy(squeeze(DMleft(:,ceil(S2/2),:))',[0.5 0.5 0.9])
spy(squeeze(DMright(:,ceil(S2/2),:))',[0.9 0.5 0.5])
ki=find(Left & Right);
Inter=zeros(size(G));
Inter(ki)=1;
spy(squeeze(Inter(:,ceil(S2/2),:))',[0.2 0.8 0.2])
if DF==4,
    ptsl1=[xi(vl_left_proj) zi(vl_left_proj)];
    ptsl2=[xi(dm_left_proj) zi(dm_left_proj)];
    ptsr1=[xi(vl_right_proj) zi(vl_right_proj)];
    ptsr2=[xi(dm_right_proj) zi(dm_right_proj)];
    for i=1:length(vl_left_proj);
        pt1=ptsl1(i,:);
        pt2=ptsl2(i,:);
        arc(pt1,pt2,'left',[0.8 0.8 0.8])
    end
    for i=1:length(vl_right_proj);
        pt1=ptsr1(i,:);
        pt2=ptsr2(i,:);
        arc(pt1,pt2,'right',[0.8 0.8 0.8])
    end
end
axis square
set(gca,'XTick',[])
set(gca,'YTick',[])    
title('SCN slice')
if hint~=hext,
    subplot(1,3,3)
    scnplot(H2);
    title('Heterogeneity')
end