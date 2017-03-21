% stats_scn script to compute statistics of the coupled system


%% GEOMETRY

% N=size(K,2);
lscn=G(find(Left));
rscn=G(find(Right));
vlleft=G(find(VLleft));
vlright=G(find(VLright));
dmleft=G(find(DMleft));
dmright=G(find(DMright));
if strcmp(geom,'SCN'),
    dim=3;
else
    dim=2;
end
colormap jet

%% Intrinsic periods

figure(3); clf;
% set(gcf,'WindowStyle','docked')
subplot(221)
xPer=linspace(min(Per./sc),max(Per./sc),10);
hvll=hist(Per./sc(vlleft),xPer);
hvlr=hist(Per./sc(vlright),xPer);
hdml=hist(Per./sc(dmleft),xPer);
hdmr=hist(Per./sc(dmright),xPer);
bar(xPer,[hvll' hvlr' hdml' hdmr'],1.5)
title('Intrinsic periods of the neurons')
xlabel(['Number of neurons N=' num2str(N)])
set(gca,'Xtick',16:4:32);
subplot(222)
C=zeros(size(G));
k=find(G);
C(k)=Per./sc;
scnplot(C,dim);
title('Intrinsic periods of the neurons')

%% connectivity

connectivity = (K>0)*ones(N,1)/N;
subplot(223)
if max(connectivity)>min(connectivity);
    xcon=linspace(min(connectivity),max(connectivity),10);
    hvll=hist(connectivity(vlleft),xcon);
    hvlr=hist(connectivity(vlright),xcon);
    hdml=hist(connectivity(dmleft),xcon);
    hdmr=hist(connectivity(dmright),xcon);
    bar(xcon,[hvll' hvlr' hdml' hdmr'],1.5)
    title('Connectivity of the neurons')
    xlabel(['Number of neurons N=' num2str(N)])
    subplot(224)
    C=zeros(size(G));
    k=find(G);
    C(k)=connectivity;
    scnplot(C,dim);
    title('Connectivity of the neurons')
end


%% DYNAMICS

resol=0.1;
xint=[0:resol:tend];
yint=deval(sol,xint);
if N>1,
    mf=mean(yint(1:D:end,:));
else
    mf=yint(1:D:end,:);
end
    
%% Time evolution

figure(4); clf;
% set(gcf,'WindowStyle','docked')
subplot(221)
% plot(xint,yint(10*D+1,:),xint,yint(11*D,:),xint,mf);
% title('Y1, V, and MF')
plot(xint,mf)
title('Global Mean Field')
set(gca,'XTick',0:48:tend)
subplot(222)
plot(xint,mean(yint(D*lscn,:)),'b')
hold on
plot(xint,mean(yint(D*rscn,:)),'r')
title('Left and Right Mean Fields')
legend('L','R')
set(gca,'XTick',0:48:tend)
subplot(223)
outv=ceil(length(lscn)*rand(10,1));
plot(xint,yint(D*lscn(outv),:))
title('VI Left SCN (10 randomly selected)')
set(gca,'XTick',0:48:tend)
subplot(224)
outv=ceil(length(rscn)*rand(10,1));
plot(xint,yint(D*rscn(outv),:))
title('VI Right SCN (10 randomly selected)')
set(gca,'XTick',0:48:tend)

%% Amplitude and Phase

amplitude = max(yint,[],2)-min(yint,[],2);
figure(5); clf;
% set(gcf,'WindowStyle','docked')
subplot(221)
xamp=linspace(0,max(amplitude(1:D:end)),10);
hvll=hist(amplitude(D*vlleft-D+1)',xamp);
hvlr=hist(amplitude(D*vlright-D+1)',xamp);
hdml=hist(amplitude(D*dmleft-D+1)',xamp);
hdmr=hist(amplitude(D*dmright-D+1)',xamp);
bar(xamp,[hvll' hvlr' hdml' hdmr'],1.5)
title('Amplitude of the neurons in X')
xlabel(['Number of neurons N=' num2str(N)])
subplot(222)
C=zeros(size(G));
k=find(G);
C(k)=amplitude(1:D:end);
scnplot(C,dim);
title('Amplitude of the neurons in X')

%% phase

% maxmf=find(mf>[mf(2:end) Inf] & mf>[Inf mf(1:end-1)]);
% tmaxmf=xint(maxmf);
% lmmf=length(tmaxmf);
% for i=1:N,
%     maxi=find(yint(D*i,:)>[yint(D*i,2:end) Inf] & yint(D*i,:)>[Inf yint(D*i,1:end-1)]);
%     tmaxi=xint(maxi);
%     lpd=min(lmmf,length(tmaxi));
%     phasediff(i)=mod(mean(tmaxmf(1:lpd)-tmaxi(1:lpd))+12,24)-12;
% end
minmf=find(mf(:)'<[mf(2:end) Inf] & mf(:)'<[Inf mf(1:end-1)]);
tminmf=xint(minmf);
lmmf=length(tminmf);
for i=1:N,
    mini=find(yint(D*(i-1)+1,:)<[yint(D*(i-1)+1,2:end) Inf] & yint(D*(i-1)+1,:)<[Inf yint(D*(i-1)+1,1:end-1)]);
    tmini=xint(mini);
    lpd=min(lmmf,length(tmini));
    phasediff(i)=mod(mean(tminmf(1:lpd)-tmini(1:lpd))+12,24)-12;
end

subplot(223)
xpd=linspace(min(phasediff),max(phasediff),10);
hvll=hist(phasediff(vlleft)',xpd);
hvlr=hist(phasediff(vlright)',xpd);
hdml=hist(phasediff(dmleft)',xpd);
hdmr=hist(phasediff(dmright)',xpd);
bar(xpd,[hvll' hvlr' hdml' hdmr'],1.5)
title('Phase difference of the neurons in X')
xlabel(['Number of neurons N=' num2str(N)])
subplot(224)
C=zeros(size(G));
k=find(G);
C(k)=phasediff;
scnplot(C,dim);
title('Phase difference of the neurons in X')

figure(6); clf;
% set(gcf,'WindowStyle','docked')
subplot(121)
% plot(Per./sc(lscn),phasediff(lscn),'.b')
% hold on;
% plot(Per./sc(rscn),phasediff(rscn),'.r')
plot(Per./sc(vlleft),phasediff(vlleft),'.','Color',[0.1 0.1 0.7])
hold on;
plot(Per./sc(vlright),phasediff(vlright),'.','Color',[0.7 0.1 0.1])
plot(Per./sc(dmleft),phasediff(dmleft),'.','Color',[0.5 0.5 0.9])
plot(Per./sc(dmright),phasediff(dmright),'.','Color',[0.9 0.5 0.5])
line=[Per-6 0; Per+6 0];
plot(line(:,1),line(:,2),'k--');
hold on;
line=[Per -Per/2; Per Per/2];
plot(line(:,1),line(:,2),'k--');
axis tight
xlabel('Period (h)');
ylabel('Phase difference with respect to the mean field (h)');
title('Phase difference vs intrinsic period');
subplot(122)
C=zeros(size(G));
k=find(G);
C(k)=phasediff;
scnplot(C,dim);

%% Period after coupling

figure(7); clf;
% set(gcf,'WindowStyle','docked')
for i=1:N,
    per(i)=periode(xint,yint(D*(i-1)+1,:));
end
subplot(221)
hist(per)
xper=linspace(min(per),min([max(per) 35]),10);
hvll=hist(per(vlleft)',xper);
hvlr=hist(per(vlright)',xper);
hdml=hist(per(dmleft)',xper);
hdmr=hist(per(dmright)',xper);
bar(xper,[hvll' hvlr' hdml' hdmr'],1.5)
title('Periods of the neurons')
xlabel(['Number of neurons N=' num2str(N)])
subplot(222)
C=zeros(size(G));
k=find(G);
C(k)=per;
scnplot(C,dim);

%% Phase difference

subplot(223)
hvll=hist(phasediff(vlleft)',xpd);
hvlr=hist(phasediff(vlright)',xpd);
hdml=hist(phasediff(dmleft)',xpd);
hdmr=hist(phasediff(dmright)',xpd);
bar(xpd,[hvll' hvlr' hdml' hdmr'],1.5)
title('Phase difference of the neurons in X')
xlabel(['Number of neurons N=' num2str(N)])
subplot(224)
C=zeros(size(G));
k=find(G);
C(k)=phasediff;
scnplot(C,dim);
title('Phase difference of the neurons in V')

% Amplitude vs phase and period

% figure(8); clf;
% % set(gcf,'WindowStyle','docked')
% subplot(211)
% plot(phasediff(vlleft),amplitude(D*vlleft-D+1),'.','Color',[0.1 0.1 0.7])
% hold on;
% plot(phasediff(vlright),amplitude(D*vlright-D+1),'.','Color',[0.7 0.1 0.1])
% plot(phasediff(dmleft),amplitude(D*dmleft-D+1),'.','Color',[0.5 0.5 0.9])
% plot(phasediff(dmright),amplitude(D*dmright-D+1),'.','Color',[0.9 0.5 0.5])
% title('Amplitude vs phase relationship')
% ylabel('Amplitude in X')
% xlabel('Phase (h)')
% subplot(212)
% plot(Per./sc(vlleft),amplitude(D*vlleft-D+1),'.','Color',[0.1 0.1 0.7])
% hold on;
% plot(Per./sc(vlright),amplitude(D*vlright-D+1),'.','Color',[0.7 0.1 0.1])
% plot(Per./sc(dmleft),amplitude(D*dmleft-D+1),'.','Color',[0.5 0.5 0.9])
% plot(Per./sc(dmright),amplitude(D*dmright-D+1),'.','Color',[0.9 0.5 0.5])
% title('Amplitude vs period relationship')
% ylabel('Amplitude in X')
% xlabel('period (h)')
figure(8); clf;
% set(gcf,'WindowStyle','docked')
subplot(211)
plot(phasediff(vlleft),amplitude(D*vlleft-D+1),'.','Color',[0.1 0.1 0.7])
hold on;
plot(phasediff(vlright),amplitude(D*vlright-D+1),'.','Color',[0.7 0.1 0.1])
plot(phasediff(dmleft),amplitude(D*dmleft-D+1),'.','Color',[0.5 0.5 0.9])
plot(phasediff(dmright),amplitude(D*dmright-D+1),'.','Color',[0.9 0.5 0.5])
title('Amplitude vs phase relationship')
ylabel('Amplitude in X')
xlabel('Phase (h)')
subplot(212)
plot(Per./sc(vlleft),amplitude(D*vlleft-D+1),'.','Color',[0.1 0.1 0.7])
hold on;
plot(Per./sc(vlright),amplitude(D*vlright-D+1),'.','Color',[0.7 0.1 0.1])
plot(Per./sc(dmleft),amplitude(D*dmleft-D+1),'.','Color',[0.5 0.5 0.9])
plot(Per./sc(dmright),amplitude(D*dmright-D+1),'.','Color',[0.9 0.5 0.5])
title('Amplitude vs period relationship')
ylabel('Amplitude in X')
xlabel('period (h)')

figure(9); clf;
set(gcf,'Position',[100 400 1000 125])
fr=[1:48]';
pos=[0.1+0.033*(mod(fr-1,24)) 0.55*(fr<=24)+0.1*(fr>24) 0.03+0*fr 0.3+0*fr];
xcc=1:1:48;
ycc=deval(sol,xcc);
xp=1;
yccs=ycc(xp:D:end,:)/max(max(ycc(xp:D:end,:),[],2)-min(ycc(xp:D:end,:),[],2));
yccs=yccs-min(min(yccs));
% yccs=ycc(xp:D:end,:)./repmat(max(ycc(xp:D:end,:),[],2)-min(ycc(xp:D:end,:),[],2),1,size(ycc,2));
% yccs=yccs-repmat(min(yccs,[],2),1,size(yccs,2));,
lx=length(xcc);
for i=1:48,
    subplot('position',pos(i,:))
    C=zeros(size(G));
    k=find(G);
    C(k)=yccs(:,i);
    scnplot2(C,2,2);
    ylim=get(gca,'YLim');
    text(0,ylim(1)-0.15*(ylim(2)-ylim(1)),[num2str(xcc(i)) ' h'],'FontS',7)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    %axis square
end


%% Time evolution

figure(10); clf;
% set(gcf,'WindowStyle','docked')
xsurf=linspace(0,tend,100);
ysurf=deval(xsurf,sol);
skp=D*ceil(N/50);
%plotNeuron=1:ceil(N/50):N;
plotNeuron=1:N;
[region lreg] = findregion(plotNeuron,vlleft,vlright,dmleft,dmright);
[sortReg, regIndex]=sort(region);
sortNeuron=plotNeuron(regIndex);
[xm, ym]=meshgrid(xsurf,plotNeuron);
subplot('position',[0.2300    0.1285    0.7050    0.7465])
surf(xm,ym,ysurf(D*sortNeuron-D+1,:),'EdgeColor','none')
axis([0 tend 1 plotNeuron(end)])
set(gca,'YTick',[])
set(gca,'XAxisLocation','top')
xlabel('Time (h)')
view(2)
subplot('position',[0.1300    0.1285    0.05    0.7465])
lregsc = lreg/length(region);
fill([0 0 1 1],[0 lregsc(1) lregsc(1) 0],[0.1 0.1 0.7]);
hold on
fill([0 0 1 1],lregsc(1)+[0 lregsc(2) lregsc(2) 0],[0.7 0.1 0.1]);
fill([0 0 1 1],sum(lregsc(1:2))+[0 lregsc(3) lregsc(3) 0],[0.2 0.8 0.2]);
fill([0 0 1 1],sum(lregsc(1:3))+[0 lregsc(4) lregsc(4) 0],[0.5 0.5 0.9]);
fill([0 0 1 1],sum(lregsc(1:4))+[0 lregsc(5) lregsc(5) 0],[0.9 0.5 0.5]);
fill([0 0 1 1],sum(lregsc(1:5))+[0 lregsc(6) lregsc(6) 0],[0.2 0.8 0.2]);
set(gca,'XTick',[],'YTick',[])
axis([0 1 0 1]);
view(2)
ylabel('Neuron label')
if L0>0,
    subplot('position',[0.2300    0.0685    0.7050    0.02])
    LDplot([0 tend],0,light,0,1,'scaled')
    set(gca,'XTick',[],'YTick',[])
    axis tight
end

figure(11); clf;
% set(gcf,'WindowStyle','docked')
xsurf=linspace(0,tend,600);
ysurf=deval(xsurf,sol);
[sortSc, indSc]=sort(sc);
[xm, ym]=meshgrid(xsurf,1:N);
subplot('position',[0.1300    0.2185    0.7750    0.7065])
surf(ym',xm',ysurf(D*(indSc-1)+1,:)','EdgeColor','none')
axis([1 N 0 tend])
set(gca,'XTick',[])
ylabel('Time (h)')
view(2)
subplot('position',[0.1300    0.1100    0.7750    0.05])
plot(1./sortSc);
set(gca,'YTick',[])
xlabel('Period scaling index')

%%% Global stats

%%% Order parameter

b=yint(1:D:end,:);     % to compute the order param on the whole time serie
M=mean(b,1);
num=mean(M.^2)-mean(M).^2;
B=mean(b.^2,2)-mean(b,2).^2;
den=mean(B);
R=num/den;
fprintf('Order parameter: R=%g \n',R);
bl=yint(D*lscn-D+1,:);
Ml=mean(bl,1);
numl=mean(Ml.^2)-mean(Ml).^2;
Bl=mean(bl.^2,2)-mean(bl,2).^2;
denl=mean(Bl);
Rl=numl/denl;
br=yint(D*rscn-D+1,:);
Mr=mean(br,1);
numr=mean(Mr.^2)-mean(Mr).^2;
Br=mean(br.^2,2)-mean(br,2).^2;
denr=mean(Br);
Rr=numr/denr;
fprintf('Order parameter in the left SCN: Rleft=%g \n',Rl);
fprintf('Order parameter in the right SCN: Rright=%g \n',Rr);

%%% Period

mfper=periode(xint,mf);
fprintf('Mean field period: P=%g \n',mfper);

