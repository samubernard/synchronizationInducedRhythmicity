% goG1DM

N=10;
ic=rand(N,1);
sc=1+0.05*randn(N,1);

tspan=[0 3820];

sol=dde23(@g1d,7,ic,tspan,[],sc,.1,1);
figure(1); clf;
plot(sol.x,sol.y(1:10,:))
hold on
plot([192 192],[0 2.5])
plot([384 384],[0 2.5])
plot([576 576],[0 2.5])
axis([tspan(1) tspan(2) .3 2.5])
xlabel('Time (h)')
ylabel('concentration')
title('Goodwin 1D with delay n=100')
set(gca,'XTick',tspan(0):48:tspan(2))