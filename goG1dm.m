% goG1DM

N=2;
ic=rand(N,1);
sc=1+0.05*randn(N,1);

sol=dde23(@g1dm,7,ic,[0 768],[],[0.95 1.05]',.1,1);
figure(1); clf;
plot(sol.x,sol.y(1:2,:))
hold on
plot([192 192],[0 2.5])
plot([384 384],[0 2.5])
plot([576 576],[0 2.5])
axis([0 768 .3 2.5])
xlabel('Time (h)')
ylabel('concentration')
title('Goodwin 1D with delay n=100')
set(gca,'XTick',0:48:768)