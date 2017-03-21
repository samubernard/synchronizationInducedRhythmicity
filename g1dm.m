function y = g1dm(t,x,xtau,sc,kappa,opt)
% G1DM Goodwin 1D with delay and MM degradation

f=1./(1+xtau.^1.5);
m=x./(1+x);
if opt==0,
    k=kappa*mean(x);
    L=0;
else
    k=kappa*mean(x)*((t<=192 | t>=384)*(t<576));
    L=kappa*(t>=576);
end

y = sc.*(f - m) + k + L;