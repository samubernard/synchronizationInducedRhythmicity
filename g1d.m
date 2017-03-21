function y = g1d(t,x,xtau,sc,kappa,opt)
% G1DM Goodwin 1D with delay 

f=1./(1+xtau.^4);
m=x;
if opt==0,
    k=kappa*mean(x);
    L=0;
else
    k=kappa*mean(x)*((t<=192 | t>=384)*(t<576));
    L=2*kappa*(t>=576);
end

y = sc.*(f - m) + k + L;