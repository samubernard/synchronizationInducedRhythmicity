% run_BWntStoch.m

N=10;
tend=240;
ic=repmat([0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35 1 1]',N,1);

spread=0.05;
sc = 1./(1+spread*randn(N,1)); 

kappa=0.6;
K=kappa/N*ones(1,N);
noise=[];
% noise(11:12:12*N)=1;
noise(12:12:12*N)=1;
noise=0.05*noise';

[t,y] = sdeEuler(@BWntStoch,[0 48],ic,noise,sc,K);

[t,y] = sdeEuler(@BWntStoch,[0 tend],y(end,:)',noise,sc,K);

%% ANALYSIS
b=y(:,1:12:end)';     % to compute the order param on the whole time serie
M=mean(b,1);
num=mean(M.^2)-mean(M).^2;
B=mean(b.^2,2)-mean(b,2).^2;
den=mean(B);
R=num/den;
fprintf('Order parameter: R=%g \n',R);

% variability of the mean field 
f=mean(y(:,1:12:end)');
[per, perstd, pmaxf]=periode(t,f);

% variability of individual coupled oscillators
perindn=[];
for i=1:N,
    fn=y(:,12*(i-1)+1)';
    [pern(i), perstdn(i), pmaxfind]=periode(t,fn);
    perindn=[perindn diff(t(pmaxfind))];
    plot([t(pmaxfind(1:10))' t(pmaxfind(1:10))']',repmat([i-0.45; i+0.45],1,10),'b','LineW',2)
end

% variability of individual uncoupled oscillators
icu=[0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35 1 1]';
[tu,yu] = sdeEuler(@BWntStoch,[0 48],icu,noise(1:12,:),1,kappa);
[tu,yu] = sdeEuler(@BWntStoch,[0 tend],yu(end,:)',noise(1:12,:),1,kappa);
fu=yu(:,1:12:end)';
[peru, perstdu, pmaxfu]=periode(tu,fu);

% save noise_effet_sims
